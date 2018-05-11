# This script builds an -extended- charon file from a BLAST-10-custom file and an OTU_Samples_Counts file
# It is extended beacause there is a final column of the OTU

# BLAST-10-custom : ' -outfmt "10 qseqid sseqid sscinames staxids evalue score bitscore" '

args <- commandArgs (TRUE)

# Ensure 4 parameters
stopifnot (length (args) == 4)

FILE_BLAST <- args[1]    # FILE_BLAST <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/BlastSansEnv/blast_genbankxNT.out"
FILE_SAMPLE <- args[2]   # FILE_SAMPLE <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/OTU_size/Simons.combo.nonchimeras.clean.OTU.samples_counts.txt"
FILE_OUTPUT <- args[3]   # FILE_OUTPUT <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/BlastSansEnv/Simons.combo.charon_OTU.csv"
TAXDIR <- args[4]        # TAXDIR <- "/media/Keaolani/ReferenceDatabase/TAXO/"


#############
# R Options #
#############

options(stringsAsFactors = FALSE)


###############
# Dependecies #
###############

suppressMessages (library("pracma"))
suppressMessages (library("CHNOSZ"))

# Load taxonomic database into memory for faster access
# This will be used to get higher level taxonomic information for 
# a given scientific name
ncbi_nodes <- getnodes (taxdir=TAXDIR)
ncbi_names <- getnames (taxdir=TAXDIR)

#############
# Load Data #
#############

data <- read.table (file = FILE_BLAST, header = FALSE, sep = ',', stringsAsFactors = FALSE)

colnames (data)[1:4] <- c ("OTU_SEQID", "SSEQID", "SSCINAMES", "STAXIDS")
blast <- data[,c ("OTU_SEQID", "SSEQID", "SSCINAMES", "STAXIDS")]

samples <- read.table (file = FILE_SAMPLE, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "SAMPLE_ID", "COUNT")

#######################
# Parse BLAST results #
#######################

# Only keep one TAXID
blast$STAXIDS <- sapply (X = blast$STAXIDS, FUN = function (x) gsub (pattern = ";.*", replacement = "", x = x))

# Filter top hit for each sequence
blast <- blast[!duplicated (blast$OTU_SEQID),]

########################
# Assign NCBI Taxon ID #
########################

# Get scientific name from TAXID
scinames <- sciname (id = blast$STAXIDS, taxdir = TAXDIR, names = ncbi_names)
scinames <- gsub (x = scinames, pattern = " ", replacement = "_")

# Replace SCINAME with new list
blast$SSCINAMES <- scinames

#################
# Parse Samples #
#################

charon <- merge (x = blast, y = samples, by = "OTU_SEQID")
print ( head (samples))
print ( nrow (blast))

print (blast)

print ( nrow (charon))

charon_Plus_OTU <- cbind.data.frame (charon$QSEQID, charon$SSEQID, charon$SSCINAMES, charon$STAXIDS, charon$SAMPLE_ID, charon$COUNT, charon$OTU_SEQID)
colnames (charon_Plus_OTU) <- c ("QUERY_SEQID", "REFERENCE_SEQID", "SCINAME", "NCBI_TAXID", "SAMPLE_ID", "COUNT", "OTU_SEQID")
charon <- charon_Plus_OTU
charon$OTU_SEQID <- NULL

write.table (x = charon, file = FILE_OUTPUT, row.names = FALSE, col.names = FALSE, sep = ",", quote= FALSE)
