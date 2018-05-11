# This script builds an -extended- charon file from a VSearch output formatted with the "blast6" style. 
# Requires that vsearch output as well an OTU_Samples_Counts file.
# It is considered extended because there is a final column containing the OTU

args <- commandArgs (TRUE)

# Ensure 4 parameters
stopifnot (length (args) == 4)

FILE_VSEARCH <- args[1]            # FILE_VSEARCH <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/PRES_vsearch/OTU.vsearched.blast6out.corrected.csv"
FILE_SAMPLE <- args[2]             # FILE_SAMPLE <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/OTU_size/FIX/Simons.combo.nonchimeras.clean.OTU.samples_counts.BEST.csv"
FILE_OUTPUT_CHARON <- args[3]      # FILE_OUTPUT_CHARON <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/Simons.combo.vsearch.charon.csv"
FILE_OUTPUT_CHARON_OTU <- args[4]  # FILE_OUTPUT_CHARON_OTU  <- "/media/Wapuilani/evan/Charybdis_Runs/Combined/Simons.combo.vsearch.charon_OTU.csv"

TAXDIR <- "/work/hobi/GCL/db/TAXO"

#############
# R Options #
#############

options(stringsAsFactors = FALSE)


###############
# Dependecies #
###############

suppressMessages (library("stringr"))

#############
# Load Data #
#############

data <- read.table (file = FILE_VSEARCH, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames (data)[1:2] <- c ("OTU_SEQID", "TSEQID")
colnames (data)[13] <- "REF"
vsearch <- data[,c ("OTU_SEQID", "TSEQID", "REF")]

samples <- read.table (file = FILE_SAMPLE, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "SAMPLE_ID", "COUNT")

#########################
# Parse Vsearch results #
#########################

# Extract Scinames and TAXIDs from Target SEQID
#seqid <- vsearch$TSEQID
#seqid_split <- str_split_fixed (seqid, "\\|", 3)
#scinames <- seqid_split[,2]
#taxids <- seqid_split[,3]
#taxids <- gsub (x = taxids, pattern = "\\.[0-9]*", replacement = "")

seqid <- vsearch$REF
seqid_split <- str_split_fixed (seqid, ",", 2)
scinames <- seqid_split[,1]
taxids <- seqid_split[,2]

vsearch$REF <- NULL

# VSEARCH Error correction: If taxid is not numeric, replace with 0
taxids[which(sapply (X = sapply (X = taxids, as.numeric), is.na))] <- 0

# Add Sciname and TAXID as columns in 
vsearch <- cbind.data.frame (vsearch, scinames, taxids)

# Update column names
colnames (vsearch)[1:4] <- c ("OTU_SEQID", "TSEQID", "SCINAME", "TAXID")

######################
# Merge by OTU_SEQID #
######################

charon <- merge (x = samples, y = vsearch, by = "OTU_SEQID")
charon_Plus_OTU <- cbind.data.frame (charon$QSEQID, charon$TSEQID, charon$SCINAME, charon$TAXID, charon$SAMPLE_ID, charon$COUNT, charon$OTU_SEQID)
colnames (charon_Plus_OTU) <- c ("QUERY_SEQID", "REFERENCE_SEQID", "SCINAME", "NCBI_TAXID", "SAMPLE_ID", "COUNT", "OTU_SEQID")
charon <- charon_Plus_OTU
charon$OTU_SEQID <- NULL

# Write Charon
write.table (x = charon, file = FILE_OUTPUT_CHARON, row.names = FALSE, col.names = FALSE, sep = ",", quote= FALSE)

# Write Charon_OTU
write.table (x = charon_Plus_OTU, file = FILE_OUTPUT_CHARON_OTU, row.names = FALSE, col.names = FALSE, sep = ",", quote= FALSE)

