# Generates charon file from Vsearch's 'blast6' output. 
# Assumes that Vsearch database sequences headers are of format:
#         ACCESSION|GI|TAXID

# Arguments

args <- commandArgs (TRUE)

FILE_VSEARCH <- args[1]
FILE_SAMPLE <- args[2]
FILE_OUTPUT <- args[3]
TAXDIR <- args[4]

# Global R options

options(stringsAsFactors = FALSE)

# Dependencies
suppressMessages (library("pracma"))
suppressMessages (library("CHNOSZ"))

# Load taxonomic database into memory for faster access
ncbi_nodes <- getnodes (taxdir=TAXDIR)
ncbi_names <- getnames (taxdir=TAXDIR)

# Load data
data <- read.table (file = FILE_VSEARCH, 
	header = FALSE, sep = '\t', stringsAsFactors = FALSE)

samples <- read.table (file = FILE_SAMPLE, header = FALSE, sep = ',', 
	stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "SAMPLE_ID", "COUNT")

# Extract needed data from Vsearch results

OTU_SEQID <- data$V1      # OTU sequence header
VSEARCH_SEQID <- data$V2  # Vsearch db sequence header

# Parse Vsearch header, add to charon
sepseqs <- lapply(VSEARCH_SEQID, function (x) as.list(strsplit(x, "\\|")[[1]]))
charon <- data.frame ("ACCESSION", "GI", "TAXID")
names(charon) <- c("ACCESSION", "GI", "TAXID")
charon = charon[FALSE, ]
for (seq in sepseqs){
	charon[nrow(charon) + 1, ] = list(seq[[1]], seq[[2]], seq[[3]])
}

# Add scientific names to charon
scinames <- sciname (id = charon$TAXID, taxdir = TAXDIR, names = ncbi_names)
scinames <- gsub (x = scinames, pattern = " ", replacement = "_") 
charon$SCINAME <- scinames

# Add sequence headers to charon
charon$OTU_SEQID <- OTU_SEQID
charon$VSEARCH_SEQID <- VSEARCH_SEQID

# Add OTU counts to charon
charon <- merge (x = charon, y = samples, by = "OTU_SEQID")

charon <- charon[, c ("OTU_SEQID", "VSEARCH_SEQID", "SCINAME", "TAXID", 
                      "SAMPLE_ID", "COUNT")]

write.table (x = charon, file = FILE_OUTPUT,
	 row.names = FALSE, col.names = FALSE, sep = ",", quote= FALSE)





