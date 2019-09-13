# Generates charon from file Ecotag results, 
# where the assigned fasta has been converted to 
# CSV format using 'obitab'

# Arguments
args <- commandArgs (TRUE)

FILE_ECOTAG <- args[1]
FILE_SAMPLE <- args[2]
FILE_OUTPUT <- args[3]

# Global R options
options(stringsAsFactors = FALSE)

# Load data
data <- read.table(file = FILE_ECOTAG, row.names = NULL,
	header = TRUE, sep = ",", stringsAsFactors = FALSE)
samples <- read.table (file = FILE_SAMPLE, 
	header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "SAMPLE_ID", "COUNT")

# Filter unassigned sequences
data <- data[!is.na(data$best_match),]

# Extract data needed from Ecotag results
OTU_SEQID <- data$id
ECOTAG_SEQID <- data$best_match
SCINAME <- data$scientific_name

# Parse Ecotag header
sepseqs <- lapply(ECOTAG_SEQID,
	function (x) as.list(strsplit(x, "\\|" )[[1]]))
charon <- data.frame ("ACCESSION", "GI", "TAXID")
names(charon) <- c("ACCESSION", "GI", "TAXID")
charon = charon[FALSE, ]
for (seq in sepseqs){ 
	print(seq)
	charon[nrow(charon) + 1, ] = list(seq[[1]], seq[[2]], seq[[3]])
	print("---")
}

# Add additional columns to charon
charon$SCINAME = SCINAME 
charon$OTU_SEQID = OTU_SEQID
charon$ECOTAG_SEQID = ECOTAG_SEQID
# Add OTU counts to charon
charon <- merge (x = charon, y = samples, by = "OTU_SEQID")

# Order the columns
charon <- charon[, c ("OTU_SEQID", "ECOTAG_SEQID", "SCINAME", "TAXID",
	"SAMPLE_ID", "COUNT")]

write.table (x = charon, file = FILE_OUTPUT, 
	row.names = FALSE, col.names = FALSE, sep = ",", quote= FALSE)

