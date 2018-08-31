# This script builds and -extended- charon file from SAP output and an OTU_Samples_Counts file.
# It is extended because there is an OTU column (counts).
# Also, has additional SAP info: higher level ranks and their posterior probabilities.


args <- commandArgs (TRUE)

# Ensure  parameters
stopifnot (length (args) == 4)

# Get all assignments
FILE_SAP_ASSIGNMENTS <- args[1]
# Get posterior probabilities of assignments
FILE_SAP_PROBABILITIES <- args[2]
# Counts per OTU
FILE_SAMPLE <- args[3]
# Write Charon file to
FILE_OUTPUT <- args[4]

# R Options
options (stringsAsFactors = FALSE)

# Dependencies
library ("pracma")
library ("myTAI") # Sadly, script currently accesses interfor for sciname -> TAXID

# Load data
assignments <- read.table (file = FILE_SAP_ASSIGNMENTS, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (assignments) <- c ("FILE", "CUTOFF", "DETAIL", "OTU_SEQUENCE", 
			     "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", 
			     "HOMOLOGUE_COUNT", "MIN_FREQ_HOMOLOGUE", "MIN_TAXON_PROB")

probs <- read.table (file = FILE_SAP_PROBABILITIES, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (probs) <- c ("FILE", "ID", "RANK", "TAXON", "POSTERIOR_PROBABILITY")

samples <- read.table (file = FILE_SAMPLE, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "SAMPLE_ID", "COUNT")
samples$OTU_SEQID <- trimws (samples$OTU_SEQID)

# Functions
getLowestAssignment <- function (ranks){
	tail(ranks[!is.na(ranks)], n=1)
}

getTAXID <- function (sciname){

	tryCatch(taxonomy (organism = sciname, db = 'ncbi', output = 'taxid')[1,1], error = function(e) {""})
}
	

# Build Charon columns
data = cbind.data.frame (assignments$OTU_SEQUENCE) #, assignments$OTU_SEQUENCE)
colnames (data) <- c ("OTU_SEQID")
data$OTU_SEQID <- trimws (data$OTU_SEQID)
data$REFERENCE_SEQID <- sapply (1:nrow(data), function (x) { NA })

data$SCINAME <- apply (assignments[,c ("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")], 1, function(x) getLowestAssignment(x))
data$SCINAME <- gsub (" sp.*", "", data$SCINAME)
data$NCBI_TAXID <- sapply (data$SCINAME, function (x) getTAXID (x))

write.table (x = data, file = "TAXIDS.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)

data$SCINAME <- gsub (" ", "_", data$SCINAME)
# Merge with sample_counts file
charon <- merge (x = samples, y = data, by = "OTU_SEQID")
charon_Plus_OTU <- cbind.data.frame (charon$QSEQID, charon$REFERENCE_SEQID, charon$SCINAME, charon$NCBI_TAXID, charon$SAMPLE_ID, charon$COUNT, charon$OTU_SEQID)
# Cannot remember why, but eventually dropped the need for 'OTU_SEQID'..
charon <-charon_Plus_OTU
charon$OTU_SEQID <- NULL

# Write
write.table (x = charon, file = FILE_OUTPUT, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
