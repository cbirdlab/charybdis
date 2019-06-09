# Append Ecotag scores to a CrittersVsTubes (CVT)

# Script requires 4 parameters

args <- commandArgs (TRUE)

stopifnot (length (args) == 4)

CVT_FILE <- args[1]
BLAST_FILE <- args[2]
SOURCE_NAME <- args[3]
OUTPUT_FILE <- args[4]

CVT <- read.csv(file = CVT_FILE, header = TRUE,
	stringsAsFactors = FALSE)
ECOTAG <- read.csv (file = BLAST_FILE, header = TRUE, row.names = NULL,
	stringsAsFactors = FALSE)

ECOTAG <- ECOTAG[, c ("id", "best_identity")]
colnames (ECOTAG) <- c ("OTU_SEQID", 
	paste0("BLAST_IDENTITY_", SOURCE_NAME))


CVT_ECOTAG <- merge (x = ECOTAG, y = CVT, by = "OTU_SEQID")

write.csv (x = CVT_ECOTAG, quote = FALSE, file = OUTPUT_FILE)

