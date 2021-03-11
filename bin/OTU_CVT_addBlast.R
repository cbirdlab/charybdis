# Append Blast scores to a CrittersVsTubes (CVT)

# Script requires 4 parameters

args <- commandArgs (TRUE)

# Ensure 4 parameters
stopifnot (length (args) == 5)

CVT_FILE <- args[1]     #CVT_FILE  <- "out/Drymon-Fall2018.OTUvsTubes.blast.SP.csv"
BLAST_FILE <- args[2]   # BLAST_FILE <- "out/Drymon-Fall2018.OTU.blastOut"
SOURCE_NAME <- args[3]  # SOURCE_NAME <- "NCBI"
OUTPUT_FILE <- args[4]  # OUTPUT_FILE <- "out/Drymon-Fall2018.OTUvsTubes.blast.SPA.csv"
COMMON_NAMES_FILE <- args[5] # COMMON_NAMES_FILE <- "common_names.txt"


CVT <- read.csv (file = CVT_FILE, header = TRUE, stringsAsFactors = FALSE)
COMMON_NAMES <- read.csv(COMMON_NAMES_FILE, header = FALSE, stringsAsFactors = FALSE)
names(COMMON_NAMES) <- "CommonName"
CVT<-data.frame(cbind(CVT[,1:3],COMMON_NAMES, CVT[,4:length(names(CVT))]))

BLAST <- read.csv (file = BLAST_FILE, header = FALSE, stringsAsFactors = FALSE)

blast_columns <- c ("OTU_SEQID", "SSEQID", "SSCINAMES", "STAXIDS", paste0("BLAST_E-SCORE_", SOURCE_NAME),
                    paste0("BLAST_SCORE_", SOURCE_NAME), paste0("BLAST_BIT_SCORE_", SOURCE_NAME),
                    paste0("BLAST_IDENTITY_", SOURCE_NAME), paste0("BLAST_QCOV_", SOURCE_NAME) )

colnames (BLAST)[1:9] <- blast_columns
BLAST <- BLAST[,blast_columns]

# Only keep one TAXID
BLAST$STAXIDS <- sapply (X = BLAST$STAXIDS, FUN = function (x) gsub (pattern = ";.*", replacement = "", x = x))

# Filter top hit for each sequence
BLAST <- BLAST[!duplicated (BLAST$OTU_SEQID),]

CVT_BLAST <- merge (x = BLAST, y = CVT, by.x = "OTU_SEQID", by.y = "OTU_SEQID")

CVT_BLAST$SSCINAMES = NULL
CVT_BLAST$STAXIDS = NULL
write.csv (x = CVT_BLAST, quote = FALSE, file = OUTPUT_FILE)
