# Append Blast scores to a CrittersVsTubes (CVT)

# Script requires 4 parameters

args <- commandArgs (TRUE)

# Ensure 4 parameters
stopifnot (length (args) == 4)

CVT_FILE <- args[1]     # CVT_FILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.OTU.CVT"
BLAST_FILE <- args[2]   # BLAST_FILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.OTU.blastOut"
SOURCE_NAME <- args[3]  # SOURCE_NAME <- "NCBI-custom"
OUTPUT_FILE <- args[4]  # OUTPUT_FILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.OTU.CVT-EXT"

CVT <- read.csv (file = CVT_FILE, header = TRUE, stringsAsFactors = FALSE)
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

#clean up this file
CVT_BLAST$X <- NULL

col_names <- colnames(CVT_BLAST)
startcol <- which(col_names == "superorder")
endcol <- which(col_names == "species.group")
CVT_BLAST[,startcol:endcol] <- NULL

write.csv (x = CVT_BLAST, quote = FALSE, file = paste(OUTPUT_FILE,"cleaner.csv", sep="."))

                      
