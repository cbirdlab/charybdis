# Add Vsearch scores to OTUvsTubes

args <- commandArgs (TRUE)

stopifnot (length (args) == 4)

CVT_FILE <- args[1]
VSEARCH_FILE <- args[2]  # B6 output
SOURCE_NAME <- args[3] # EX: NCBI-NT
OUTPUT_FILE <- args[4]

CVT <- read.csv (file = CVT_FILE, header = TRUE, stringsAsFactors = FALSE)
VSEARCH <- read.csv (file = VSEARCH_FILE, header = FALSE, stringsAsFactors = FALSE, sep = '\t')

vsearch_columns <- c("OTU_SEQID", "SSEQID", 
                     paste0("VSEARCH_IDENTITY_", SOURCE_NAME), 
                     paste0("VSEARCH_ALI_LENGTH_", SOURCE_NAME), 
                     paste0("VSEARCH_NUM_MISMATCH_", SOURCE_NAME), 
                     paste0("VSEARCH_NUM_GAP_OPENS_", SOURCE_NAME), 
                     paste0("VSEARCH_START_POS_Q_", SOURCE_NAME), 
                     paste0("VSEARCH_END_POS_Q_", SOURCE_NAME), 
                     paste0("VSEARCH_START_POS_T_", SOURCE_NAME), 
                     paste0("VSEARCH_END_POS_T_", SOURCE_NAME), 
                     paste0("VSEARCH_EVALUE_", SOURCE_NAME), 
                     paste0("VSEARCH_BITSCORE_", SOURCE_NAME)) 

colnames (VSEARCH)[1:12] <- vsearch_columns
VSEARCH <- VSEARCH[,vsearch_columns]


CVT_VSEARCH <- merge (x = VSEARCH, y = CVT, by.x = "OTU_SEQID", by.y = "OTU_SEQID") 

write.csv (x = CVT_VSEARCH, quote = FALSE, file = OUTPUT_FILE)

