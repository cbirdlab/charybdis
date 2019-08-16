# Modified the column names for samples to also include the descriptive sample name, 
# rather than just the unique identification value


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CANNOT have column name starting with 'X' after the tube names

# Format of PRED_FILE:


args <- commandArgs(TRUE)

stopifnot(length(args) == 4)

CVT_FILE  <- args[1]
PRED_FILE <- args[2]
OUT_FILE  <- args[3]
SAMPLES_FILE <- args[4]

CVT  <- read.csv (file = CVT_FILE,  header = TRUE, stringsAsFactors = FALSE)
PRED <- read.csv (file = PRED_FILE, header = TRUE, stringsAsFactors = FALSE)
SAMPLES <- read.csv (file = SAMPLES_FILE, header = FALSE, stringsAsFactors = FALSE)

colnames <- colnames(CVT)
numSamples = nrow(SAMPLES)
startPos <- which(colnames == "Sample") + 1
if (identical(startPos, numeric(0))) {
	startPos <- which(colnames == "sample") + 1
}

stopPos  <- startPos + numSamples - 1
names <- colnames[startPos:stopPos]
PRED$Sample <- make.names(PRED$Sample)

predNames <- list()
newCols <- list()
for (t in names){
        prow = PRED[PRED$Sample == t, ]
        predNames <- c(predNames, prow$Description)
	newCols <- c(newCols, paste(prow$Description, t, sep="."))
}
predNames <- unlist(predNames)
newCols <- unlist(newCols)

newColsFMT <- colnames(CVT)
newColsFMT[startPos:stopPos] <- newCols

CVT_PRED <- CVT
colnames(CVT_PRED) <- newColsFMT

write.csv (x = CVT_PRED, quote = FALSE, file = OUT_FILE)
