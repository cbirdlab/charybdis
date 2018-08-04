# Modified the column names for samples to also include the descriptive sample name, 
# rather than just the unique identification value


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!
# CANNOT have column name starting with 'X' after the tube names

# Format of PRED_FILE:


args <- commandArgs(TRUE)

stopifnot(length(args) == 3)

CVT_FILE  <- args[1]
PRED_FILE <- args[2]
OUT_FILE  <- args[3]

CVT  <- read.csv (file = CVT_FILE,  header = TRUE, stringsAsFactors = FALSE)
PRED <- read.csv (file = PRED_FILE, header = TRUE, stringsAsFactors = FALSE)

colnames <- colnames(CVT)
count = 0
lowest = 0
for(c in colnames){
        if (grepl("^X", c)){
                lowest = count
        }
        count = count + 1
}

startPos <- which(colnames == "sample") + 1
stopPos  <- lowest + 1

names <- colnames[startPos:stopPos]
PRED$Tube <- make.names(PRED$Tube)

predNames <- list()
newCols <- list()
for (t in names){
        prow = PRED[PRED$Tube == t, ]
        predNames <- c(predNames, prow$Description)
	newCols <- c(newCols, paste(prow$Description, t, sep=" - "))
}
predNames <- unlist(predNames)
newCols <- unlist(newCols)

newColsFMT <- colnames(CVT)
newColsFMT[startPos:stopPos] <- newCols

CVT_PRED <- CVT
colnames(CVT_PRED) <- newColsFMT

write.csv (x = CVT_PRED, quote = FALSE, file = OUT_FILE)
