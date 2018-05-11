args <- commandArgs(trailingOnly = TRUE)
file_BLAST_GIs <- args[1]
file_markers <- args[2]
file_COI <- args[3]
file_ITS <- args[4]
file_16S <- args[5]

stringsAsFactors = "FALSE"

blast <- read.csv (file = file_BLAST_GIs, header = FALSE, stringsAsFactors = FALSE)
colnames (blast) <- c ("QSEQID", "GI")

markers <- read.csv (file = file_markers, header = FALSE, stringsAsFactors = FALSE)
colnames (markers) <- c ("GI", "MARKER")
markers <- markers

data <- merge (x = blast, y = markers, by = "GI")

coi <- data[data$MARKER == 1, ]
its <- data[data$MARKER == 2, ]
sixteenS <- data[data$MARKER == 3, ]
other <- data[data$MARKER == 0, ]

coi_seqids <- coi$QSEQID
its_seqids <- its$QSEQID
sixteenS_seqids <- sixteenS$QSEQID


write.table (x = coi_seqids, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",",
  file = file_COI)

write.table (x = its_seqids, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",",
  file = file_ITS)

write.table (x = sixteenS_seqids, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",",
  file = file_16S)

