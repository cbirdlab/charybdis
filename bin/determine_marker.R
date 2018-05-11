args = commandArgs(trailingOnly=TRUE)
file_BLAST <- args[1] # file_BLAST <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out2/sample1.blastOut.fix"
file_COI <- args[2]
file_ITS <- args[3]
file_16S <- args[4]

library("stringr")
stringsAsFactors = "FALSE"

blast <- read.csv (file = file_BLAST, stringsAsFactors = FALSE, header = FALSE)
colnames (blast) <- c ("QSEQID", "SSEQID", "SSCINAMES", "STAXIDS", "STITLE")

# Filter top hit for each sequence
blast <- blast[!duplicated (blast$QSEQID),]
blast$SSEQID <- str_split_fixed (blast$SSEQID, "\\|", 4)[,2]

# Only keep one TAXID
blast$STAXIDS <- sapply (X = blast$STAXIDS, FUN = function (x) gsub (pattern = ";.*", replacement = "", x = x))

# Categorize
COI_terms <- c ("coi", "cytochrome", "oxidase", "cox1")
ITS_terms <- c ("internal", "transcribed", "spacer", "ITS1", "ITS2", "ITS3")
sixteenS_terms <- c ("16S")

blast_COI <- blast[grepl (pattern = paste (COI_terms, collapse = '|'), 
		x = blast$STITLE, ignore.case = TRUE), ]
blast_ITS <- blast[grepl (pattern = paste (ITS_terms, collapse = '|'),
		x = blast$STITLE, ignore.case = TRUE), ]
blast_16S <- blast[grepl (pattern = paste (sixteenS_terms, collapse = '|'),
		x = blast$STITLE, ignore.case = TRUE), ]

write.table (x = blast_COI$QSEQID, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",",
	file = file_COI)
write.table (x = blast_ITS$QSEQID, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",",
	file = file_ITS)
write.table (x = blast_16S$QSEQID, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",",
	file = file_16S)
