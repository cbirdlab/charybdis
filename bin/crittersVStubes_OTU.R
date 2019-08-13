####################
# Input Parameters #
####################

args <- commandArgs (TRUE)

# Ensure 7 parameters
stopifnot (length (args) == 7)

# Name of input file, minus the extension
# Will be used as the prefix for all output data in order to 
# associate each step of the pipeline with the input file
# PREFIX <- "jones.COI"
PREFIX <- args[1]

# Samples File
# SAMPLES_FILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/in/jones.16S.samples.txt"
SAMPLES_FILE <- args[2]

# Total Counts File
# TOTAL_COUNTS_FILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.total_counts.csv"
TOTAL_COUNTS_FILE <- args[3]

# Charon file from which scientific names may be pulled
# CHARON_FILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.OTU.charon"
CHARON_FILE <- args[4]

# Samples file for attaching OTU information. 
# FILE_SAMPLE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.full.nonchimeras.clean.OTU.cluster.size_fix"
FILE_SAMPLE <- args[5]

# Output file name
# OUTFILE <- "/media/Wapuilani/evan/Charybdis_Runs/MaryJones/out/jones.16S.OTU.CVT"
OUTFILE <- args[6]

# Directory containing the NCBI flatfile taxonomic database
# TAXDIR <- "/media/Wapuilani/Databases/NCBI_TAXO/"
TAXDIR <- args[7]

#############
# R Options #
#############

options(stringsAsFactors = FALSE)

###############
# Dependecies #
###############

# Load libraries
suppressMessages (library ("CHNOSZ"))  # For NCBI database query
suppressMessages (library ("pracma"))  # For string manipulation functions
suppressMessages (library ("furrr"))  # For string manipulation functions

# Load taxonomic database into memory for faster access
# This will be used to get higher level taxonomic information for 
# a given scientific name
ncbi_nodes <- getnodes (taxdir=TAXDIR)
ncbi_names <- getnames (taxdir=TAXDIR)


####################
# Define Functions #
####################

# Function: expandTaxonRankName
# Purpose: Takes in a list of basic taxonomic rank names
#   and returns a list of those same names, but each also has
#   the "infra", "sub" and "super" prefixes
expandTaxonRankName<- function (x){
  c (strcat ( c ("infra", x)), strcat ( c ("sub", x)), x, strcat ("super", x) )
}


####################
# Read, Parse Data #
####################

# Obtain a list of taxonomic ranks in the form of "infraspecies, subspecies, species, superspecies, ......., superkingdom"
taxonomicRanksOfInterestCompressed <- c ("species", "genus", "family", "order", "class", "phylum", "kingdom")
#taxonomicRanksOfInterest <- sapply (X = taxonomicRanksOfInterestCompressed, function (x) expandTaxonRankName (x))
#taxonomicRanksOfInterest <- as.vector (taxonomicRanksOfInterest)
taxonomicRanksOfInterest <- taxonomicRanksOfInterestCompressed

samples <- unlist (read.table (SAMPLES_FILE, header = FALSE, stringsAsFactors = FALSE))
samples_names <- make.names (samples)
charon <- read.table (CHARON_FILE, sep = ",", stringsAsFactors = FALSE)
charon$V5 <- make.names (charon$V5)
seq_stats <- read.table (TOTAL_COUNTS_FILE, header = FALSE, sep = ":")


# Build table sequence IDs (query and target) and scientific name   
seqidQ_seqidT_sciname_taxid <- charon[,1:4]
colnames (seqidQ_seqidT_sciname_taxid) <- c ("SEQID_QUERY", "SEQID_TARGET", "SCINAME", "NCBI_TAXID")
# Compress into uniques
seqidQ_seqidT_sciname_taxid <- unique (seqidQ_seqidT_sciname_taxid[c ("SEQID_QUERY", "SEQID_TARGET", "SCINAME", "NCBI_TAXID")])

# Build empty table such that each column corrosponds to a sample
# and the number of rows is the numbers of rows in charon
sampleInit <- matrix (data = 0, nrow = nrow (seqidQ_seqidT_sciname_taxid), ncol = length (samples))
colnames (sampleInit) <- samples_names

# Combine into empty cVt
CVT_Expanded <- cbind.data.frame (seqidQ_seqidT_sciname_taxid, sampleInit)

# Assign (SeqQ, SeqT) counts
assign <- function (seqidQ, seqidT, taxid, sampleID, count) {
  CVT_Expanded[CVT_Expanded[,1]==seqidQ & CVT_Expanded[,2]==seqidT & CVT_Expanded[,4]==taxid, sampleID] <<- 
      CVT_Expanded[CVT_Expanded[,1]==seqidQ & CVT_Expanded[,2]==seqidT & CVT_Expanded[,4]==taxid, sampleID] + count
}

apply (X = charon, MARGIN = 1, function (x) {assign (x[[1]], x[[2]], as.numeric (x[[4]]), x[[5]], as.numeric (x[[6]]) ) })

cVT_Expanded_Backup <- CVT_Expanded

##############################
# Section for OTU attachment #
##############################
# Read the map of SEQID <-> OTU_SEQID
samples <- read.table (file = FILE_SAMPLE, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "COUNT", "SAMPLE_ID")
OTU_Seqid_map <- samples[ c ("OTU_SEQID", "QSEQID")]
# Add OTU_SEQID to each row
CVT_Expanded_OTU <- merge (y = CVT_Expanded, x = OTU_Seqid_map, by.y = "SEQID_QUERY", by.x = "QSEQID")
# Collapse by OTU 
CVT <- CVT_Expanded_OTU[, names (CVT_Expanded_OTU) != "SEQID_QUERY"]
CVT <- CVT[, names (CVT) != "QSEQID"]
CVT <- CVT[, names (CVT) != "SEQID_TARGET"]
CVT <- aggregate ( . ~ SCINAME + NCBI_TAXID + OTU_SEQID, data = CVT, FUN = sum )


######################
# End OTU attachment #
######################

# Now, need to compress the data such that the query sequence is not considered and the total counts for 
# each target is summed for each sample. 
##CVT <- CVT_Expanded[, names (CVT_Expanded) != "SEQID_QUERY"]
##CVT <- CVT[, names (CVT) != "SEQID_TARGET"]
##CVT <- aggregate ( . ~ SCINAME + NCBI_TAXID, data = CVT, FUN = sum )

critter_totals <- apply (X = CVT[, 4:ncol (CVT)], MARGIN = 1, sum)
CVT <- cbind.data.frame (CVT$OTU_SEQID, CVT$NCBI_TAXID, CVT$SCINAME, critter_totals, CVT[, 4:ncol (CVT)])
colnames (CVT)[1:4] <- c ("OTU_SEQID", "NCBI_TAXID", "SCINAME", "TOTAL")



################################################################
# Initialize and populate table of rank names for each critter #
################################################################

# For each unique scientific name, get associated taxonomic rank
ranks <- sapply (X = CVT$NCBI_TAXID, function (x) getrank (x, nodes=ncbi_nodes))

# Change "no rank" to "clade"
ranks <- sapply (X = ranks, function (x) { if (!is.na (x)){  if (x == "no rank") {"clade/other"} else {x}} else {"NA"} })
attributes (ranks)$names = NULL

# Initialize rank matrix
rankMat <- matrix (data = "", nrow = nrow (CVT), ncol = length (taxonomicRanksOfInterest))
# Use the list of taxonomic rank types ("infraspecies", .... "superkingdom") to label the columns
colnames (rankMat) <- taxonomicRanksOfInterest

# Combine the critter names, rank of critter name, empty rank matrix. and tube assingments into one table
full <- cbind.data.frame (CVT[, 1:2], rankMat, CVT[, 4:ncol (CVT)])
fullBackup <- full

# For each critter name, traverse the phylogenic tree (via NCBI database) to get higher rank classifications
# These are returned in form of numeric TAXID
# This is a list of lists, CEB updated to run in parallel using furrr
plan(multiprocess)
higherTaxa <- future_map(CVT$NCBI_TAXID, function (x) {  tryCatch ({
  allparents (id = x, taxdir = TAXDIR, nodes = ncbi_nodes)
}, warning = function (w){
  
}, error = function (e){
  c ("NA")
}, finally = {
}
)})

# For each TAXID returned, obtain the corrosponding rank
for (i in 1:(length (higherTaxa))){
  rankOfh <- sapply (X = higherTaxa[i], function (x) getrank (id = x, taxdir = TAXDIR, nodes = ncbi_nodes))
  # Use the "names" attribute of each list to associate the rank and TAXID
  attributes (higherTaxa[[i]])$names = rankOfh
}

# Use the taxonimic rank and the TAXID as coordinates to assign the scientific name
# in the appropriate field
for (i in 1:(length (higherTaxa))){
  for (j in 1:(length (higherTaxa[[i]]))){
    if (is.na (attributes (higherTaxa[[i]][j])$names)){ 
      break
    }
    
    if (attributes (higherTaxa[[i]][j])$names != "no rank"){  # Skip "no rank"
      colIdx <- attributes (higherTaxa[[i]][j])$names
      full[i,colIdx] <- as.character (sciname (id = as.numeric (higherTaxa[[i]][j]), taxdir = TAXDIR, names = ncbi_names))
    }
  }
}

# Reverse order in which taxonomic ranks are displayed in table
#fulltest <- full [c (1, 2, ((length (taxonomicRanksOfInterest) + 2):3), ((length (taxonomicRanksOfInterest) + 3) : ( length (samples_names) + (length (taxonomicRanksOfInterest) + 3))))]



######################################################
#### EXTRA: Use BOLD to fill in taxonomic "holes" ####
##### Because NCBI lacks some needed information #####
######################################################
# Less throroughly documented since it is a hacky solution 
# that is likely to be superceded

BOLD_API_STR = "http://www.boldsystems.org/index.php/API_Public/specimen?format=tsv&taxon="
suppressMessages (library ("bold"))

# Should always have the higher taxon for any given taxon
# for these normal categories
taxonomicRanksOfInterestCompressed

# Function to get the lowest available type of taxon 
getLowestTax <- function (taxNameRow){
  # the minimum index where the field is not an empty string
  idx <- min(which(taxNameRow != "" ))
  # the column name at that index  
  
  if (is.infinite (idx))
  {
    -1
  }
  else
  {
    idx
  } 
}


# Subset the full data set to get just the taxonomic info 
taxMat <- full[,3:(length (taxonomicRanksOfInterest) + 2)]
taxMat <- taxMat [ , taxonomicRanksOfInterestCompressed]

# Find the lowest available taxon in each row in order to determine which 
# higher level taxons are blank
lowest <- apply (X = taxMat, MARGIN = 1, function (x) getLowestTax (x))
lenMax <- length (taxonomicRanksOfInterestCompressed)

# Find the "holes" by checking for empty strings that occur
# after the lowest available taxon for a row
holes <- list () 
for (i in 1:nrow(taxMat)){
  if (lowest[[i]] > 0){
    holes[[i]] <- which (taxMat[i,][lowest[[i]]:lenMax] == "")    
  }
  else
  {
    holes[[i]] <- c (0)
  }
}


# Pull the taxonomic info for each row that has "holes"
for (i in 1:(nrow (taxMat)-1)){
  
  if (lowest [[i]] > 0 && length (holes[[i]]) > 0){
    scinameLow <- taxMat[i, lowest[[i]]]
    #print (scinameLow)
    
    res <- bold_tax_name (scinameLow)
    if ("taxon" %in% colnames (res)){
      boldTaxid = res[1,]$taxid
      
      res2 <- bold_tax_id (boldTaxid, includeTree = TRUE)
      #print (res2)
      selectRows <- taxonomicRanksOfInterestCompressed[holes[[i]]]
      
      want <- subset (res2, tax_rank %in% selectRows)
      
      if (nrow(want) > 0){
        for (k in 1:(nrow (want))){
          full[i, want[k,]$tax_rank] = want[k,]$taxon
        }        
      }
      else
      {
        # Could not find desired taxon 
      }
    }
    else
    {
      # Could not find desired taxon 
    }
  }
}

CVT <- full[with (full, order (full$TOTAL, decreasing = TRUE)),]


OUTFILE_STR <- OUTFILE
write.csv (CVT, OUTFILE_STR, row.names = FALSE)
