###########
# Options #
###########
args <- commandArgs (TRUE)
# Ensure >6 parameters 
# (only six needed, but 7 for older scripts that call this)
stopifnot (length (args) >= 6)

# Project prefix
PREFIX <- args[1]
# Samples File
SAMPLES_FILE <- args[2]
# Total Counts File
TOTAL_COUNTS_FILE <- args[3]
# Charon file from which scientific names may be pulled
CHARON_FILE <- args[4]
# Samples file for attaching OTU information. 
FILE_SAMPLE <- args[5]
# Output file name
OUTFILE <- args[6]

options(stringsAsFactors = FALSE)

# Load libraries
suppressMessages (library ("taxizedb")) # For NCBI database query
suppressMessages (library ("pracma", warn.conflicts = FALSE))  # For string manipulation functions
suppressMessages (library ("furrr", warn.conflicts = FALSE))  # For parallelization
suppressMessages (library ("tidyr", warn.conflicts = FALSE)) # for reshaping data

#############
# Functions #
#############

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
taxonomicRanksOfInterest <- taxonomicRanksOfInterestCompressed

samples <- unlist (read.table (SAMPLES_FILE, header = FALSE, stringsAsFactors = FALSE))
samples_names <- make.names (samples)
charon <- read.table (CHARON_FILE, sep = ",", stringsAsFactors = FALSE)
print (charon)
charon$V5 <- make.names (charon$V5)
colnames (charon) <- c ("QSEQID", "GenBankID", "SCINAME", "NCBI_TAXID", "SAMPLE_ID", "COUNT")
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

#########################
# Count OTUs per Sample #
#########################
# Read the map of SEQID <-> OTU_SEQID
samples <- read.table (file = FILE_SAMPLE, header = FALSE, sep = ',', stringsAsFactors = FALSE)
colnames (samples) <- c ("OTU_SEQID", "QSEQID", "SAMPLE_ID", "COUNT")
OTU_Seqid_map <- samples[ c ("OTU_SEQID", "QSEQID")]
# Merge charon and OTU_Seqid_map
charon_OTU_expanded <- merge (y = charon, x = OTU_Seqid_map, by.y = "QSEQID", by.x = "QSEQID")
charon_OTU_expanded <- charon_OTU_expanded[, names (charon_OTU_expanded) != "QSEQID"]
charon_OTU_expanded <- charon_OTU_expanded[, names (charon_OTU_expanded) != "GenBankID"]
charon_OTU <- aggregate( . ~ SCINAME + NCBI_TAXID + SAMPLE_ID + OTU_SEQID, data = charon_OTU_expanded, FUN = sum)

# Collapse by OTU 
CVT <- charon_OTU %>% spread(SAMPLE_ID, COUNT) 
CVT[is.na(CVT)] <- 0

# Now, need to compress the data such that the query sequence is not considered and the total counts for 
# each target is summed for each sample. 
critter_totals <- apply (X = CVT[, 4:ncol (CVT)], MARGIN = 1, sum)
CVT <- cbind.data.frame (CVT$OTU_SEQID, CVT$NCBI_TAXID, CVT$SCINAME, critter_totals, CVT[, 4:ncol (CVT)])
colnames (CVT)[1:4] <- c ("OTU_SEQID", "NCBI_TAXID", "SCINAME", "TOTAL")

#####################
# Get OTU taxonomy  #
#####################

# For each unique scientific name, get associated taxonomic rank
### OLD TAXO     ranks <- sapply (X = CVT$NCBI_TAXID, function (x) getrank (x, nodes=ncbi_nodes))
ranks <- sapply (X = CVT$NCBI_TAXID, function (x) taxid2rank (x))

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
plan(multicore)
higherTaxa <- future_map(CVT$NCBI_TAXID, function (x) {  tryCatch ({
  ##### OLD TAXO    allparents (id = x, taxdir = TAXDIR, nodes = ncbi_nodes)
  classification(x, db='ncbi')
}, warning = function (w){
  
}, error = function (e){
  c ("NA")
}, finally = {
}
)})


# Assign selected ranks, if available
# (Should be made dynamic...)
for (i in 1:(length (higherTaxa))){
    ht = higherTaxa[i][[1]][[1]]

    if (length(ht[ht$rank=="species", ]$name) > 0) {
        full[i, "species"] = ht[ht$rank=="species", ]$name
    }
    if (length(ht[ht$rank=="genus", ]$name) > 0) {
        full[i, "genus"] = ht[ht$rank=="genus", ]$name
    }
    if (length(ht[ht$rank=="family", ]$name) > 0) {
        full[i, "family"] = ht[ht$rank=="family", ]$name
    }
    if (length(ht[ht$rank=="order", ]$name) > 0) {
        full[i, "order"] = ht[ht$rank=="order", ]$name
    }
    if (length(ht[ht$rank=="class", ]$name) > 0) {
        full[i, "class"] = ht[ht$rank=="class", ]$name
    }
    if (length(ht[ht$rank=="phylum", ]$name) > 0) {
        full[i, "phylum"] = ht[ht$rank=="phylum", ]$name
    }
    if (length(ht[ht$rank=="kingdom", ]$name) > 0) {
        full[i, "kingdom"] = ht[ht$rank=="kingdom", ]$name
        }
}

full_ordered <- full[with (full, order (full$TOTAL, decreasing = TRUE)),]
write.csv (full_ordered, OUTFILE, row.names = FALSE)

################################################
# EXIT EARLY - make sure taxonomy changes work #
# Remaining script fills in missing taxonomy,  #
# but hoping that new NCBI data reliable       #
################################################ 

quit()



######################################################
#### EXTRA: Use BOLD to fill in taxonomic "holes" ####
##### Because NCBI lacks some needed information #####
######################################################
# Less throroughly documented since it is a hacky solution 
# that is likely to be superceded

BOLD_API_STR = "http://www.boldsystems.org/index.php/API_Public/specimen?format=tsv&taxon="
suppressMessages (library ("bold" , warn.conflicts = False))

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
    
    #res <- bold_tax_name (scinameLow)
    
    #replacing previous line with this one that should prevent crashing.  delete this comment if script runs
    res <- tryCatch( bold_tax_name(scinameLow), 
      error=function(cond) {
        message("error triggered by bold_tax_name()")
        message(paste(i,"th row with missing taxon was not fixed"))
        message(cond)
        message("")
        return("NA")
      },
      warning=function(cond) {
        message("warning triggered by bold_tax_name()")
        message(paste(i,"th row with missing taxon was not fixed"))
        message(cond)
        message("")
        return("NA")
      },
      finally={
      }
    )
    
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
