rm(list=ls())
## Code to process Charybdis blast results
## The input for the "read_Charibdis" function is the csv output file from Charibdis.
## The second input "year" is to allow for grouping of multiple datasets by year they were run
## The example at the bottom assumes that column names are a composit of 3 things. THe predator species, the sample number, and the subsample ID. These are all seperated by a "." with an "X" preceeding the sample number.



library(tidyverse)
library(magrittr)


read_Charibdis<-function(file_name, year){
  #Function to read in standard output file from Charybdis blast results
  read_csv(file_name) %>%
    dplyr::select(-X1, -X, -suborder:-SEQUENCE, -TOTAL, -sample, SEQUENCE) %>%
    gather(ID, count, -OTU_SEQID:-kingdom, -SEQUENCE) %>%
    mutate(year = year) %>%
    dplyr::select(year, ID, OTU_SEQID, count, BLAST_IDENTITY_NCBI_BLAST, BLAST_QCOV_NCBI_BLAST, species:kingdom, everything()) %>%
    dplyr::rename(BLAST_IDENTITY = BLAST_IDENTITY_NCBI_BLAST,
                  BLAST_COVERAGE = BLAST_QCOV_NCBI_BLAST) %>%
    mutate(count = as.numeric(count)) %>%
    filter(count > 0)
}


filterIDs <- function(x, family_threshold = 97, order_threshold = 85){
  #Filter blast results by BLAST_Identity_score
  x %>%
    mutate(species = if_else(BLAST_IDENTITY < family_threshold, NA_character_, species),
           genus = if_else(BLAST_IDENTITY < family_threshold, NA_character_, genus),
           family = if_else(BLAST_IDENTITY < family_threshold, NA_character_, family),
           order = if_else(BLAST_IDENTITY < order_threshold, NA_character_, order),
           class = if_else(BLAST_IDENTITY < order_threshold, NA_character_, class),
           phylum = if_else(BLAST_IDENTITY < order_threshold, NA_character_, phylum),
           kingdom = if_else(BLAST_IDENTITY < order_threshold, NA_character_, kingdom))
}

## Example to read in 2 files and join together with year distinguishing the two files
blast_res <- read_Charibdis("FILE_1", 2016) %>%
  bind_rows(read_Charibdis("FILE_2", 2018)) %>%
  mutate(ID = str_replace_all(ID, '\\.\\.','\\.')) %>%
  separate(ID, into = c('Predator_Species','Sample','Subsample'), sep = '\\.', extra = 'merge', fill = 'right') %>%
  mutate(Sample = str_remove(Sample, "^X")) %>%
  filterIDs(family_threshold = 97, order_threshold = 85) %T>%
  write_csv('Full_metabarcoding_Blast_results.csv')

## Example taking the processed data and outputting a list of the unique species identified
unique_species <- blast_res %>%
  dplyr::select(species:kingdom, NCBI_TAXID) %>%
  filter(!is.na(species)) %>%
  distinct %T>%
  write_csv('Unique_Species_ID.csv')

