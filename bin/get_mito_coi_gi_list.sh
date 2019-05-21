#### How to create 'mitochondial_coi.NCBI_nucl.gi'

## When using an entrex query to pull data from NCBI,
## there is often a default maximum size of 20.
## Thus, you will need to run it once to find the count so that
## you can set the default size manually to pull all records

# Get default number of records by using the query online, at the following URL
### 		https://www.ncbi.nlm.nih.gov/nuccore/?term=mitochondria%20or%20cytochrome%20or%20coi%20or%20co1%20or%20cox1%20or%20coxi%20or%20%22mitochondrial%20genome%22%20or%20%22mitochondria%20genome%22
# Find the total number of records
# Look for "Items: 1 to 20 of <TOTAL_COUNT>"


## Get GI list from NCBI using wget with an entrez query
## Replace <TOTAL_COUNT> with that from the website
TOTAL_COUNT=$1

wget 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=mitochondria or cytochrome or coi or co1 or cox1 or coxi or "mitochondrial genome" or "mitochondria genomea"&retmax='"$TOTAL_COUNT"  -O mitochondrial_coi.NCBI_nucl.temp

## Format as simply a list of GIs

grep "<Id>" mitochondrial_coi.NCBI_nucl.temp | \
grep -o -e "[0-9]*" > mitochondrial_coi.NCBI_nucl.gi

## Remove temp file

rm mitochondrial_coi.NCBI_nucl.temp
