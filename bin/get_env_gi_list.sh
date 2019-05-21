## When using an entrex query to pull data from NCBI,
## there is often a default maximum size of 20.
## Thus, you will need to run it once to find the count so that
## you can set the default size manually to pull all records

# Get default number of records by using the query online, at the following URL
### https://www.ncbi.nlm.nih.gov/nuccore/?term=%22environmental+samples%22%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+sp%5BTitle%5D
# Find the total number of records
# Look for "Items: 1 to 20 of <TOTAL_COUNT>"


## Get GI list from NCBI using wget with an entrez query
## Replace <TOTAL_COUNT> with that from the website
TOTAL_COUNT=$1

wget 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="environmental samples"[organism]%20or%20metagenomes[orgn]%20or%20sp[Title]&retmax='"$TOTAL_COUNT"  -O env.NCBI_nucl.temp

## Format as simply a list of GIs

grep "<Id>" env.NCBI_nucl.temp | \
grep -o -e "[0-9]*" > env.NCBI_nucl.gi

## Remove temp file

rm env.NCBI_nucl.temp
