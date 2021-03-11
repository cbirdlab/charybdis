#!/bin/bash 

#SBATCH --job-name=AddBlast
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

>&2 echo ""
>&2 echo `date` "Running AddBlast.slurm"

if [ -v "$SLURM_JOB_ID" ]; then
	module load R
	module load parallel
fi

PREFIX=$1
INDIR=$2
OUTDIR=$3
DB_NAME=$4
GCL_BIN=$5
TAXON_DIR=$6
THREADS=$7
METHOD=$8

#get common names from taxonomic db
>&2 echo `date` "     Extract common names from NCBI tax db"
OvT_FILE=$OUTDIR""/$PREFIX"".OTUvsTubes.$METHOD.SP.csv
get_common_name(){
	PATTERN=$1
	TAXON_DIR=$2
	RESULT=$(grep -P -m 1 "^${PATTERN}\t.*genbank common name" ${TAXON_DIR}/names.dmp | cut -f3)
	if [ "$RESULT" == "" ]; then
		RESULT=NA
	fi
	echo $RESULT
}
export -f get_common_name
awk -F "\"*,\"*" '{print $3}' $OvT_FILE  | tail -n +2 | parallel --no-notice -j $THREADS -k "get_common_name {} $TAXON_DIR" > common_names.txt

# Ensure no commas since CSV
sed -i -e 's/,/_/g' common_names.txt

>&2 echo `date` "     Run OTU_CVT_addBlast.R"

>&2 cat $GCL_BIN""/OTU_CVT_addBlast.R

Rscript $GCL_BIN""/OTU_CVT_addBlast.R \
	$OUTDIR""/$PREFIX"".OTUvsTubes.$METHOD.S.csv \
	$OUTDIR""/$PREFIX"".OTU.$METHOD.Out \
	$DB_NAME \
	$OUTDIR""/$PREFIX"".OTUvsTubes.$METHOD.SA.csv \
	common_names.txt

rm -f common_names.txt
