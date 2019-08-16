#!/bin/bash 

#SBATCH --job-name=AddBlast
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

if [ -v "$SLURM_JOB_ID" ]; then
	module load R
fi

PREFIX=$1
INDIR=$2
OUTDIR=$3
DB_NAME=$4
GCL_BIN=$5
TAXON_DIR=$6
THREADS=$7

#get common names from taxonomic db
OvT_FILE=$OUTDIR""/$PREFIX"".OTUvsTubes.blast.SP.csv
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

Rscript $GCL_BIN""/OTU_CVT_addBlast.R \
	$OUTDIR""/$PREFIX"".OTUvsTubes.blast.SP.csv \
	$OUTDIR""/$PREFIX"".OTU.blastOut \
	$DB_NAME \
	$OUTDIR""/$PREFIX"".OTUvsTubes.blast.SPA.csv \
	common_names.txt
