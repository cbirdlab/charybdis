#!/bin/bash 

#SBATCH --job-name=AddVsearch
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

PREFIX=$1
INDIR=$2
OUTDIR=$3

echo "IDENTITY,ALI_LENGTH,NUM_MISMATCH,NUM_GAP_OPENS,START_POS_Q,END_POS_Q,START_POS_T,END_POS_T,EVALUE,BITSCORE" \
	> $OUTDIR""/$PREFIX"".OTU.vsearchStats.csv

## Bad: does not consider the order. 
#awk '{print $3","$4","$5","$6","$7","$8","$9","$10","$11","$12}' \
#	$OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
#	>> $OUTDIR""/$PREFIX"".OTU.vsearchStats.csv

awk -F',' '{print $1}' $OUTDIR""/$PREFIX"".OTUvsTubes.vsearch.SP.csv \
	| sed -e 's/"//g' \
	| while read line; do grep $line \
		$OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
		| awk '{print $3","$4","$5","$6","$7","$8","$9","$10","$11","$12}' \
		>> $OUTDIR""/$PREFIX"".OTU.vsearchStats.csv ; \
	done

paste -d',' $OUTDIR""/$PREFIX"".OTUvsTubes.vsearch.SP.csv \
	$OUTDIR""/$PREFIX"".OTU.vsearchStats.csv \
	> $OUTDIR""/$PREFIX"".OTUvsTubes.vsearch.SPA.csv 


