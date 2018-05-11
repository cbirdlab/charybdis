#!/bin/bash 

#SBATCH --job-name=AddBlast
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

PREFIX=$1
INDIR=$2
OUTDIR=$3

echo "IDENTITY,EVALUE,SCORE,BITSCORE,QCOV,QCOVHSP,QCOVUS" \
	> $OUTDIR""/$PREFIX"".OTU.blastStats.csv

awk -F',' '{print $1}' $OUTDIR""/$PREFIX"".OTUvsTubes_S.csv \
	|  sed -e 's/"//g' \
	| while read line; do \
		grep -m 1 $line $OUTDIR""/$PREFIX"".OTU.blastOut \
		| awk -F',' '{print $8","$5","$6","$7","$9","$10","$11}' \
		>> $OUTDIR""/$PREFIX"".OTU.blastStats.csv  ; \
	done

paste -d',' $OUTDIR""/$PREFIX"".OTUvsTubes_S.csv \
	$OUTDIR""/$PREFIX"".OTU.blastStats.csv \
	> $OUTDIR""/$PREFIX"".OTUvsTubes_SB.csv 


