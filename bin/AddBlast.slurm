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

Rscript $GCL_BIN""/OTU_CVT_addBlast.R \
	$OUTDIR""/$PREFIX"".OTUvsTubes.blast.SP.csv \
	$OUTDIR""/$PREFIX"".OTU.blastOut \
	$DB_NAME \
	$OUTDIR""/$PREFIX"".OTUvsTubes.blast.SPA.csv
