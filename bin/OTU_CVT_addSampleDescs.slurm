#!/bin/bash

#SBATCH --job-name=AddSampleDescriptions
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

if [ -v "$SLURM_JOB_ID" ]; then
	module load R
fi

PREFIX=$1
INDIR=$2
OUTDIR=$3
GCL_BIN=$4
METHOD=$5

Rscript $GCL_BIN""/OTU_CVT_addSampleDescs.R \
	$OUTDIR""/$PREFIX.OTUvsTubes.$METHOD"".S.csv \
	$INDIR""/$PREFIX.sampledescs.csv \
	$OUTDIR""/$PREFIX.OTUvsTubes.$METHOD"".SP.csv \
	$INDIR""/$PREFIX.samples.txt

