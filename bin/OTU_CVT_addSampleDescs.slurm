#!/bin/bash

#SBATCH --job-name=AddSampleDescriptions
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

PREFIX=$1
INDIR=$2
OUTDIR=$3
METHOD=$4

Rscript $GCL_BIN""/OTU_CVT_addSampleDescs.R \
	$OUTDIR""/$PREFIX.OTUvsTubes.$METHOD"".S.csv \
	$INDIR""/$PREFIX.sampledesc.csv \
	$OUTDIR""/$PREFIX.OTUvsTubes.$METHOD"".SP.csv \

