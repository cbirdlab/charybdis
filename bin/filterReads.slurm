#!/bin/bash

#SBATCH --job-name=filterReads
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

if [ -v "$SLURM_JOB_ID" ]; then
	module load obitools
	module load parallel
	module load vsearch
	module load fastx_toolkit
fi

PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
MARKER_LEN_LOWER=$5
MARKER_LEN_HIGHER=$6
CHIMERA_DB=$7
GCL_BIN=$8

$GCL_BIN""/ObiToolsPipeline.sh \
	$PREFIX $INDIR $OUTDIR $CHUNKS \
	$MARKER_LEN_LOWER  $MARKER_LEN_HIGHER $CHIMERA_DB $GCL_BIN

