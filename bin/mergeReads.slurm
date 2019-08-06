#!/bin/bash

#SBATCH --job-name=mergeReads
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

if [ -v "$SLURM_JOB_ID" ]; then
  module load obitools
  module load parallel
fi

echo; echo `date` Running mergeReads.* ...

PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
GCL_BIN=$5

$GCL_BIN""/mergeReads.sh $PREFIX $INDIR $OUTDIR $CHUNKS $GCL_BIN
