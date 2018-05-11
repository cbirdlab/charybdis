#!/bin/bash

#SBATCH --job-name=ClusterOTU
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load bio-misc 
module load gsl
module load obitools

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64

PREFIX=$1
BP=$2
INDIR=$3
OUTDIR=$4
CHUNKS=$5
GCL_BIN=$6

$GCL_BIN""/ClusterOTU.sh \
	$PREFIX $BP $INDIR $OUTDIR $CHUNKS $GCL_BIN
