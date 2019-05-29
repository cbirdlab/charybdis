#!/bin/bash

#SBATCH --job-name=VsearchAssignment
#SBATCH --time=96:00:00
#SBATCH -p normal 
#SBATCH --nodes=1

PREFIX=$1
INDIR=$2
OUTDIR=$3
ID=$4
DB=$5
GCL_BIN=$6
TAXDIR=$7

module load vsearch

# VSEARCH
bash $GCL_BIN""/Vsearch.sh \
	$PREFIX $INDIR $OUTDIR $ID $DB

module unload vsearch
module load R

# Convert to Charon format
Rscript $GCL_BIN""/vsearch_B6out_to_charon.R \
    $OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
    $OUTDIR""/$PREFIX"".full.nonchimeras.clean.OTU.cluster.size_fix \
    $OUTDIR""/$PREFIX"".OTU.vsearch.charon \
    $TAXDIR
