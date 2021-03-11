#!/bin/bash 

#SBATCH --job-name=AddVsearch
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load R

PREFIX=$1
INDIR=$2
OUTDIR=$3 
DB_NAME=$4
GCL_BIN=$5 

Rscript $GCL_BIN""/OTU_CVT_addVsearch.R \
    $OUTDIR""/$PREFIX"".OTUvsTubes.vsearch.S.csv \
    $OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
    $DB_NAME \
    $OUTDIR""/$PREFIX"".OTUvsTubes.vsearch.SA.csv

