#!/bin/bash

#SBATCH --job-name=AddEcotag
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load R

PREFIX=$1
INDIR=$2
OUTDIR=$3
DB_NAME=$4 
GCL_BIN=$5


Rscript $GCL_BIN""/OTU_CVT_addEcotag.R \
    $OUTDIR""/$PREFIX"".OTUvsTubes.ecotag.SP.csv \
    $OUTDIR""/$PREFIX"".OTU.ecotag.csv \
    $DB_NAME \
    $OUTDIR""/$PREFIX"".OTUvsTubes.ecotag.SPA.csv 





