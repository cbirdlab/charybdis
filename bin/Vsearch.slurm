#!/bin/bash

#SBATCH --job-name=VsearchAssignment
#SBATCH --time=96:00:00
#SBATCH -p normal 
#SBATCH --nodes=1

module load vsearch
module load R

PREFIX=$1
INDIR=$2
OUTDIR=$3
ID=$4
DB=$5
GCL_BIN=$6

TAXDIR="/work/hobi/GCL/db/TAXO/"

# VSEARCH
bash $GCL_BIN""/Vsearch.sh \
	$PREFIX $INDIR $OUTDIR $ID $DB

# Convert to Charon format
bash $GCL_BIN""/vsearch_getTAXIDfromBOLDseqid.sh $OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
	$TAXDIR | sed -e 's/ /_/g'  \
	> $OUTDIR/$PREFIX"".OTU.vsearchTAXIDS.csv

paste $OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
	$OUTDIR""/$PREFIX"".OTU.vsearchTAXIDS.csv  \
	> $OUTDIR""/$PREFIX"".OTU.vsearchB6out-fix.txt

Rscript $GCL_BIN""/vsearch_blast6custom_to_charon_OTU.R \
	$OUTDIR""/$PREFIX"".OTU.vsearchB6out-fix.txt \
	$OUTDIR""/$PREFIX"".full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR""/$PREFIX"".OTU.vsearch.charon \
	$OUTDIR""/$PREFIX"".trash

