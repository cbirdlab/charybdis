#!/bin/bash

#SBATCH --job-name=OTUvsTube
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load parallel
module load R
module load obitools

PREFIX=$1
INDIR=$2
OUTDIR=$3
TAXON_DIR=$4
GCL_BIN=$5

Rscript $GCL_BIN""/crittersVStubes_OTU.R \
	$PREFIX \
	$INDIR/$PREFIX.samples.txt \
	$OUTDIR/$PREFIX.total_counts.csv \
	$OUTDIR/$PREFIX.OTU.charon \
	$OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX.OTUvsTubes.csv \
	$TAXON_DIR

# Add sequences
awk -F',' '{print $1}' $OUTDIR/$PREFIX.OTUvsTubes.csv | \
	sed 's/"//g' \
	> $OUTDIR/$PREFIX.OTU.queries.txt
echo "SEQUENCE" > $OUTDIR""/$PREFIX"".OTU.seqeuences.txt


cat $OUTDIR""/$PREFIX"".OTU.queries.txt \
	| while read line; do \
		grep -A 1 $line $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta ; \
	done \
	| grep -v '^>' \
	>> $OUTDIR""/$PREFIX"".OTU.seqeuences.txt


paste -d',' $OUTDIR""/$PREFIX"".OTUvsTubes.csv $OUTDIR""/$PREFIX"".OTU.seqeuences.txt \
	> $OUTDIR""/$PREFIX"".OTUvsTubes_S.csv 
