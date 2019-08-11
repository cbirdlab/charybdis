#!/bin/bash

#SBATCH --job-name=BlastMeta
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

if [ -v "$SLURM_JOB_ID" ]; then
	module load blast+
	module load parallel
	module load bio-misc
	module load genometools
	module load R
fi

PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
THRESHOLD=$5
BLAST_DB=$6
GCL_BIN=$7
TAXON_DIR=$8
BLAST_IGNORE=$9

# BLAST
$GCL_BIN""/BlastMeta.sh \
	$PREFIX $INDIR $OUTDIR $CHUNKS \
	$THRESHOLD $BLAST_DB $BLAST_IGNORE

# Combine into single FASTA
rm $OUTDIR""/$PREFIX.OTU.blastOut
cat $OUTDIR""/$PREFIX.OTU.blastOut.* \
	> $OUTDIR/$PREFIX.OTU.blastOut

# Convert to Charon format
Rscript $GCL_BIN""/blast_10custom_to_charon_OTU.R \
	$OUTDIR/$PREFIX.OTU.blastOut \
	$OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX.OTU.blast.charon \
	$TAXON_DIR
