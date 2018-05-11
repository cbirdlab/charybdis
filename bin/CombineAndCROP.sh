#!/bin/bash

	# This script is a work in progress
	# The CROP section is good	
	# Paths need to be updated for HPC


#############
# ARGUMENTS #
#############

INPUT_DATA_FILE=$1 # The FASTA
OUTDIR=$2 # Pipeline output directory
PREFIX=$3 # Project name
CHUNKS=$4 # Number of parallel units
BP=$5 # Base pairs in query region

# Parallel Index Sequence
seq 1 $CHUNKS > $OUTDIR/$PREFIX.loop.dat

#################
# Combine Fasta #
#################

# In order to effectively group OTUs, need to allow the
# clustering program to see all sequences

awk -F',' '{print $2}' $INPUT_DATA_FILE | xargs cat \
    > $OUTDIR/$PREFIX.combo.nonchimeras.clean.fasta

########
# CROP #
########

NUMSEQS=$(grep -c '^>' $OUTDIR/$PREFIX.combo.nonchimeras.clean.fasta)
B_PARAM=$(( $NUMSEQS / 50))
Z_PARAM=$(( 150000 / $BP ))
Z_PARAM=$( echo "$Z_PARAM - .1 * $Z_PARAM" | bc)
Z_PARAM=$(echo "$Z_PARAM" | python -c "print round(float(raw_input()))")
Z_PARAM=${Z_PARAM%.0}
E_PARAM=$((Z_PARAM * 10))

CROPLinux -i $OUTDIR/$PREFIX.combo.nonchimeras.clean.fasta \
	  -o $OUTDIR/$PREFIX.combo.nonchimeras.clean.OTU \
	  -b $B_PARAM \
	  -e $E_PARAM \
	  -z $Z_PARAM \
	  -m $CHUNKS \
	  -r 5 -s 



# Split for parallel 
#gt splitfasta -numfiles $CHUNKS -force \
#	$OUTDIR/$PREFIX.combo.nonchimeras.clean.OTU.cluster.fasta


