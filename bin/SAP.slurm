#!/bin/bash

#SBATCH --job-name=SAP
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load sap

PREFIX=$1
INDIR=$2
OUTDIR=$3
THRESHOLD=$4
SAPDB=$5

# SAP
$GCL_BIN""/SAP.sh \
        $PREFIX $INDIR $OUTDIR $THRESHOLD $SAPDB

# Convert to Charon format
Rscript $GCL_BIN""/SAP_to_charon_OTU.R \
	$OUTDIR""/SAP-$PREFIX""/assignments.csv \
	$OUTDIR""/SAP-$PREFIX""/taxon_probabilities.csv \
	$OUTDIR""/$PREFIX"".full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR""/$PREFIX"".OTU.sap.charon
