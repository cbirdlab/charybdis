#!/bin/bash

#SBATCH --job-name=EcotagAssignment
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

PREFIX=$1
INDIR=$2
OUTDIR=$3
ID=$4
DB_ECO=$5
DB_FASTA=$6
GCL_BIN=$7
TAXDIR=$8
CHUNKS=$9

module load obitools

# ECOTAG
bash $GCL_BIN""/Ecotag.sh \
	$PREFIX $INDIR $OUTDIR $ID $DB_ECO $DB_FASTA $CHUNKS

# Convert assigned FASTA to CSV
obitab $OUTDIR""/$PREFIX"".OTU.ecotag.fasta \
	-a species_list \
	--output-field-separator=, \
	| sed -e '1s/:[^,]*,/,/g' -e 's/ /_/g' \
	> $OUTDIR""/$PREFIX"".OTU.ecotag.csv

# Convert to charon format
Rscript $GCL_BIN""/ecotag_tab_to_charon.R \
    $OUTDIR""/$PREFIX"".OTU.ecotag.csv \
    $OUTDIR""/$PREFIX"".full.nonchimeras.clean.OTU.cluster.size_fix \
    $OUTDIR""/$PREFIX"".OTU.ecotag.charon
