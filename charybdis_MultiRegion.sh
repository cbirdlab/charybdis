#/bin/bash

PREFIX=Jones
INDIR=		# !!! Put your input directory
OUTDIR= 	# !!! Put your output directory
CHUNKS=30
BLAST_DB="/work/GenomicSamples/GCL/db/NCBI_BLAST_DBs/coi_its_16s"
NASTY_GIS="/work/GenomicSamples/GCL/db/GI_lists/nasty.NCBI_nucl.gi"


/work/GenomicSamples/ekrell/GCL_charybdis/bin/AssignToMarkers.sh \
	$PREFIX $INDIR $OUTDIR $CHUNKS $BLAST_DB $NASTY_GIS

PREFIX_COI=jones.COI
MARKER_LEN_LOWER_COI=200
MARKER_LEN_HIGHER_COI=222
PREFIX_ITS=jones.ITS
MARKER_LEN_LOWER_ITS=285
MARKER_LEN_HIGHER_ITS=345
PREFIX_16S=jones.16S
MARKER_LEN_LOWER_16S=425
MARKER_LEN_HIGHER_16S=475

CHIMERA_DB="/work/GenomicSamples/GCL/db/BOLD_FULL.fasta"
STAT_PATH="/work/GenomicSamples/ekrell/GCL_charybdis/bin/ObiToolsPipeline_stats.sh"

cat $OUTDIR/*.COI.ann.fasta > $OUTDIR/$PREFIX_COI.ann.fasta

cat $OUTDIR/*.ITS.ann.fasta > $OUTDIR/$PREFIX_ITS.ann.fasta

cat $OUTDIR/*.16S.ann.fasta > $OUTDIR/$PREFIX_16S.ann.fasta

/work/GenomicSamples/ekrell/GCL_charybdis/bin/ObiToolsPipeline.sh \
	$PREFIX_COI $INDIR $OUTDIR $CHUNKS \
	$MARKER_LEN_LOWER_COI $MARKER_LEN_HIGHER_COI \
	$CHIMERA_DB $STAT_PATH 
/work/GenomicSamples/ekrell/GCL_charybdis/bin/ObiToolsPipeline.sh \
	$PREFIX_ITS $INDIR $OUTDIR $CHUNKS \
	$MARKER_LEN_LOWER_ITS $MARKER_LEN_HIGHER_ITS \
	$CHIMERA_DB $STAT_PATH
/work/GenomicSamples/ekrell/GCL_charybdis/bin/ObiToolsPipeline.sh \
	$PREFIX_16S $INDIR $OUTDIR $CHUNKS \
	$MARKER_LEN_LOWER_16S $MARKER_LEN_HIGHER_16S \
	$CHIMERA_DB $STAT_PATH


BASEPAIRS_COI=211
BASEPAIRS_ITS=314
BASEPAIRS_16S=450
/work/GenomicSamples/ekrell/GCL_charybdis/bin/ClusterOTU.sh \
	$PREFIX_COI $BASEPAIRS_COI $INDIR $OUTDIR $CHUNKS
/work/GenomicSamples/ekrell/GCL_charybdis/bin/ClusterOTU.sh \
	$PREFIX_ITS $BASEPAIRS_ITS $INDIR $OUTDIR $CHUNKS
/work/GenomicSamples/ekrell/GCL_charybdis/bin/ClusterOTU.sh \
	$PREFIX_16S $BASEPAIRS_16S $INDIR $OUTDIR $CHUNKS


# BLAST OTUs
/work/GenomicSamples/ekrell/GCL_charybdis/bin/BlastMeta.sh \
	$PREFIX_COI $INDIR $OUTDIR $CHUNKS \
	97 $BLAST_DB $NASTY_GIS
/work/GenomicSamples/ekrell/GCL_charybdis/bin/BlastMeta.sh \
	$PREFIX_ITS $INDIR $OUTDIR $CHUNKS \
	97 $BLAST_DB $NASTY_GIS
/work/GenomicSamples/ekrell/GCL_charybdis/bin/BlastMeta.sh \
	$PREFIX_16S $INDIR $OUTDIR $CHUNKS \
	97 $BLAST_DB $NASTY_GIS

# Combine into single FASTA for each 
cat $OUTDIR/$PREFIX_COI.OTU.blastOut.* \
	> $OUTDIR/$PREFIX_COI.OTU.blastOut
cat $OUTDIR/$PREFIX_ITS.OTU.blastOut.* \
	> $OUTDIR/$PREFIX_ITS.OTU.blastOut
cat $OUTDIR/$PREFIX_16S.OTU.blastOut.* \
	> $OUTDIR/$PREFIX_16S.OTU.blastOut

TAXON_DIR=/work/GenomicSamples/GCL/db/TAXO

# Convert to CHARON
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/blast_10custom_to_charon_OTU.R \
	$OUTDIR/$PREFIX_COI.OTU.blastOut \
	$OUTDIR/$PREFIX_COI.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX_COI.OTU.charon \
	$TAXON_DIR
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/blast_10custom_to_charon_OTU.R \
	$OUTDIR/$PREFIX_ITS.OTU.blastOut \
	$OUTDIR/$PREFIX_ITS.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX_ITS.OTU.charon \
	$TAXON_DIR
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/blast_10custom_to_charon_OTU.R \
	$OUTDIR/$PREFIX_16S.OTU.blastOut \
	$OUTDIR/$PREFIX_16S.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX_16S.OTU.charon \
	$TAXON_DIR

# Create OTUs VS Tubes data sheet
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/crittersVStubes_OTU.R \
	$PREFIX_COI \
	$INDIR/$PREFIX_COI.samples.txt \
	$OUTDIR/$PREFIX_COI.total_counts.csv \
	$OUTDIR/$PREFIX_COI.OTU.charon \
	$OUTDIR/$PREFIX_COI.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX_COI.OTUvsTubes.csv \
	$TAXON_DIR
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/crittersVStubes_OTU.R \
	$PREFIX_ITS \
	$INDIR/$PREFIX_ITS.samples.txt \
	$OUTDIR/$PREFIX_ITS.total_counts.csv \
	$OUTDIR/$PREFIX_ITS.OTU.charon \
	$OUTDIR/$PREFIX_ITS.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX_ITS.OTUvsTubes.csv \
	$TAXON_DIR
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/crittersVStubes_OTU.R \
	$PREFIX_16S \
	$INDIR/$PREFIX_16S.samples.txt \
	$OUTDIR/$PREFIX_16S.total_counts.csv \
	$OUTDIR/$PREFIX_16S.OTU.charon \
	$OUTDIR/$PREFIX_16S.full.nonchimeras.clean.OTU.cluster.size_fix \
	$OUTDIR/$PREFIX_16S.OTUvsTubes.csv \
	$TAXON_DIR



# Add BLAST information and sequence
sed '$!N;s/\n/ /' $OUTDIR/$PREFIX_COI.full.nonchimeras.clean.OTU.cluster.fasta \
	| sed -e 's/>//' -e 's/ /,/' $OUTDIR/$PREFIX_COI.full.nonchimeras.clean.OTU.cluster.csv
sed '$!N;s/\n/ /' $OUTDIR/$PREFIX_ITS.full.nonchimeras.clean.ITS.cluster.fasta \
	| sed -e 's/>//' -e 's/ /,/' $OUTDIR/$PREFIX_ITS.full.nonchimeras.clean.OTU.cluster.csv
sed '$!N;s/\n/ /' $OUTDIR/$PREFIX_16S.full.nonchimeras.clean.16S.cluster.fasta \
	| sed -e 's/>//' -e 's/ /,/' $OUTDIR/$PREFIX_16S.full.nonchimeras.clean.OTU.cluster.csv

Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/OTU_CVT_addBlast.R \
	$OUTDIR/$PREFIX_COI.OTUvsTubes.csv \
	$OUTDIR/$PREFIX_COI.OTU.blastOut \
	NCBI \
	$OUTDIR/$PREFIX_COI.OTUvsTubes.blast.csv \
	$OUTDIR/$PREFIX_COI.full.nonchimeras.clean.OTU.cluster.csv
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/OTU_CVT_addBlast.R \
	$OUTDIR/$PREFIX_ITS.OTUvsTubes.csv \
	$OUTDIR/$PREFIX_ITS.OTU.blastOut \
	NCBI \
	$OUTDIR/$PREFIX_ITS.OTUvsTubes.blast.csv \
	$OUTDIR/$PREFIX_ITS.full.nonchimeras.clean.OTU.cluster.csv
Rscript /work/GenomicSamples/ekrell/GCL_charybdis/bin/OTU_CVT_addBlast.R \
	$OUTDIR/$PREFIX_16S.OTUvsTubes.csv \
	$OUTDIR/$PREFIX_16S.OTU.blastOut \
	NCBI \
	$OUTDIR/$PREFIX_16S.OTUvsTubes.blast.csv \
	$OUTDIR/$PREFIX_16S.full.nonchimeras.clean.OTU.cluster.csv



