#!/bin/bash

PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
BLAST_DB=$5
NASTY_GIS=$6

SAMPLES_FILE=$INDIR/$PREFIX.samples.txt

seq 1 $CHUNKS > $OUTDIR/$PREFIX.loop.dat

###
echo "illuminapairedend"
### merge forward and reverse reads
###

CMD="illuminapairedend --score-min=40"
CMD="$CMD $INDIR/{}_forward.fastq"
CMD="$CMD -r $INDIR/{}_reverse.fastq"
CMD="$CMD > $OUTDIR/{}.paired.fastq"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

###
echo "obiconvert"
### convert FASTQ to FASTA
###

CMD="obiconvert $OUTDIR/{}.paired.fastq"
CMD="$CMD --nuc --fasta-output"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/{}.ali.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

###
echo "obigrep"
###  remove unaligned sequences
###

CMD="obigrep -p 'mode!=\"joined\"' $OUTDIR/{}.ali.fasta"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/{}.ali.len.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

###
echo "BLAST"
### Assign the reads to determine marker
###

CMD="blastn -db $BLAST_DB"
CMD="$CMD -query $OUTDIR/{}.ali.len.fasta"
CMD="$CMD -outfmt '10 qseqid sseqid sscinames staxids stitle'"
CMD="$CMD -negative_gilist $NASTY_GIS"
CMD="$CMD > $OUTDIR/{}.blastOut"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

CMD="awk -F',' '{print \$1\",\"\$2\",\"\$3\",\"\$4\",\"\$5}' $OUTDIR/{}.blastOut"
CMD="$CMD > $OUTDIR/{}.blastOut.fix"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

###
echo "determine_marker.R"
### Get GIs from BLAST results
###

CMD="Rscript /media/Keaolani/charybdis/pipeline_hpc/determine_marker.R"
CMD="$CMD $OUTDIR/{}.blastOut.fix"
CMD="$CMD $OUTDIR/{}.seqids_COI.txt $OUTDIR/{}.seqids_ITS.txt $OUTDIR/{}.seqids_16S.txt"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

###
echo "obigrep"
### split FASTA by sequence IDs
###

CMD="obigrep --id-list $OUTDIR/{}.seqids_COI.txt"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/{}.ali.len.fasta"
CMD="$CMD > $OUTDIR/{}.COI.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

CMD="obigrep --id-list $OUTDIR/{}.seqids_ITS.txt"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/{}.ali.len.fasta"
CMD="$CMD > $OUTDIR/{}.ITS.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

CMD="obigrep --id-list $OUTDIR/{}.seqids_16S.txt"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/{}.ali.len.fasta"
CMD="$CMD > $OUTDIR/{}.16S.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

##
echo "obiannotate"
## attach sample and marker information to each sequence
##

CMD="obiannotate -S sample:'{}' -S marker:'COI'"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/{}.COI.fasta"
CMD="$CMD > $OUTDIR/{}.COI.ann.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

CMD="obiannotate -S sample:'{}' -S marker:'ITS'"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/{}.ITS.fasta"
CMD="$CMD > $OUTDIR/{}.ITS.ann.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

CMD="obiannotate -S sample:'{}' -S marker:'16S'"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/{}.16S.fasta"
CMD="$CMD > $OUTDIR/{}.16S.ann.fasta"

cat $SAMPLES_FILE | parallel --no-notice  -j $CHUNKS $CMD

echo "done: AssignByMarkers"
