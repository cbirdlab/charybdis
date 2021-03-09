#!/bin/bash

PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
GCL_BIN=$5

SAMPLES_FILE=$INDIR/$PREFIX.samples.txt
seq 1 $CHUNKS > $OUTDIR/$PREFIX.loop.dat

###
echo "fastq-splitter.pl:"
echo "   Splits forward,  reverse reads into chunks for parallel --no-notice "
###
$GCL_BIN/fastq-splitter.pl --n-parts $CHUNKS --measure count --outdir $OUTDIR $INDIR""/$PREFIX""_forward.fastq &
$GCL_BIN/fastq-splitter.pl --n-parts $CHUNKS --measure count --outdir $OUTDIR $INDIR""/$PREFIX""_reverse.fastq &
wait

###
echo "illuminapairedend:"
echo "   Joins forward and reverse reads"
###
CMD="illuminapairedend --score-min=40"
CMD="$CMD $OUTDIR""/$PREFIX""_forward.part-{}.fastq"
CMD="$CMD -r $OUTDIR/$PREFIX""_reverse.part-{}.fastq"
CMD="$CMD > $OUTDIR/$PREFIX"".paired.part-{}.fastq"
cat $OUTDIR/$PREFIX.loop.dat | parallel --no-notice  -j $CHUNKS $CMD

###
echo "obiconvert:"
echo "  convert FASTQ to FASTA..."
###
CMD="obiconvert $OUTDIR""/$PREFIX"".paired.part-{}.fastq"
CMD="$CMD --nuc --fasta-output"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR""/$PREFIX"".illumina.part-{}.fasta"
cat $OUTDIR/$PREFIX.loop.dat | parallel --no-notice  -j $CHUNKS $CMD

###
echo "obigrep"
echo "  remove unaligned sequences..."
###
CMD="obigrep -p 'mode!=\"joined\"' $OUTDIR/$PREFIX.illumina.part-{}.fasta"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/$PREFIX.ali.part-{}.fasta"
cat $OUTDIR/$PREFIX.loop.dat | parallel --no-notice  -j $CHUNKS $CMD

###
echo "obigrep"
echo "   Remove seqs whose alignmentscore is < minimum"
MIN_ALIGN=20
CMD="obigrep -p 'ali_length>=$MIN_ALIGN' $OUTDIR/$PREFIX.ali.part-{}.fasta"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/$PREFIX.ali.2.part-{}.fasta"
cat $OUTDIR/$PREFIX.loop.dat | parallel --no-notice  -j $CHUNKS $CMD

###
echo "ngsfilter"
echo "  Assign sample ID to sequences based on barcode..."
###
CMD="ngsfilter -e 2 -t $INDIR/$PREFIX.barcodes.txt"
CMD="$CMD -u $OUTDIR/$PREFIX.unidentified.part-{}.fasta"
CMD="$CMD $OUTDIR/$PREFIX.ali.2.part-{}.fasta"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/$PREFIX.ann.part-{}.fasta"
cat $OUTDIR/$PREFIX.loop.dat | parallel --no-notice  -j $CHUNKS $CMD

###
echo "cat"
echo "   Concatenate into single FASTA"
###
cat $OUTDIR""/$PREFIX"".ann.part-*.fasta > $OUTDIR""/$PREFIX"".ann.fasta
