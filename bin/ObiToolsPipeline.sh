#!/bin/bash

#########      ############      ##########      ##########      ##########
# SPLIT # ---> # PARALLEL # ---> # CONCAT # ---> # DELETE # ---> # SERIAL # ---|
#########      ############      ##########      ##########      ##########    |
   # ^--------------------------------------------------------------------------

#################
# CONSTANT DATA #
#################

CHIMERA_DB=$7
GCL_PATH=$8

#############
# ARGUMENTS #
#############

PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
MARKER_LEN_LOWER=$5
MARKER_LEN_HIGHER=$6

SAMPLES_FILE=$INDIR/$PREFIX.samples.txt

seq 1 $CHUNKS > $OUTDIR/$PREFIX.loop.dat


###
echo "obigrep"
echo "   Split by samples"
###
CMD="obigrep -a sample:^{}$"
CMD="$CMD $OUTDIR/$PREFIX.ann.fasta"
CMD="$CMD --without-progress-bar > $OUTDIR/$PREFIX.ann.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

###
echo "obiannotate"
echo "   Keep only specified attributes"
###
CMD="obiannotate -k sample -k count -k seq_length"
CMD="$CMD -k ali_length -k seq_ab_match"
CMD="$CMD -k seq_a_mismatch -k seq_a_deletion -k seq_a_insertion"
CMD="$CMD -k seq_b_mismatch -k seq_b_deletion -k seq_b_insertion"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/$PREFIX.ann.{.}.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

###
echo "obiuniq"
echo "   Keep only unique sequences"
###
CMD="obiuniq -m sample $OUTDIR/$PREFIX.ali.assigned.ann.{.}.fasta"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.uniq.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

###
echo "obiclean"
echo "   filter PCR errors"
###
CMD="obiclean -r 0.5 -d 1 -H"
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.{.}.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD


echo "obigrep"
echo "   length filter"
###
CMD="obigrep -l $MARKER_LEN_LOWER -L $MARKER_LEN_HIGHER"
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.{.}.fasta"
CMD="$CMD --without-progress-bar"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

#################################################
echo "vsearch - UCHIME (chimera detection )"
#################################################


# The next grouping of lines are to prepare the data for vsearch's chimera detection
# This is because the denovo implementation requires abundance counts using a 'size' attribute.
# This can be pulled directly from the obiclean_count.
# The only issue is that vsearch does not allow for spaces in the sequence ID
# so we remove them. Obitools require spaces, so they will need to be added in before further obi processing.
# It is my hope that someday obitools will not require spaces, since it is already delmited by semicolons.

CMD="obiannotate -S size:'\"%d\" % obiclean_count[\"XXX\"]'"
CMD="$CMD --without-progress-bar"
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.{.}.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

CMD='sed -e "s\ \; \\"'
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.{.}.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace1.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

CMD='sed -e "s\ \\\g"'
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace1.{.}.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

CMD="vsearch --uchime_denovo $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.fasta"
CMD="$CMD --chimeras $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.chimeras.fasta"
CMD="$CMD --nonchimeras $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.nonchimeras.fasta"
CMD="$CMD &> $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.chimeraReport"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

#Put spaces back in for obitools compatability
CMD='sed -e "s\;\; \g" -e "s/;//"'
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.chimeras.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.{.}.chimeras.clean.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD

CMD='sed -e "s\;\; \g" -e "s/;//"'
CMD="$CMD $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.size.nospace.{.}.nonchimeras.fasta"
CMD="$CMD > $OUTDIR/$PREFIX.{.}.nonchimeras.clean.fasta"
cat $SAMPLES_FILE | parallel -j $CHUNKS $CMD


##################################################
echo "Concatenate samples into single fasta"
##################################################
rm $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta
cat $OUTDIR/$PREFIX.*.nonchimeras.clean.fasta > \
	$OUTDIR/$PREFIX.full.nonchimeras.clean.fasta
fasta_formatter -i $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta -t \
	| sed -e 's/[^= ]*=[^;]*;//g' -e 's/>//g' -e 's/[[:blank:]]\+/,/g' \
	> $OUTDIR/$PREFIX.full.nonchimeras.clean.csv

# Get sequence counts for various steps
$GCL_PATH/ObiToolsPipeline_stats.sh $PREFIX $INDIR $OUTDIR > $OUTDIR/$PREFIX.seq_stats.txt

# Get csv of total number of seqs per sample after obigrep-len
awk '/\[obigrep-len\]/{flag=1;next}/\[chimeras\]/{flag=0}flag' $OUTDIR/$PREFIX.seq_stats.txt > \
	$OUTDIR/$PREFIX.total_counts.csv

echo "done: Quality Filter"
