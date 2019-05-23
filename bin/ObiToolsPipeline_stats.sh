### 
### Author: Evan Krell
### Purpose: Collects information related to the various
###     outputs of ObiToolsPipeline.sh

#############
# Arguments #
#############

PREFIX=$1
INDIR=$2
OUTDIR=$3

SAMPLES_FILE=$INDIR/$PREFIX.samples.txt

##########################
# Step 0: Before filters #
##########################

COUNT_FORWARD=$(wc -l $INDIR/$PREFIX.forward.fastq | awk '{print $1}')
COUNT_REVERSE=$(wc -l $INDIR/$PREFIX.reverse.fastq | awk '{print $1}')
COUNT_SAMPLES=$(wc -l $INDIR/$PREFIX.samples.txt | awk '{print $1}')

# Write
echo "[raw data]"
echo "NUM_FORWARD_READS:$COUNT_FORWARD"
echo "NUM_REVERSE_READS:$COUNT_REVERSE"
echo "NUM_SAMPLES:$COUNT_SAMPLES"

#####################################
# Step 1: illuminapairedend/obigrep #
#####################################
# Get number of sequences for total data


COUNT=0
while read i; do
	COUNT=$(( $COUNT + $(grep -c '^+' $OUTDIR/$PREFIX.ali.part-$i.fastq) ))
done < $OUTDIR/$PREFIX.loop.dat

# Write
echo "[illuminapairedend]"
echo "NUM_SEQS:$COUNT"

####################
# Step 2: ngfilter #
####################
# Get number of sequences for total data

COUNT=0
while read i; do
	COUNT=$(( $COUNT + $(grep -c '^+' $OUTDIR/$PREFIX.ngsfilter.part-$i.fastq) ))
done < $OUTDIR/$PREFIX.loop.dat

COUNT_ASSIGNED=$COUNT

COUNT=0
while read i; do
	COUNT=$(( $COUNT + $(grep -c '^+' $OUTDIR/$PREFIX.unidentified.part-$i.fastq) ))
done < $OUTDIR/$PREFIX.loop.dat
COUNT_UNASSIGNED=$COUNT

# Write
echo "[ngsfilter]"
echo "NUM_SEQS_ASSIGNED:$COUNT_ASSIGNED"
echo "NUM_SEQS_UNASSIGNED:$COUNT_UNASSIGNED"

###################
# Step 3: obiuniq #
###################

echo "[obiuniq]"
while read i; do
	COUNT=$(obistat --without-progress-bar $OUTDIR/$PREFIX.ali.assigned.ann.uniq.$i.fasta \
		| tail -n 1 | awk '{print $2}')
	echo "$i:$COUNT"
done < $SAMPLES_FILE

####################
# Step 3: obiclean #
####################

echo "[obiclean]"
while read i; do
	COUNT="$(obigrep $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.$i.fasta \
		| obistat --without-progress-bar -c 'obiclean_count["XXX"]' \
		 | tail -n +1 | awk '{print $1*$2}' | awk '{s+=$1} END {print s}')"
	echo "$i:$COUNT"
done < $SAMPLES_FILE

###########################
# Step 4: obigrep, length #
###########################

echo "[obigrep-len]"
while read i; do
	COUNT="$(obigrep $OUTDIR/$PREFIX.ali.assigned.ann.uniq.clean.len.$i.fasta \
		| obistat --without-progress-bar -c 'obiclean_count["XXX"]' \
		 | tail -n +1 | awk '{print $1*$2}' | awk '{s+=$1} END {print s}')"
	echo "$i:$COUNT"
done < $SAMPLES_FILE


####################
# Step 5: chimeras #
####################

echo "[chimeras]"
while read i; do
	COUNT="$(obigrep $OUTDIR/$PREFIX.$i.chimeras.clean.fasta \
		| obistat --without-progress-bar -c 'obiclean_count["XXX"]' \
		 | tail -n +1 | awk '{print $1*$2}' | awk '{s+=$1} END {print s}')"
	echo "$i:$COUNT"
done < $SAMPLES_FILE

########################
# Step 6: non-chimeras #
########################

echo "[non-chimeras]"
while read i; do
	COUNT="$(obigrep $OUTDIR/$PREFIX.$i.nonchimeras.clean.fasta \
		| obistat --without-progress-bar -c 'obiclean_count["XXX"]' \
		 | tail -n +1 | awk '{print $1*$2}' | awk '{s+=$1} END {print s}')"
	echo "$i:$COUNT"
done < $SAMPLES_FILE




