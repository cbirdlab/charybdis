PREFIX=$1
BP=$2
INDIR=$3
OUTDIR=$4
CHUNKS=$5
GCL_BIN=$6

mkdir $OUTDIR/CROP_$PREFIX

NUMSEQS=$(grep -c '^>' $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta)
B_PARAM=$(( $NUMSEQS / 50))
Z_PARAM=$(( 150000 / $BP ))
Z_PARAM=$( echo "$Z_PARAM - .1 * $Z_PARAM" | bc)
Z_PARAM=$(echo "$Z_PARAM" | python -c "print round(float(raw_input()))")
Z_PARAM=${Z_PARAM%.0}
E_PARAM=$((Z_PARAM * 10))

####
echo "CROP (Clustering reads into OTUs)"
####

CROPLinux -i $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta \
	-o $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU \
	-b $B_PARAM \
	-e $E_PARAM \
	-z $Z_PARAM \
	-m $CHUNKS \
	-r 5 -s

####
echo "OTU size correction"
####

cd $OUTDIR/CROP_$PREFIX
$GCL_BIN""/CROP_size_fix.sh \
	$OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.list $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta \
	> $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.size_fix

echo "done: Cluster into OTUs"


