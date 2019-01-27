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


SOURCE_PATH=$(pwd)
cd $GCL_BIN
GCL_BIN_ABS=$(pwd)
cd $SOURCE_PATH
cd $OUTDIR
OUTDIR_ABS=$(pwd)
echo "OUTDIR_ABS: $OUTDIR_ABS"
cd $OUTDIR_ABS/CROP_$PREFIX
$GCL_BIN_ABS""/CROP_size_fix.sh \
	$OUTDIR_ABS/$PREFIX.full.nonchimeras.clean.OTU.cluster.list $OUTDIR_ABS/$PREFIX.full.nonchimeras.clean.fasta \
	> $OUTDIR_ABS/$PREFIX.full.nonchimeras.clean.OTU.cluster.size_fix
cd $SOURCE_PATH

echo "done: Cluster into OTUs"


