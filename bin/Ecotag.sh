PREFIX=$1
INDIR=$2
OUTDIR=$3
ID=$4
DB_ECO=$5
DB_FASTA=$6
CHUNKS=$7

# Hack: copy over the fasta.. in case this is being done in parallel with Blast which also
#       splits files. This way, splitfile output will have different name -> no overwriting
seq 1 $CHUNKS > $OUTDIR/$PREFIX.loop.eco.dat
cp $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta \
   $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta.eco
gt splitfasta -force -numfiles $CHUNKS $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta.eco
LIMIT=$(ls $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta.eco.* | wc -l)
seq 1 $LIMIT > $OUTDIR/$PREFIX.loop2.eco.dat 

CMD="ecotag -d $DB_ECO"
CMD="$CMD -R $DB_FASTA"
CMD="$CMD -m $ID"
CMD="$CMD -r $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta.eco.{}"
CMD="$CMD > $OUTDIR/$PREFIX.OTU.ecotag.fasta.{}"

cat $OUTDIR/$PREFIX.loop2.eco.dat | parallel --no-notice -j $CHUNKS $CMD

cat $OUTDIR/$PREFIX.OTU.ecotag.fasta.{} > $OUTDIR/$PREFIX.OTU.ecotag.fasta
