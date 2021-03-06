
PREFIX=$1
INDIR=$2
OUTDIR=$3
CHUNKS=$4
THRESHOLD=$5
BLASTDB=$6
GI_IGNORE=$7
BLAST_TASK=$8

seq 1 $CHUNKS > $OUTDIR/$PREFIX.loop.dat
gt splitfasta -force -numfiles $CHUNKS $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta
LIMIT=$(ls $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta.* | wc -l)
seq 1 $LIMIT > $OUTDIR/$PREFIX.loop2.dat


CMD="blastn -db $BLASTDB"
CMD="$CMD -query $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta.{}"
CMD="$CMD -task $BLAST_TASK"
CMD="$CMD -negative_gilist $GI_IGNORE"
CMD="$CMD -outfmt '10 qseqid sseqid sscinames staxids evalue score"
CMD="$CMD bitscore pident qcovs qcovhsp qcovus'"
CMD="$CMD > $OUTDIR/$PREFIX.OTU.$BLAST_TASK.Out.{}"

cat $OUTDIR/$PREFIX.loop2.dat | parallel --no-notice  -j $CHUNKS $CMD
