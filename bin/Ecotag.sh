PREFIX=$1
INDIR=$2
OUTDIR=$3
ID=$4
DB_ECO=$5
DB_FASTA=$6

ecotag -d $DB_ECO \
       -R $DB_FASTA \
       -m $ID \
       -r $OUTDIR""/$PREFIX"".full.nonchimeras.clean.OTU.cluster.fasta \
       > $OUTDIR""/$PREFIX"".OTU.ecotag.fasta
