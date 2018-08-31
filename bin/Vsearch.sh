PREFIX=$1
INDIR=$2
OUTDIR=$3
ID=$4
DB=$5

vsearch --usearch_global $OUTDIR""/$PREFIX"".full.nonchimeras.clean.OTU.cluster.fasta \
	--db $DB \
	--alnout $OUTDIR""/$PREFIX"".OTU.vsearchAlnout.txt \
	--userfields query+target+id+alnlen+qcov \
	--userout $OUTDIR""/$PREFIX"".OTU.vsearchUserout.txt \
	--blast6out $OUTDIR""/$PREFIX"".OTU.vsearchB6out.txt \
	--id $ID


