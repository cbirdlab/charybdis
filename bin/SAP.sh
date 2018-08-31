# Runs SAP using a local SAP-formatted FASTA as database
# For setting up a compatible FASTA, see: https://github.com/cbirdlab/makesapdb
# See the 'caveats' section of that repo's README.md to understand why
# this is not being done in parallel
# Parallel is possible, but is a huge hassle.
### If that caveat is no longer listed, then it should be safe to make parallel. 


PREFIX=$1
INDIR=$2
OUTDIR=$3
THRESHOLD=$5
SAPDB=$6

sap --database $SAPDB --project $OUTDIR""/SAP-$PREFIX"" $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.fasta
