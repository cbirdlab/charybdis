#!/bin/bash

# Pipeline for generic metabarcoding
# For more details, see the GCL Metabarcoding manual 
# (should be included wherever this file was found)

# Prerequisites:
# Set up working directory structure




# Define options
PREFIX=""       # -p  Pipeline finds inputs and writes outputs using project name
INDIR=""        # -i  Where input files are stored
OUTDIR=""       # -o  Where to write output files
CHUNKS=1        # -n  Number of instances of parallel execution
BLAST_DB=""     # -b  Path to BLAST database
BLAST_IGNORE="" # -d  List of NCBI TAXIDs to ignore from BLAST database
VSEARCH_DB=""   # -v  Path to VSEARCH database
CHIMERA_DB=""   # -c  Path to Chimera database
BASEPAIRS=-1    # -x  Number of basepairs in target genetic region
GCL_BIN=""      # -g  Where the pipeline scripts are stored
TAXON_DIR=""    # -t  Directory with NCBI taxonomic database


# Parse options
while getopts ":p:i:o:n:b:d:v:f:c:s:x:g:t:" opt; do
	case $opt in
		p) 
			PREFIX=$OPTARG
			;;
		i)
			INDIR=$OPTARG
			;;
		o)
			OUTDIR=$OPTARG
			;;
		n)
			CHUNKS=$OPTARG
			;;
		b) 
			BLAST_DB=$OPTARG
			;;
		d)
			BLAST_IGNORE=$OPTARG
			;;	
		v)
			VSEARCH_DB=$OPTARG
			;;
		c)
			CHIMERA_DB=$OPTARG
			;;
		x)
			BASEPAIRS=$OPTARG
			;;
		g)
			GCL_BIN=$OPTARG
			GCL_BIN=${GCL_BIN%/}
			;;
		t)
			TAXON_DIR=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac
done


# Ensure options
if [ "$PREFIX" == "" ]
then
	echo "Requires '-p': 'project prefix'"
	exit 1
else
	echo "Project prefix: $PREFIX"
fi

if [ "$INDIR" == "" ]
then
	echo "Requires '-i': 'input directory'"
	exit 1
else
	echo "Input directory: $INDIR"
fi

if [ "$OUTDIR" == "" ]
then
	echo "Requires '-o': 'output directory'"
	exit 1
else
	echo "Output directory: $OUTDIR"
fi

if [ "$GCL_BIN" == "" ]
then
	echo "Requires '-g': 'directory where GCL metabarcoding scripts stored (bin)'"
	exit 1
else
	echo "GCL bin: $GCL_BIN"
fi

if [ "$BLAST_DB" == "" ] && [ "$VSEARCH_DB" == "" ]
then
	echo "Requires either '-b': 'BLAST path' OR '-v': 'VSEARCH path'"
	exit 1
else
	if [ "$BLAST_DB" != "" ] && [ "$VSEARCH_DB" != "" ]
	then
		echo "Cannot have both '-b' and '-v'"
		exit 1
	elif [ "$BLAST_DB" != "" ]
	then
		echo "BLAST path: $BLAST_DB"
		if [ "$BLAST_IGNORE" == "" ]
		then
			echo "Requires '-d': 'File of NCBI TAXIDs to ignore'"
			exit 1
		else
			echo "BLAST_IGNORE: $BLAST_IGNORE"
		fi
	elif [ "$VSEARCH_DB" != "" ]
	then
		echo "VSEARCH path: $VSEARCH_DB"
	fi
fi

if [ "$CHIMERA_DB" == "" ]
then
	echo "Requires '-c': 'chimera FASTA'"
	exit 1
else
	echo "Chimera FASTA: $CHIMERA_DB"
fi

if [ "$BASEPAIRS" == -1 ]
then
	echo "Requires '-x': 'number of basepairs'"
	exit 1
else
	BASEPAIRS_LOWER=$(($BASEPAIRS - 15))
	BASEPAIRS_HIGHER=$(($BASEPAIRS + 15)) 
	echo "Basepairs: $BASEPAIRS"
fi

if [ "$TAXON_DIR" == "" ]
then
	echo "Requires '-t': 'NCBI taxonomy database directory'"
	exit 1
else
	echo "Taxonomy Database: $TAXON_DIR"
fi


echo "Parallel instances: $CHUNKS"

####
# Make 'samples.txt'
####

awk '{print $2}' $INDIR""/$PREFIX"".barcodes.txt \
	> $INDIR""/$PREFIX"".samples.txt

####
# Merge Reads
####	
JOB_ID1=$(sbatch $GCL_BIN""/mergeReads.slurm \
	$PREFIX $INDIR $OUTDIR $CHUNKS $GCL_BIN | grep -oh "[0-9]*")
echo Submitted job: $JOB_ID1

####
# Filter Reads
####
JOB_ID2=$(sbatch --dependency=afterany:$JOB_ID1 \
	$GCL_BIN""/filterReads.slurm \
	$PREFIX $INDIR $OUTDIR $CHUNKS \
	$BASEPAIRS_LOWER $BASEPAIRS_HIGHER $CHIMERA_DB $GCL_BIN | grep -oh "[0-9]*")
echo Submitted job: $JOB_ID2

####
# Cluster into OTUs
####
JOB_ID3=$(sbatch --dependency=afterany:$JOB_ID2 \
	$GCL_BIN""/ClusterOTU.slurm \
	$PREFIX $BASEPAIRS $INDIR $OUTDIR $CHUNKS $GCL_BIN | grep -oh "[0-9]*")
echo Submitted job: $JOB_ID3

###
# BLAST 
###
if [ "$BLAST_DB" != "" ]
then
	JOB_ID4=$(sbatch --dependency=afterany:$JOB_ID3 \
		$GCL_BIN""/BlastMeta.slurm \
		$PREFIX $INDIR $OUTDIR $CHUNKS 0.7 \
		$BLAST_DB $GCL_BIN $TAXON_DIR $BLAST_IGNORE \
		| grep -oh "[0-9]*")
	echo Submitted job: $JOB_ID4
fi

###
# VSEARCH
###
if [ "$VSEARCH_DB" != "" ]
then
	JOB_ID4=$(sbatch --dependency=afterany:$JOB_ID3 \
		$GCL_BIN""/Vsearch.slurm \
		$PREFIX $INDIR $OUTDIR 0.7 $VSEARCH_DB \
		| grep -oh "[0-9]*")
	echo Submitted job: $JOB_ID4
fi

###
# OTUs & Sample counts
###
JOB_ID5=$(sbatch --dependency=afterany:$JOB_ID4 \
	$GCL_BIN""/OTUvsTube.slurm \
	$PREFIX $INDIR $OUTDIR $TAXON_DIR $GCL_BIN \
	| grep -oh "[0-9]*")
echo Submitted job: $JOB_ID5

###
# BLAST, add stat columns
###
if [ "$BLAST_DB" != "" ]
then
	echo klsdjfhsklafj
fi

###
# VSEARCH, add stat columns
###
if [ "$VSEARCH_DB" != "" ]
then
	JOB_ID6=$(sbatch --dependency=afterany:$JOB_ID4 \
		$GCL_BIN/AddVsearch.slurm \
		$PREFIX $INDIR $OUTDIR \
		| grep -oh "[0-9]*")
	echo Submitted job: $JOB_ID6
fi

