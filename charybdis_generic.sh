#!/bin/bash

# Pipeline for generic metabarcoding
# For more details, see the GCL Metabarcoding manual
# (should be included wherever this file was found)

# Prerequisites:
# Set up working directory structure

# Define options
PREFIX=""          # -p  Pipeline finds inputs and writes outputs using project name
INDIR=""           # -i  Where input files are stored
OUTDIR=""          # -o  Where to write output files
CHUNKS=1           # -n  Number of instances of parallel execution
BLAST_DB=""        # -b  Path to BLAST database
BLAST_IGNORE=""    # -d  List of NCBI TAXIDs to ignore from BLAST database
VSEARCH_DB=""      # -v  Path to VSEARCH database
SAP_DB=""          # -s  Path to SAP database (SAP-formatted FASTA)
CHIMERA_DB=""      # -c  Path to Chimera database
BASEPAIRS=-1       # -x  Number of basepairs in target genetic region
GCL_BIN=""         # -g  Where the pipeline scripts are stored
TAXON_DIR=""       # -t  Directory with NCBI taxonomic database
ECOTAG_DB=""       # -e  Path to ECOTAG database (not the fasta)
ECOTAG_FASTA_DB="" # -f Path to ECOTAG database (fasta)

# Parse options
while getopts ":p:i:o:n:b:d:v:f:c:s:x:g:t:e:" opt; do
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
		s)
			SAP_DB=$OPTARG
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
		e)
			ECOTAG_DB=$OPTARG
			;;
		f)
			ECOTAG_FASTA_DB=$OPTARG
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


if [ "$BLAST_DB" != "" ]
then
	echo "BLAST path: $BLAST_DB"
	if [ "$BLAST_IGNORE" == "" ]
	then
		echo "Requires '-d': 'File of NCBI TAXIDs to ignore'"
		exit 1
	else
		echo "BLAST_IGNORE: $BLAST_IGNORE"
	fi
fi

if [ "$ECOTAG_DB" != "" ]
then
	echo "ECOTAG database path: $ECOTAG_DB"
	if [ "$ECOTAG_FASTA_DB" == "" ]
	then
		echo "Using '-e' required '-f': 'Path to ECOTAG FASTA'"
	else
		echo "ECOTAG_FASTA_DB: $ECOTAG_FASTA_DB"
	fi
fi

if [ "$VSEARCH_DB" != "" ]
then
	echo "VSEARCH path: $VSEARCH_DB"
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

echo "GCL Charybdis metabarcoding pipeline initiated"
echo "Run with options:"
echo "    $@"


# Make 'samples.txt'
awk '{print $2}' $INDIR""/$PREFIX"".barcodes.txt \
	> $INDIR""/$PREFIX"".samples.txt

# Merge Reads
JOB_ID1=$($GCL_BIN""/sbatch $GCL_BIN""/mergeReads.slurm \
	$PREFIX $INDIR $OUTDIR $CHUNKS $GCL_BIN | grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
echo Submitted job: $JOB_ID1

# Filter Reads
JOB_ID2=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID1 \
	$GCL_BIN""/filterReads.slurm \
	$PREFIX $INDIR $OUTDIR $CHUNKS \
	$BASEPAIRS_LOWER $BASEPAIRS_HIGHER $CHIMERA_DB $GCL_BIN | grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
echo Submitted job: $JOB_ID2

# Cluster into OTUs
JOB_ID3=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID2 \
	$GCL_BIN""/ClusterOTU.slurm \
	$PREFIX $BASEPAIRS $INDIR $OUTDIR $CHUNKS $GCL_BIN | grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
echo Submitted job: $JOB_ID3

# BLAST
if [ "$BLAST_DB" != "" ]
then
	# Taxonomic Assignment with BLAST
	JOB_ID4B=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID3 \
		$GCL_BIN""/BlastMeta.slurm \
 		$PREFIX $INDIR $OUTDIR $CHUNKS 0.7 \
		$BLAST_DB $GCL_BIN $TAXON_DIR $BLAST_IGNORE \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID4B

     	# Create OTUvsTubes
	JOB_ID5B=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID4B \
		$GCL_BIN""/OTUvsTube.slurm \
		$PREFIX $INDIR $OUTDIR $TAXON_DIR $GCL_BIN blast\
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID5B

	# Add descriptive names to OTUvsTubes
	JOB_ID6B=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID5B \
		$GCL_BIN""/OTU_CVT_addSampleDescs.slurm \
		$PREFIX $INDIR $OUTDIR $GCL_BIN blast \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID6B

	# Add BLAST information
	JOB_ID7B=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID6B \
		$GCL_BIN/AddBlast.slurm \
		$PREFIX $INDIR $OUTDIR NCBI $GCL_BIN $TAXON_DIR $CHUNKS \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID7B

fi


# VSEARCH
if [ "$VSEARCH_DB" != "" ]
then

	# Taxonomic Assignment with VSEARCH
	JOB_ID4V=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID3 \
		$GCL_BIN""/Vsearch.slurm \
		$PREFIX $INDIR $OUTDIR 0.7 $VSEARCH_DB $GCL_BIN $TAXON_DIR \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID4V

	# Create OTUvsTubes
	JOB_ID5V=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID4V \
		$GCL_BIN""/OTUvsTube.slurm \
		$PREFIX $INDIR $OUTDIR $TAXON_DIR $GCL_BIN vsearch \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID5V

	# Add descriptive names to OTUvsTubes
	JOB_ID6V=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID5V \
		$GCL_BIN""/OTU_CVT_addSampleDescs.slurm \
		$PREFIX $INDIR $OUTDIR $GCL_BIN vsearch \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID6V

	# Add VSEARCH information
	JOB_ID7V=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID6V \
		$GCL_BIN/AddVsearch.slurm \
        $PREFIX $INDIR $OUTDIR NCBI $GCL_BIN \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID7V
fi

# ECOTAG
if [ "$ECOTAG_DB" != "" ]
then

	# Taxonomic Assignment with ECOTAG
	JOB_ID4E=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID3 \
		$GCL_BIN""/Ecotag.slurm \
		$PREFIX $INDIR $OUTDIR 0.7 $ECOTAG_DB $ECOTAG_FASTA_DB $GCL_BIN $TAXON_DIR \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID4E

	# Create OTUvsTubes
	JOB_ID5E=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID4E \
		$GCL_BIN""/OTUvsTube.slurm \
		$PREFIX $INDIR $OUTDIR $TAXON_DIR $GCL_BIN ecotag \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID5E

	# Add descriptive names to OTUvsTubes
	JOB_ID6E=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID5E \
		$GCL_BIN""/OTU_CVT_addSampleDescs.slurm \
		$PREFIX $INDIR $OUTDIR $GCL_BIN ecotag \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID6E

	# Add ECOTAG information
	JOB_ID7E=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID6E \
		$GCL_BIN/AddEcotag.slurm \
        $PREFIX $INDIR $OUTDIR NCBI $GCL_BIN \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID7E
fi

if [ "$SAP_DB" != "" ]
then
	echo "SAP not yet implemented"

	# Taxonomic assignment with SAP
	JOB_ID4S=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID3 \
		$GCL_BIN""/SAP.slurm \
		$PREFIX $INDIR $OUTDIR 0.95 $SAP_DB \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID4S

	# Create OTUvsTubes
	JOB_ID5S=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID4S \
		$GCL_BIN""/OTUvsTube.slurm \
		$PREFIX $INDIR $OUTDIR $TAXON_DIR $GCL_BIN sap \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID5S

	# Add predator names to OTUvsTubes
	JOB_ID6S=$($GCL_BIN""/sbatch --dependency=afterany:$JOB_ID5S \
		$GCL_BIN""/OTU_CVT_addSampleDescs.slurm \
		$PREFIX $INDIR $OUTDIR sap \
		| grep -oh "[0-9]*" | grep -oh '^[^ ]* ')
	echo Submitted job: $JOB_ID6S

	# Classify SAP results into categories

	# Add SAP information
fi

# Combine method-specific OTUvTube tables into single mega table

# Clean up



