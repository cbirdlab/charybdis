#!/bin/bash

#SBATCH --job-name=ObiToolsPipeline
#SBATCH --time=24:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mail-user=evan.krell@tamucc.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load bio-misc
module load obitools
module load parallel
module load vsearch

# Add your data files here
PREFIX=$1
MARKER_LEN_LOWER=$2
MARKER_LEN_HIGHER=$3
INDIR=/work/GenomicSamples/GCL/maryJones/in
OUTDIR=/work/GenomicSamples/GCL/maryJones/out
CHUNKS=20

#/work/GenomicSamples/GCL/maryJones/scripts/ObiToolsPipeline.sh \
#$PREFIX $INDIR $OUTDIR $CHUNKS $MARKER_LEN_LOWER $MARKER_LEN_HIGHER






