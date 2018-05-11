#!/bin/bash

#SBATCH --job-name=AssignToMarkers
#SBATCH --time=24:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mail-user=evan.krell@tamucc.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load bio-misc
module load obitools
module load parallel
module load blast+
module load R

# Arguments
PREFIX=Jones
INDIR=/work/GenomicSamples/GCL/maryJones/in
OUTDIR=/work/GenomicSamples/GCL/maryJones/out
CHUNKS=20

/work/GenomicSamples/GCL/maryJones/scripts/AssignToMarkers.sh \
$PREFIX $INDIR $OUTDIR $CHUNKS


