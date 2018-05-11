#!/bin/bash

#SBATCH --job-name=CROP
#SBATCH --time=24:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mail-user=evan.krell@tamucc.edu
##SBATCH --mail-type=begin
##SBATCH --mail-type=end

module load bio-misc
module load obitools

# Arguments
PREFIX=$1
BP=$2
INDIR=/work/GenomicSamples/GCL/maryJones/in
OUTDIR=/work/GenomicSamples/GCL/maryJones/out
CHUNKS=20

mkdir $OUTDIR/CROP_$PREFIX


#Set up CROP parameters
NUMSEQS=$(grep -c '^>' $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta)
B_PARAM=$(( $NUMSEQS / 50))
Z_PARAM=$(( 150000 / $BP ))
Z_PARAM=$( echo "$Z_PARAM - .1 * $Z_PARAM" | bc)
Z_PARAM=$(echo "$Z_PARAM" | python -c "print round(float(raw_input()))")
Z_PARAM=${Z_PARAM%.0}
E_PARAM=$((Z_PARAM * 10))

#Run CROP
CROPLinux -i $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta \
          -o $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU \
          -b $B_PARAM \
          -e $E_PARAM \
          -z $Z_PARAM \
          -m $CHUNKS \
          -r 5 -s


#Fix SIZE value
cd $OUTDIR/CROP_$PREFIX

/work/GenomicSamples/GCL/maryJones/scripts/CROP_size_fix.sh \
$OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.list $OUTDIR/$PREFIX.full.nonchimeras.clean.fasta \
> $OUTDIR/$PREFIX.full.nonchimeras.clean.OTU.cluster.size_fix








