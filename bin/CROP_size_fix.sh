## This script matches up OTU cluster size with the actual size
## This is an issue because the obitools-produced fasta file already may include merged sequences

## EX: If your fasta has the following sequences with merged size attribute:
##     SEQ1:10  SEQ2:7  SEQ3:2
## And if the are clustered as a single OTU such that:
##     SEQ1 = {SEQ1, SEQ2, SEQ3}
## Then Crop will consider the size of SEQ1 as:
##     | SEQ1 | = 3
## When the actual size should include the pre-existing merged size attributes:
##    | SEQ1 | = 10 + 7 + 2 = 19

## Resulting csv should be of format:
##    OTU1,size1
##    OTU2,size2
##    OTU3,size3


CROP_LIST=$1
FASTA=$2

# Make files that contain all sequence IDs in each OTU
while read p; do
    FHandle=`echo $p | awk '{print $1}'`
    echo $FHandle > $FHandle.cluster.list
    echo $p | awk '{print $2}' | sed -e 's/,/\n/g' >> $FHandle.cluster.list
done <$CROP_LIST

# Separate into uniq fastas per csv
#while read p; do
#    FHandle=`echo $p | awk '{print $1}'`
#    SIZE=`obigrep --id-list $FHandle.cluster.list $FASTA --without-progress-bar | grep -e 'size=[[0-9]*' -o | sed -e 's/size=//' | awk '{ SUM += $1} END {print SUM}'`
#    echo "$FHandle,$SIZE"
#done <$CROP_LIST

# Build files where each OTU has a file containing samples and counts
while read p; do
    FHandle=`echo $p | awk '{print $1}'`
    obigrep --id-list $FHandle.cluster.list $FASTA --without-progress-bar | \
        # Grep for the Sample ID and the size 
        grep -e 'merged_sample={[^:]*:' -e 'size=[0-9]*' -o | awk '!(NR%2){print$0p}{p=$0}' | \
	# Sed to clean up extra characters
	sed -e "s/merged_sample={'//" -e "s/':size=/,/" > $FHandle.cluster.samples.txt
done <$CROP_LIST 

# This builds the TABLE
#while read p; do
#    FHandle=`echo $p | awk '{print $1}'`
#    obigrep --id-list $FHandle.cluster.list $FASTA --without-progress-bar | \
#        # Grep for the Sample ID and the size 
#        grep -e '^>[^ ]* ' -e 'merged_sample={[^:]*:' -e 'size=[0-9]*' -o  | awk '{ORS = NR%3 ? "," : "\n"}1' | \
#	# Sed to clean up extra characters
#	sed -e "s/merged_sample={'//" -e "s/size=//" -e "s/'://" -e "s/ //g" -e 's/>//' | sed -e "s/^/$FHandle,/"
#done <$CROP_LIST 

>SEQID_SAMPLEID.temp
>SIZE.temp
# This builds the TABLE
while read p; do
    FHandle=`echo $p | awk '{print $1}'`
    obigrep --id-list $FHandle.cluster.list $FASTA --without-progress-bar | \
	# Grep for OTU Sequence ID and the Sample ID
    	grep -e '>[^ ]* ' -e 'merged_sample={[^:]*:' -o | awk '{ORS = NR%2 ? "," : "\n"}1' | \
		sed -e "s/merged_sample={'//" -e "s/'://" -e 's/merged_sample={//' -e 's/:$//' -e "s/^>/$FHandle,/g" -e 's/ //' >> SEQID_SAMPLEID.temp
    obigrep --id-list $FHandle.cluster.list $FASTA --without-progress-bar | \
        	grep -e 'size=[0-9]*' -o | sed -e 's/size=//'  >> SIZE.temp
done <$CROP_LIST 

paste SEQID_SAMPLEID.temp SIZE.temp -d ','


