#!/bin/bash
DB="/work/GenomicSamples/GCL/db/NCBI_BLAST_DBs/coi_its_16s"
GI_list=$1

while read gi; do
        echo -n "$gi,"
        ret=0
    
        STR=$(blastdbcmd -entry $gi -db $DB -outfmt \"%t\")

        if [[ $STR =~ .*COI.* ]] || [[ $STR =~ .*cox1.* ]] || [[ $STR =~ .*COX1.* ]]
        then
                echo "1"
        elif [[ $STR =~ .*ITS.* ]] || [[ $STR =~ .*spacer.* ]] || [[ $STR =~ .*SPACER.* ]]
        then
                echo "2"
        elif [[ $STR =~ .*16S.* ]]
        then
                echo "3"
        else
                echo "0"
        fi

done <$GI_list
