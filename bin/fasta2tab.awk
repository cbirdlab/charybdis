#!/usr/bin/awk -f

# Source: terdon
# https://bioinformatics.stackexchange.com/a/2651

{
        if (substr($1,1,1)==">")
    		if (NR>1)
                    	printf "\n%s\t", substr($0,2,length($0)-1)
    		else 
    			printf "%s\t", substr($0,2,length($0)-1)
            else 
                    printf "%s", $0
}END{printf "\n"}
