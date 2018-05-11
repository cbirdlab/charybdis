#!/bin/bash

VSEARCH_B6=$1
TAXDIR=$2

awk '{print $2}' $VSEARCH_B6 | awk -F'|' '{print $2}' | sed -e 's/_/ /g' | \
	while read line; do \
	echo -n $line"",; \
	grep -m 1 "$line" $TAXDIR""/names.dmp | awk '{print $1}'
	done


