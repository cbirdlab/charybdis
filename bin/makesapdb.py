#!/usr/bin/python3
import csv, argparse, sys, os, json, re
from ete3 import NCBITaxa

program = "makesapdb"
version = "0.1.0"
author = "Evan Krell, Darcy Jones"
date = "23 June 2019"
email = "evan.krell@tamucc.edu"
blurb = (
    "{program}\n"
    "{:<}"
    )
short_blurb = (
    "Convert a NCBI-formatted BLAST database to SAP formatted FASTA database"
    ).format(**locals())
license = (
    '{program}-{version}\n'
    '{short_blurb}\n\n'
    'Copyright (C) {date},  {author}'
    '\n\n'
    'This program is free software: you can redistribute it and/or modify '
    'it under the terms of the GNU General Public License as published by '
    'the Free Software Foundation, either version 3 of the License, or '
    '(at your option) any later version.'
    '\n\n'
    'This program is distributed in the hope that it will be useful, '
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the '
    'GNU General Public License for more details.'
    '\n\n'
    'You should have received a copy of the GNU General Public License '
    'along with this program. If not, see <http://www.gnu.org/licenses/>.'
    ).format(**locals())

def main (infile):

    taxon_ranks = [ "species", "genus", "family", "order", "class", "phylum" ]

    sapTempPath = infile

    ncbi = NCBITaxa()

    with open(sapTempPath, 'r') as fh:
        sapTemp = fh.read().splitlines()

    goodline = True
    for line in sapTemp:
        # If line is a header, manipulate it
        if line[0] == '>':
            goodline = True
            # Extract TAXID
            ID = re.search('\\|[^\\|;]*;', line).group(0)
            TAXID = re.search('[\d]+', ID).group(0)
            try:
                lineage = ncbi.get_lineage(TAXID)
                lineageTaxons = ncbi.get_rank(lineage)
                lineageKeys = list(lineageTaxons.keys())
                lineageNames = ncbi.get_taxid_translator(lineageKeys)
                lineageTaxons = dict((v,k) for k,v in lineageTaxons.items())

                # Generate string containing lineage information
                lineageString = ""
                for tr in taxon_ranks: 
                    try: 
                        tax = lineageTaxons[tr]
                        sciname = lineageNames[tax]
                        lineageString = lineageString + "{}: {}, ".format(tr, sciname)
                    except KeyError:
                        pass 
                lineageString = lineageString[:-2]
                # Update header to be SAP compatable
                line = line.replace("____", lineageString)
            except ValueError:
                 goodline = False

        if goodline:
            print(line)



if __name__ == '__main__':
    argParser = argparse.ArgumentParser()
    argParser.add_argument(
        "-i", "--infile", 
        dest = 'infile', 
        help = ("Path to input file (temp SAP FASTA).")
    )

    args = argParser.parse_args()
    main(**args.__dict__)    


