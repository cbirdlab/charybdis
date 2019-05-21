# Charybdis
Metabarcoding pipeline and various utility scripts. 

## Overview

Metabarcoding is a method assessing the biodiversity within mixed DNA samples.
Large amounts of environmental DNA are amplified using high-throughput sequencing.
Various filtering and taxonomic assignment techniques are used to identify organisms. 

For more information on metabarcoding:
[Matthieu Leray, Joy Y Yang, Christopher P Meyer, Suzanne C Mills, Natalia Agudelo, Vincent Ranwez, Joel T Boehm and Ryuji J Machida. A new versatile primer set targeting a short fragment of the mitochondrial COI region for metabarcoding metazoan diversity: application for characterizing coral reef fish gut contents. 14 March 2013.](https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34)

For work published using charybdis:
[J. Marcus Drymon Pearce T. Cooper Sean P. Powers Molly M. Miller Sharon Magnuson Evan Krell Chris Bird. Genetic identification of depredator species in commercial and recreational fisheries. 15 April 2019.](https://afspubs.onlinelibrary.wiley.com/doi/10.1002/nafm.10292) 

Charybdis involves numerous scripts in a variety of programming languages. A charybdis project is done by creating a pipeline of these scripts and setting the parameters of each. Because the needs of individual projects differs, it is difficult to create a one-size-fits-all program. Therefore, the use of small, modular programs appears to be the most effective. This manual describes each script in the charybdis pipeline, demonstrates an example pipeline, and offers guidance on setting up databases and running on an HPC.
Multiple assignment methods are supported and my be run simultaneously. For example, a user can supply options to use both BLAST and VSEARCH. At the assignment step, the pipeline splits and both methods are run in parallel. Each method produces a number of resulting files. These files are differentiated by including the assignment method used in the path name. Taxonomic assignment is not required. If no assignment databases are specified, the pipeline can still be used for just the filtering and clustering of reads.

## Manual

A complete manual (PDF) is located at [charybdis_man.pdf](charybdis_man.pdf).

## Quick Start

The following tutorial is for COI barcoding using a local BLAST nucleotide database.
The computing environment may be a standard Linux workstation/laptop or an HPC using SLURM for job scheduling. 

Using Charybdis is associated with a _project_, which has a unique name.
This is required to simplify handling the paths of files created and accessed.
In the following, <\projname\> refers to your unique project name.

When you see paths reference charybdis (i.e.: charybdis/bin),
the absolute path to charybdis is left out since we can't know 
where you downloaded it to. Similarly, \<projname\> as a directory name
could be located anywhere on your system.

<<<<<<< HEAD
# Install dependencies

Dependencies:
- R: 
- BLAST:

        sudo apt-get install r-base sudo apt install ncbi-blast+

=======
>>>>>>> e382172487ddd1222546555b16f387ff1e07326d
### Download, setup BLAST nucleotide database

Note that we download the nucleotide database because we are working with COI sequences.
See the documentation on [BLAST databases](ftp://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) for all options.

Download, decompress nucleotide (nt) database

        cd ..../charybdis/data/blastdb
        wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt* 
<<<<<<< HEAD
        for a in `ls -1 nt*.tar.gz`; do gzip -dc $a | tar xf -; done

[Optional] Update COI filter GI list.
Before running this, you need to see how many GIs exist.
Go to this [link](https://www.ncbi.nlm.nih.gov/nuccore/?term=mitochondria+or+cytochrome+or+coi+or+co1+or+cox1+or+coxi+or+mitochondrial+genome+or+mitochondria+genome) and note the number of database entries. We call that \<NUM\_COI\>.
=======
        tar xvzf nt* 

[Optional] Update COI filter.
Before running this, you need to see how many GIs exist.
Go to this link and note the number of database entries. We call that \<NUM\_COI\>.
>>>>>>> e382172487ddd1222546555b16f387ff1e07326d
Optional only because you could instead use the list provided with Charybdis. 

        cd ..
        bash ../bin/get_mito_coi_gi_list.sh <NUM_COI>
        mv mitochondrial_coi.NCBI_nucl.gi mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi

Filter BLAST database with list of GIs of COI sequences

<<<<<<< HEAD
        blastdb_aliastool -db charybdis/data/blastdb/nt -gilist \
            mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi \
            -dbtype nucl -out blastdb_coi -title "blastdb_coi"

[Optional] Update environmental sequence GI list.
Before running this, you need to see how many GIs exist.
Go to this [link](https://www.ncbi.nlm.nih.gov/nuccore/?term=%22environmental+samples%22%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+sp%5BTitle%5D) and note the number of database entries. We call that \<NUM\_ENV\>.
Optional only because you could instead use the list provided with Charybdis. 

        bash ../bin/get_env_gi_list.sh <NUM_ENV>
        mv env.NCBI_nucl.gi env.NCBI_NT_MAY2019.gi
=======
        blastdb_aliastool -db charybdis/data/blastdb/nt gilist \
            mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi \
            -dbtype nucl -out blastdb_coi -title "blastdb_coi"
        

>>>>>>> e382172487ddd1222546555b16f387ff1e07326d

### Get required data, setup directory

Create directory structure

        mkdir <projname>
        mkdir <projname>/in <projname>/out

### Move your data

You need four pieces of input data: 
- <projname>\_forward.fastq: Forward reads.
    - Example: testrun/TestData_forward.fastq
- <projname>\_reverse.fastq: Reverse reads.
    - Example: testrun/TestData_reverse.fastq
- <projname>.barcodes.txt: Barcodes that label sequence sample.
    - Example: testrun/TestData.barcodes.txt
- <projname>.sampledescs.csv: Arbitrary description of each sample. 
    - Example: testrun/TestData.sampledescs.csv

Your directory should like this:

        <projname>/ 
        ├── out
        ├── in
        │   ├── <projname>.barcodes.txt
        │   ├── <projname>_forward.fasta
        │   ├── <projname>_reverse.fastq
        └── └── <projname>.sampledescs.csv




## Pipelines Directory

This directory has example pipeline scripts which use the programs in charybdis/bin to perform metabarcoding 

  charybdis_generic.sh      The "default" metabarcoding process

  charybdis_MultiRegion.sh  A pipeline with features introduced to handle a specific problem. The sequences were a mix of COI, ITS,   
                            and 16S. We did not initially know what type each sequence was. We needed to split them by type.

## Test Run

See the testrun directory to get started


