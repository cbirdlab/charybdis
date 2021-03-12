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

Warning: charybdis is undergoing major renovations. This README is much more up to date than the manual. Once the changes have slowed down, the manual will be updated. 

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

### Install dependencies

Dependencies:
- R: 
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [parallel](https://www.gnu.org/software/parallel/)
- [GenomeTools](https://github.com/genometools/genometools)
- [Obitools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools) 
- [obitools3](https://git.metabarcoding.org/obitools/obitools3/wikis/home)
- [CROP](https://github.com/tingchenlab/CROP): Clustering Sequences for OTU Prediction
- [Vsearch](https://github.com/torognes/vsearch/archive/v2.13.4.tar.gz)

```bash
        sudo apt-get install r-base ncbi-blast+ genometools parallel 
        sudo apt-get install libcurl4-openssl-dev libxml2-dev
        sudo apt-get install python-dev
```

#### Install NCBI BLAST

```bash
	# find newest version here https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 
	# if it is newer than version installed using apt then remove apt installation of blast
	sudo apt --purge remove ncbi-blast+ 
	
	# replace URL, tar file, and directory below with newest version name
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
	tar -xf ncbi-blast-2.11.0+-x64-linux.tar.gz
	sudo cp ncbi-blast-2.11.0+/bin/* /usr/local/bin
	
```

#### Install Obitools

If you've already installed obitools and it's not working, the first step is to remove your old installations so that they don't conflict.

```bash
#this removes obitools installed via apt
sudo apt remove obitools

#if you didn't tell obitools to install somewhere other than where it installs by default, then you should be ok.  Otherwise, you're going to have to track down all the commands that it installs (maybe in /usr/local/bin ?) and `rm` them.
```

Next, make sure you have `python2.7` and `virtualenv`

```bash
sudo apt update
sudo apt upgrade
sudo apt install python2
sudo apt-get install build-essential libssl-dev libffi-dev python-dev
sudo apt install virtualenv
```

We are not enamored with obitools3, yet, so will continue using obitools.  The following steps will help you successfully install them. Steps adapted from Frederic Boyer's [Biostars post](https://www.biostars.org/p/235898/#237117). 

```bash
# create python 2.7 virtual env
cd
mkdir OBI
cd OBI
# if the next line fails, replace path to python below with correct path.  perhaps `/usr/bin/python2`
virtualenv --python=/usr/bin/python2.7 OBI-env

# download, extract source code
wget 'https://git.metabarcoding.org/obitools/obitools/repository/archive.tar.gz?ref=master'
tar -zxvf "archive.tar.gz?ref=master"

# activate the virtual env
source OBI-env/bin/activate

# install a specific (known working) version of sphinx
pip install sphinx==1.4.8

# enter the source
cd obitools-master-*

# build
python setup.py build

# install
python setup.py install

# leave virtualenv
deactivate

# copy binaries to system-wide location
cd ..
cp OBI-env/bin/* /usr/bin/
```

#### Install Obitools3

```bash
    # obitools has been upgraded to work with python3, we are in the process of making sure charybdis is compatible

    # tested this with conda on, installation of conda recommended prior to obitools 3
    sudo apt install python3-venv
    
    cd ~
    mkdir Downloads
    cd Downloads
    git clone https://git.metabarcoding.org/obitools/obitools3.git
    cd obitools3
    python3 -m venv obi3-env
    source obi3-env/bin/activate
    pip3 install cython
    python3 setup.py install
    source obi_completion_script.bash
        
    # test the install
    obi test
    
    # copy to an appropriate directory in the path
    sudo cp -r ../obitools3 /usr/local/bin
    
    # Add the following manually to ~/.bashrc using `nano ~/.bashrc`:
    source /usr/local/bin/obitools3/obi3-env/bin/activate
```

#### Install CROP

```bash
        cd ~/Downloads
        sudo apt-get install libgsl-dev
        git clone https://github.com/tingchenlab/CROP.git
        cd CROP
        make
        sudo cp CROPLinux /usr/local/bin/
```

#### Install Vsearch

```bash
    cd ~/Downloads
	wget https://github.com/torognes/vsearch/archive/v2.15.2.tar.gz
	tar xzf v2.15.2.tar.gz
	cd vsearch-2.15.2
	./autogen.sh
	./configure
	make
	sudo make install  # as root or sudo make install
```

#### Install R

There are many tutorials on installing R, here is [one](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu)


#### Install R packages: `pracma`, `CHNOSZ`, `bold`, `furrr`, `tidyr`, `future`, `taxizedb`

```R
	R
	install.packages(c('pracma', 'CHNOSZ', 'bold', 'furrr', 'tidyr', 'future', 'taxizedb'))
	# ctrl-d to exit R
```

### Install charybdis

```bash
	git clone https://github.com/cbirdlab/charybdis.git
```

### Download, setup BLAST nucleotide database

Note that we download the nucleotide database because we are working with COI sequences.
See the documentation on [BLAST databases](ftp://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html) for all options.

Download, decompress nucleotide (nt) database (Warning: >100 GB!)

```bash
        cd charybdis
        mdkir data
        cd data
        mkdir blastdb
        cd blastdb
        wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt* 
        #for a in `ls -1 nt*.tar.gz`; do gzip -dc $a | tar xf -; done      #serial
	ls -1 nt*.tar.gz | parallel --no-notice "gzip -dc {} | tar xf -"   #parallel
```

Create COI filter GI list.
Before running this, you need to see how many GIs exist.
Go to this [link](https://www.ncbi.nlm.nih.gov/nuccore/?term=mitochondria+or+cytochrome+or+coi+or+co1+or+cox1+or+coxi+or+mitochondrial+genome+or+mitochondria+genome) and note the number of database entries. We call that \<NUM\_COI\>.

        cd ..
        bash ../bin/get_mito_coi_gi_list.sh <NUM_COI>
        mv mitochondrial_coi.NCBI_nucl.gi mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi

Filter BLAST database with list of GIs of COI sequences

        blastdb_aliastool -db charybdis/data/blastdb/nt -gilist \     #in this line, nt is the file prefix, not a dir
            mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi \
            -dbtype nucl -out blastdb_coi -title "blastdb_coi"

Update environmental sequence GI list.
Before running this, you need to see how many GIs exist.
Go to this [link](https://www.ncbi.nlm.nih.gov/nuccore/?term=%22environmental+samples%22%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+sp%5BTitle%5D) and note the number of database entries. We call that \<NUM\_ENV\>.
Optional only because you could instead use the list provided with Charybdis. 

        bash ../bin/get_env_gi_list.sh <NUM_ENV>
        mv env.NCBI_nucl.gi env.NCBI_NT_MAY2019.gi   #change date as necessary

### Try the testrun before running your own data

Goto the testrun directory and run example_charybdis_run.sh
        
	cd testrun
        bash example_charybdis_run.sh

### Get required data, setup directory

Create directory structure (follows same directory structure as testrun)

        mkdir <projname>
        mkdir <projname>/out

### Move your data

You need four pieces of input data: 
- <projname>\_forward.fastq: Forward reads.
    - Example: testrun/TestData_forward.fastq
- <projname>\_reverse.fastq: Reverse reads.
    - Example: testrun/TestData_reverse.fastq
- <projname>.barcodes.txt: Barcodes that label sequence sample.
    - Example: testrun/TestData.barcodes.txt
    - **NOTE: it is critical that only standard IUPAC code be used, use N rather than I for inosines**
    - **Note:** barcodes cannot have 'I' basepair code. I replace with 'N'. 

Your project directory should like this for COI data:

        <projname>/ 
        ├── out
        ├── <projname>.barcodes.txt
        ├── <projname>_forward.fastq
        ├── <projname>_reverse.fastq
        ├── <projname>_charybdis_run.sh
        ├── chimera_db_stub.txt
        ├── env.NCBI_NT_MAY2019.gi  (get from charybdis/data in place of blast_ignore_stub.txt)

Your project directory should like this to search the whole blast database data:

        <projname>/ 
        ├── out
        ├── <projname>.barcodes.txt
        ├── <projname>_forward.fastq
        ├── <projname>_reverse.fastq
        ├── <projname>_charybdis_run.sh
        ├── chimera_db_stub.txt
        ├── blast_ignore_stub.txt

# Run pipeline

### Assignment with BLAST:

The following commands assume that your present working directory (pwd) is the project directory and you're barcoding COI

        bash charybdis_generic.sh \
	    -p <projname> \		#should match project dir name
	    -i . \
	    -o out \
	    -n 2 \    			#n is number of threads
            -x 313 \
	    -g ../bin \
            -t ../data/taxo \
            -b ../data/blastdb_coi \ 	#in this line, the blastdb_coi is a file prefix, not a dir
            -d env.NCBI_NT_MAY2019.gi \
            -c chimera_db_stub.txt  	#!!! NEED TO REMOVE THIS CHIMERA DEPENDENCY !!!
	    
The following commands assume that your present working directory (pwd) is the project directory and you're barcoding from the whole blast database

        bash charybdis_generic.sh \
	    -p <projname> \		#should match project dir name
	    -i . \
	    -o out \
	    -n 2 \    			#n is number of threads
            -x <length of target amplicon> \
	    -g ../bin \
            -t ../data/taxo \
            -b ../data/blastdb/nt \ 	#in this line, the nt is a file prefix, not a dir
            -d blast_ignore_stub.txt \
            -c chimera_db_stub.txt   	#!!! NEED TO REMOVE THIS CHIMERA DEPENDENCY !!!


### Assignment with VSEARCH:

See further details in this README for creating VSEARCH database

        bash charybdis_generic.sh \
            -p <projname> -i in -o out -n 20 \
            -x 313 -g charybdis/bin \
            -t charybdis/data/taxo \
            -v charybdis/data/vsearchdb_coi_clean.fasta \
            -c !!! NEED TO REMOVE THIS CHIMERA DEPENDENCY !!!

### Assignment with ECOTAG:

See further details in this README for creating ECOTAG database

        bash charybdis_generic.sh \
            -p <projname> -i in -o out -n 20 \
            -x 313 -g charybdis/bin \
            -t charybdis/data/taxo \
            -e charybdis/data/ecotagdb_coi_clean \
            -f charybdis/data/ecotagdb_coi_clean.fasta \
            -c !!! NEED TO REMOVE THIS CHIMERA DEPENDENCY !!!


## Pipelines 

Example pipeline scripts which use the programs in charybdis/bin to perform metabarcoding are: 

  charybdis_generic.sh      The "default" metabarcoding process

  charybdis_MultiRegion.sh  A pipeline with features introduced to handle a specific problem. The sequences were a mix of COI, ITS,   
                            and 16S. We did not initially know what type each sequence was. We needed to split them by type.

## Test Run

See the testrun directory to get started

## Bin Directory

The following are programs that are part of Charybdis.
Pipelines are assembled from these programs.

### Metabarcoding 
Those that deal directly with the metabarcoding process.

- **mergeReads.slurm:** SLURM script to merge forward, reverse reads.
- - **mergeReads.sh:** Merge forward, reverse reads.
- - **fastq-splitter.pl:** Split fastq files for parallel.
- **filterReads.slurm:** SLURM script to filter reads by quality.
- - **ObiToolsPipeline.sh:** Use obitools to filter reads.
- - **ObiToolsPipeline_stats.sh:** Get stats on filtering results.
- **ClusterOTU.slurm:** SLURM script to cluster reads into OTUs.
- - **ClusterOTU.sh:** Cluster reads into OTUs.
- - - **CROP_size_fix.sh:** OTUs have a count of reads in OTU. But needs to include counts done for noise-based cluster in filterReads.slurm.
- **BlastMeta.slurm:** SLURM script to run BLAST.
- - **BlastMeta.sh:** Run BLAST.
- - **blast_10custom_to_charon_OTU.R:** Convert BLAST results to our charon CSV format.
- **OTUvsTube.slurm:** SLURM script to convert charon CSV to a OTUvsTubes CSV format.
- - **crittersVStubes_OTU.R:** Convert charon CSV to a OTUvsTubes CSV format.
- **AddBlast.slurm:** SLURM script to add BLAST scores to OTUvsTubes.
- - **OTU_CVT_addBlast.R:** Add BLAST scores to OTUvsTubes.

### Database utilities

- **get_mito_coi_gi_list.sh:** Get list of GIs containing COI sequences. Used to generate COI-only database.
- **get_env_gi_list.sh:** Get list of GIs containing environmental sequences. used to filter env from database.

### Misc utilities

- **sbatch:** alias for SLURM's sbatch command. If real sbatch exits on system, use it. Otherwise, default to bash.

### Archived (need to be updated)

- **Vsearch.slurm**
- - **Vsearch.sh**
- - **vsearch_getTAXIDfromBOLDseqid.sh**
- - **vsearch_blast6custom_to_charon_OTU.R**
- **AddVsearch.slurm:** Vsearch output stats to OTUvTubes file.
- **vsearch_blast6custom_to_charon-BLASTDB_OTU.R:** Conversion from Vsearch's 'blast6' format to charon, when using BLAST db converted to Vsearch db.
- **SAP.slurm**
- - **SAP.sh**
- **SAP_to_charon_OTU.R**
- **AssignToMarkers.slurm**
- - **AssignToMarkers.sh**
- - - **determine_marker.R**
### Archived (replaced by other scripts)

- **CROP.slurm:** Replaced by ClusterOTU.slurm
- **CombineAndCROP.sh:** Replaced by ClusterOTU.slurm
- **determine_marker.sh:** Replaced by determine_marker.R

## Convert BLAST database to VSEARCH database

Create a VSEARCH-compatable FASTA from a BLAST database.

Format: 

	>ACCESSION|GI|TAXID
	CTTATACTTTCTAGTTGGAATCTGAACAGGACT

Assumes that BLAST database was made from NCBI's FTP-provided version,
as described in the Quick Start above.

Unlike BLAST, we cannot use an "ignore" file.
So we have to make the database with both the coi and env filtering.
This tutorial assumes the same coi, env-filtering scenario as the quick start.
Hopefully it is clear what is being done so you can adapt for your needs. 
If you are not using an ignore list, skip to the 'blastdbcmd' step. 


	cd charybdis/data

	# Get list of desired GIs
	# Here we get only coi that are non-environmental
	grep -F -v -f env.NCBI_NT_<MONTHYEARETC>.gi \
		mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi
		> coi_clean_NT_<MONTHYEARETC>.gi 

	# Create coi, env-free BLAST database
	blastdb_aliastool -db blastdb_coi \
		-gilist coi_clean_NT_<MONTHYEARETC>.gi \
		-dbtype nucl -out blastdb_coi_clean -title "blastdb_coi_clean"

	# Format contents as FASTA
	blastdbcmd -entry all -db blastdb_coi_clean -outfmt ">%a|%g|%T,%s" \
		| tr , '\n' > vsearchdb_coi_clean.fasta


## Convert BLAST database to ECOTAG database. 

Create a ECOTAG-compatable database from a BLAST database.
	
Format:

	>ACCESSION|GI|TAXID taxid=TAXID;
	CTTATACTTTCTAGTTGGAATCTGAACAGGACT
	
Unlike BLAST, we cannot use an "ignore" file.
So we have to make the database with both the coi and env filtering.
This tutorial assumes the same coi, env-filtering scenario as the quick start.
Hopefully it is clear what is being done so you can adapt for your needs. 
If you are not using an ignore list, skip to the 'blastdbcmd' step. 

	cd charybdis/data

	# Get list of desired GIs
	# Here we get only coi that are non-environmental
	grep -F -v -f env.NCBI_NT_<MONTHYEARETC>.gi \
		mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi
		> coi_clean_NT_<MONTHYEARETC>.gi 

	# Create coi, env-free BLAST database
	blastdb_aliastool -db blastdb_coi \
		-gilist coi_clean_NT_<MONTHYEARETC>.gi \
		-dbtype nucl -out blastdb_coi_clean -title "blastdb_coi_clean"

	# Format contents as FASTA
	blastdbcmd -entry all -db blastdb_coi_clean \
		 -outfmt ">%a|%g|%T taxid=%T;,%s" \
		| tr , '\n' > ecotagdb_temp_coi_clean.fasta

	# Format taxonomy database for OBITools
	obitaxonomy -t taxo/ -d taxo/

	# Add taxonomy to ecotagdb fasta
	obiaddtaxids -d taxo/ ecotagdb_temp_coi_clean.fasta \
		> ecotagdb_coi_clean.fasta

	# Convert ecotagdb fasta to native ecotag database format
	obiconvert -d taxo/ --fasta \
		--ecopcrdb-output=ecotagdb_coi_clean ecotagdb_coi_clean.fasta

And now, how to use ecotag for assignment:

	ecotag -d ecotagdb_coi_clean -R ecotagdb_coi_clean.fasta \
		-m 0.95 -r <query fasta> 

## Convert BLAST database to SAP database. 

Create a SAP-compatable database from a BLAST database.
	
Format:

	>ACCESSION|GI|TAXID ; rank1: name1, rank2: name2, ...  ; Sequence title 
	CTTATACTTTCTAGTTGGAATCTGAACAGGACT

Example: 

	>NM_176613.2|31341268|9913 ; species: Bos taurus, genus: Bos, family: Bovidae, class: Mammalia, phylum: Chordata ; Bos taurus ATP synthase membrane subunit c locus 2 (ATP5MC2), mRNA
	TGCCGTAGCCCTTCATCCCCTGAAAATGTACACTTGCGCCAAGTTCGTCTCCACCCCCTCCTTGATCAGGAGAACCTCTACAGTATT


Unlike BLAST, we cannot use an "ignore" file.
So we have to make the database with both the coi and env filtering.
This tutorial assumes the same coi, env-filtering scenario as the quick start.
Hopefully it is clear what is being done so you can adapt for your needs. 
If you are not using an ignore list, skip to the 'blastdbcmd' step. 

	cd charybdis/data

	# Get list of desired GIs
	# Here we get only coi that are non-environmental
	grep -F -v -f env.NCBI_NT_<MONTHYEARETC>.gi \
		mitochondrial_coi.NCBI_NT_<MONTHYEARETC>.gi
		> coi_clean_NT_<MONTHYEARETC>.gi 

	# Create coi, env-free BLAST database
	blastdb_aliastool -db blastdb_coi \
		-gilist coi_clean_NT_<MONTHYEARETC>.gi \
		-dbtype nucl -out blastdb_coi_clean -title "blastdb_coi_clean"

	# Update taxonomy database
	# Sadly, redundant with setting up taxonomy database previously,
	# but now python library ete3 has very specific setup it looks for. 
        ../bin/ete3_ncbi_update_taxonomy.py

	# Convert coi, env-free BLAST database to (mostly) suitable FASTA: temp SAP db
 	# Another step will replace the ____ placeholder with taxonomic lineage
	blastdbcmd -entry all -db blastdb_coi_clean \ 
	    -outfmt ">%a|%g|%T ; ____ ; %a|%g|%T####%s" | sed 's/####/\n/' > sapdb_coi_clean_temp.fasta

	# Insert taxonomic lineage into temp SAP db.
	# Note that some NCBI TAXIDs will not be found, despite using NCBI taxonomy browser..
	# In my experience, these are not useful sequences so we get an additional env/junk filter...
	../bin/makesapdb.py -i sapdb_coi_clean_temp.fasta > sapdb_coi_clean.fasta

And now, how to use SAP for assignment.

	# Run SAP
	sap --database sapdb_coi_clean.fasta --project <project name> <query_fasta>

Things to note about SAP database creation.

Occassional FASTA entries will fail during SAP database creation, but these will be skipped.
This is because of special characters present in the scientific names in the taxonomy database.
Typically useless sequences, but this could be fixed by removing special characters in makesapdb.py.
	
	
	


