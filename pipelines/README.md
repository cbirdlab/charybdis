This directory has example pipeline scripts. These scripts use the programs from GCL_charybdis/bin to assemble pipelines.
Both pipelines assume HPC environment with SLURM. 

**charybdis_generic.sh**
The "default" metabarcoding process. Start here for an overview of the GCL process.

**charybdis_MultiRegion.sh** 
A pipeline with features introduces to handle a specific problem. 
The sequences were a mix of COI, ITS, and 16S. We did not initially know what type each sequence was. 
We needed to split them by type. 
