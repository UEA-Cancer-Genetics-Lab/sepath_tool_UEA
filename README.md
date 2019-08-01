# SEPATH - Searching for Pathogens Computational Pipeline

-   [Introduction](#introduction)
    -   [Quick Tutorial](#Quick-Tutorial)
-   [System Requirements](#system-requirements)
    -   [**Dependencies**](#dependencies)
    -   [**Memory Requirements**](#memory-requirements)
    -   [**Disk Space**](#disk-space)
-   [Installation](#installation)
-   [Configuring Cluster Environment Variables & Databases](#configuring-cluster-environment-variables-databases)
    -   [Cluster Environment](#cluster-environment)
    -   [Databases & Adjustable Parameters](#databases-adjustable-parameters)
-   [Tutorials](#tutorials)
    -   [1. Obtaining a List of Raw Input Files](#obtaining-a-list-of-raw-input-files)
    -   [2. Processing BAM Files](#processing-bam-files)
    -   [3. Processing FASTQ Files](#processing-fastq-files)
-   [License](#license)

Authors:

-   Abraham Gihawi - <A.Gihawi@uea.ac.uk>
-   Daniel Brewer - <D.Brewer@uea.ac.uk>

Introduction
------------

SEPATH is a software designed to obtain accurate taxonomic classifications from within host tissue sequences. It is implemented in python 3 and relies on the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system. SEPATH was produced specifically to provide researchers with the ability to conduct ultra high throughput metagenomic studies from raw data to classification and has been benchmarked on simulated human cancer whole genome sequence datasets.

SEPATH provides users with the ability to:

-   Filter host-aligned BAM files to obtain potential non-host sequences
-   Process BAM files or paired end FASTQ files to obtain taxonomic classifications (using Kraken or mOTUs2) with additional outputs that are either assembled contigs or high quality, non-host sequencing reads

The latest updates to SEPATH can be found on [GITHUB](https://github.com/), please feel free to submit feature requests and bug reports. If you have used SEPATH please [cite us](https://google.com).

Gihawi, A. Rallapalli, G. Hurst, R. Cooper, C.S. Leggett, R.M. Brewer, D.S. **SEPATH: Benchmarking the Search for Pathogens in Human Tissue Whole Genome Sequence Data Leads to Template Pipelines**, Manuscript in preparation.

#### Quick Tutorial

A sample pipeline that doesn't make full use of cluster queueing systems has been made available for those that wish to use SEPATH with minimal to no edits required.

In the `tutorials/files/` directory we provide a sample bam files to help get you started.

Ensure you have [Singularity](https://singularity.lbl.gov/install-linux) (in alpha release for Mac at time of writing) and [Git](https://gist.github.com/derhuerst/1b15ff4652a867391f03) installed and configured.

Launch a sample version of SEPATH with:

```
git clone git@github.com:Agihawi/SEPATH.git
singularity pull library://agihawi/default/sepath
singularity shell sepath_latest.sif
./sample_sepath.py
```

This script downloads a minikraken database and standard human reference genomes to the ```./dmp/``` directory. By default, all of the output files will be put into the ```output/``` directory. When the script has stopped running, type ```exit``` to return to your normal shell.

To adapt this for your own analysis, run the pipeline with your own `.bam` files in the files/ directory. all parameters can be set in the first few lines of ```sepath.snake```. You can edit anby databases as required there. Threads are automatically detected, but you can swap this in for any string value. Please refer to our manuscript for recommendations.


System Requirements
-------------------

Due to the nature of SEPATH in providing the glue between tools, there are a few required dependencies that need to be installed and functioning correctly before high throughput analysis can be carried out.

#### **Dependencies**

As a minimum to run SEPATH, we recommend the following:

-   Access to a Linux based high performance computing cluster environment with a job submission scheduler (*i.e* [LSF](https://www.ibm.com/support/knowledgecenter/en/SSETD4/product_welcome_platform_lsf.html) / [SLURM](https://slurm.schedmd.com/))
-   [Python](https://www.python.org/downloads/) v3.5
-   [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) v3.13.3
-   [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.36 and as a result [Java](https://www.oracle.com/technetwork/java/javase/downloads/java-archive-javase8-2177648.html) v1.8.0\_51
-   [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/) part of the bbmap toolset v37.28
-   [mOTUs2](http://motu-tool.org/) v2.0.1
-   [Kraken 1](https://github.com/DerrickWood/kraken) - At the time of writing SEPATH, kraken 2 was in pre-release phase, and so has not been benchmarked. Support for kraken 2 may be added following its stable release.
-   [SPAdes](http://cab.spbu.ru/software/spades/) for metagenomic assembly v3.11.1
-   Optional: [Pysam](https://pypi.org/project/pysam/) v0.15.1 - Only required if you require pre-filtering of BAM files
-   Optional: [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda) - For easier installation of the above dependencies 

#### **Memory Requirements**

There are a few computational steps in SEPATH with high RAM requirements. If using Kraken as part of the pipeline, then 260GB RAM will be required to preload the kraken database into memory. The next most rate limiting step is required by BBDuk. Using our contaminant database with BBDuk, we recommend having access to submission queues with at least 150GB RAM.

#### **Disk Space**

Every effort has been made to minimise disk usage but ultimately the required storage will depend on the size of the input sequencing files and the number of them. From our simulations on cancer tissue, unmapped (non-human) gzip compressed sequencing files are typically &lt;6GB per sample but this is likely to vary for other samples.

Installation
------------

After downloading SEPATH, ensure that the SEPATH directory is in your PATH, you can find this out by entering echo $PATH in the console. Navigate to the directory containing SEPATH, add this to the PATH variable using (this will need to be done each time you start a new session, or alternatively you can add it to your .bashrc file to prevent you having to perform this each time):

	cd SEPATH/
	export PATH=$PATH:$(pwd)

SEPATH Dependencies can be easily installed through conda through an enviromment file by typing the following command in the SEPATH directory:

Ensure conda has channels added 

```

conda config --add channels conda_forge
conda config --add channels bioconda

```

Then create a conda environment for SEPATH

```
conda create --name sepath_conda --file bin/environment.yml
conda activate sepath_conda

```

If installing dependencies manually, please see the website of each dependency listed above for installation instructions. We provide a script *check\_dependencies.py* to check that each dependency is functioning correctly prior to SEPATH use.

This script can be used as follows:

    cd bin/
    chmod +x check_dependencies.py
    ./check_dependencies.py

If every dependency is functioning correctly you should see the message: `All seems to be in working order!` If not, you will be prompted on what dependencies require installation.

Note - Not all modes in SEPATH require all dependencies. They all require python and snakemake.

Bam filtering just requires pysam (python package) to obtain unmapped fastq files. Main bacterial classification (either from BAM input or FASTQ input) requires Trimmomatic, BBDuk, mOTUs2 Assembly/viral classification (either from BAM input or FASTQ input) requires all of the above as well as SPAdes and Kraken (but not mOTUs2)

Configuring Cluster Environment Variables & Databases
-----------------------------------------------------

#### Cluster Environment

Due to the differences in computing clusters between users, some configuration is required to match the users cluster as well as resource availability. We provide sample files for the LSF job scheduling system which SEPATH was developed on and further instructions on how to adapt SEPATH for your own computing cluster are provided here as well as this [tutorial](https://hpc-carpentry.github.io/hpc-python/17-cluster/) and the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration) if further explanation is required.

If you are using the University of East Anglia High Performance Computing Cluster (UEA HPC), we provide a dependencies script which you can source. Navigate to the directory containing SEPATH, then type: `source dependencies_UEA_HPC`

If you are using another cluster, each of the snakemake (ending in .snake) rules in SEPATH\_directory/bin/fastq\_filtering/ will require to suit your cluster. `launchsnake` variables within ouroboros scripts will require editing to suit your submission parameters. For the uea cluster, launchsnake is typically set to the below for our LSF submission cluster. Support for multiple cluster environments may be added in future upon request.

`launchsnake = ('bsub -J fastq_SEPATH -q long-eth -M 10000 -R "rusage[mem=10000]" "snakemake -s fastq_SEPATH.snake --config outdir=%s  --cluster \\"bsub -q {params.queue} {params.cluster}\\" --jobs 20 --latency-wait 1200 --timestamp"' %(outdir))`

`SEPATH.py` will also require setting the `long_job` variable. This is used to submit the ouroborous scripts that will need to last the duration of the analysis. For UEA HPC we use the default:

`long_job = 'bsub -q long-ib -M 10000 -R "rusage[mem=10000]" -J long '`

#### Databases & Adjustable Parameters

Each of the snakemake commands will have adjustable values. As well as the ability to easily edit each of the snakemake rules in the `shell:` section of each rule, Some frequently adapted variables are towards the top of the snakemake files such as minimum\_quality (see below)

-   `contaminant_db` - FASTA format of contaminant sequences, for our analysis we generated a file containing [GRCh38](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz) concatenated with [COSMIC cancer genes](https://cancer.sanger.ac.uk/cosmic/download), but feel free to substitute this with your own fasta for your contaminants (*i.e.* mouse genome).
-   `minimum_quality` - This is the Phred quality score that SEPATH parses to Trimmomatic for quality control. We advise the user adjusts this parameter according to the quality of their data. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a great tool for quick visual analysis of sequencing data quality and can provide you with some insight into how much data you will remove by setting a particular Phred score cut off. We have set the default value low to account for lower quality data.
-   `kraken_db` - The path to the kraken database (only required to set if using Kraken pipeline (assembly)). See Kraken documentation for advice in generating these files. For our analysis we use contig level and above genomes on NCBI RefSeq. Our databases may be made available in future versions of SEPATH.

Tutorials
---------

This section will provide tutorials on how to use some of the main features of SEPATH.

#### 1. Obtaining a List of Raw Input Files

The first step in any multi-step computational analysis is obtaining a curated list of raw input files. While this is simple and routine for those that know how, we provide a script for those that don't (or for the lazier among us): `return_files.py`.

As long as `return_files.py` is in your system PATH (mentioned above), the script can be launched by typing in `return_files.py`

Once this script has launched you will be prompted to enter the directory containing the files that you want to process, enter the absolute path to the files of your interest. Then enter the end of the files that you would like to find within that directory, the files will be saved to input\_files.txt. This script will provide SEPATH with the list of files that it requires (one file per line). Use of this script is demonstrated below.

    #navigate to the tutorial folder from within the SEPATH directory
    $ cd tutorial/
    $ return_files.py
    Enter directory containing files: /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/files/
    Enter extension of files to process (i.e. .fastq.gz): .bam
    4 files found, writing to input_files.txt
    $ cat input_files.txt
    /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/files/001.bam
    /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/files/002.bam
    /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/files/003.bam
    /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/files/004.bam


#### 2. Processing BAM Files

Now we have a list of input BAM files, we can launch SEPATH. We can call the command below to initiate SEPATH which will by default read input files from input\_files.txt

    SEPATH.py --mode bam_bacterial --num_dirs 2 --outdir /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/output/ 
    Running mOTUs2 bacterial classification pipeline
    Output Directory: /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/output/
    Output directory not detected, creating directory: /gpfs/home/user/scratch/Scripts/SEPATH_BETA/tutorial/output/
    mkdir: created directory `/gpfs/home/user/scratch/Scripts/SEPATH_BETA/test_dir/output/'
    Input file list: input_files.txt
    5 files in input file list
    mkdir: created directory `/gpfs/scratch/user/Scripts/SEPATH_BETA/test_dir/workdir1'
    mkdir: created directory `/gpfs/scratch/user/Scripts/SEPATH_BETA/test_dir/workdir2'
    Entering BAM  --> Bacterial classification mode
    Job <260567> is submitted to queue <long-ib>.
    Job <260568> is submitted to queue <long-ib>.

The output files of this script will appear in the output directory specified. This script will output high quality unmapped reads (paired-end reads ending .1.fastq.gz and .2.fastq.gz and single-end reads ending .0.fastq.gz) for each sample and will also output the genus level results for mOTUs2 in motu format and .biom format.

Alternatively, `SEPATH.py --mode bam_viral` can be used, which will perform a Kraken classification for assembled contigs. Kraken results will automatically have a confidence threshold of 0.2 applied and will be produced as output along with the folder generated by metaspades containing the contigs.

#### 3. Processing FASTQ Files

If not processing BAM files directly, SEPATH can be ran on FASTQ files. Launch SEPATH as above, but call either `--mode fastq_bacterial` or `--mode fastq_viral`. the input fastq files are required to be gzip compressed and named in the format: `sample000_R1.fq.gz  sample000_R2.fq.gz`

Alternative naming conventions and varying forms of data compression may be added in future releases of SEPATH. Support for single end reads will not be supported.

License
-------

Copyright (c) Abraham Gihawi, Daniel Brewer, University of East Anglia

SEPATH is a free software that utilises other open source programs. You are permitted to redistribute and/or modify SEPATH under the terms of the GNU General Public License version 3 or later. SEPATH is provided in the hope that it will be useful for research purposes but without any warranty; without even the implied warranty or fitness for a particular purpose. The authors accept no liability for inaccurate or incorrect interpretation of results. Please see the file *LICENSE.txt* for further information.
