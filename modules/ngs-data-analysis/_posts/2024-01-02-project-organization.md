---
title: Prerequisite
---

## Prerequisite
This module will come after the introductory Linux module and therefore assumes
familiarity with basic Linux command-line use. It also assumes you have an
account and are operating in the ILRI computing cluster from a local Linux
environment.

>**Note**
>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the
>hpc username that you were assigned, for example `Bio4Info$$`. Your username,
>by default, is stored in a variable called `USER`. By using it, you will not
>have to type-in your username, rather, your shell will automatically pick your
>username which is the value stored in the `USER` variable. The `$` (dollar)
>character-prefix to a variable name is used to call the value of that variable.

## Set-Up
We will use the computer lab at ILRI, which is already equipped with
Linux-operating desktop computers. Since we will be working from the remote
servers, we will not need special setup for personal laptops. However, toward
the end of the program, we can look into access to a Linux server from a Windows
PC; or how to install a Linux (sub)system for any interested persons.


## Log into the HPC
From the terminal (or equvalent tool) of your local computer, you can log into
the HPC using the folowing command line, followed by pressing <ENTER>. You will
be promted to type-in your password (the password will not be visible as you
type it; just have faith). On a Linux system, you can use the `Ctrl-Alt-T`
keyboard shortcut to open a terminal.
`ssh <user_name>@hpc.ilri.cgiar.org`

The HPC head node has 4 CPUs and we need to access more CPUs/resources in other
compute nodes.
You will have to move from the cluster's head node into the node where we will
be working from (it is called `compute05`). Use the following command; `-w`
requests (a) specific list of host(s).
```
interactive -w compute05
```

`ssh` allows you to securely connect to the remote computer over internet, while
`interactive` allows you to reserve resources to work interactively in a
specified node within the computing cluster using the `-w` flag.

>**Note**
>When running a job interactively, the time limit is 8 hours and Default number of CPU is 1.

## Project organisation

1. We then change into the `compute05` `scratch` directory to create our project
   directory. Using the option`-p` (parent) `mkdir` will create any missing
   intermediate directories.
    ```
    cd /var/scratch/
    mkdir -p $USER/AfricaCDC_training
    cd $USER/AfricaCDC_training
    ```
2. The `assets, databases, primer-schemes, scripts` directories will be linked
   to the project directory, to limit redundancy. `-s` (soft) means that we are
   creating a soft link.
    ```
    ln -s /var/scratch/global/AfricaCDC_training/[adps]* .
    ```
3. We will create the directories `data`, `results` and `genome` to store raw
   data in ```fastq``` format, output and reference genomes respectively.
   Intermediate output files per `tool/software` will be created within the
   `results` directory. We will exploit the bash array data structure to create
   all the directories at once.
    ```
    mkdir data genome results
    mkdir -p results/{fastqc,fastp,kraken,samtools,ivar,snpeff,pangolin,nextclade,multiqc,bowtie2,bedtools}
    ```
4. Change into the `data` directory, from where we will retrieve our ```fastq``` files.
    ```
    cd data
    ls
    ```
## Data retrieval and integrity checks
1. While there are specialised tools for data retieval from nucleotide sequence
   databases, universal `Unix` command (`wget`) can be used to download data
   over internet.
    ```
    wget --no-check-certificate https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars-fastqs.tar.gz
    ```
2. After downloading your data, say from a sequencing facility site, it is often
   good practice to verify that your data was not intentionally/accidentally
   tampered with. To do this, your data service provider will likely accompany
   your data with a file containing a verification code: `checksum_file`
   (***will be provided***). The `md5sum` command, using the `-c` (check) tag,
   allows for checking the integrity of a file downloaded or acquired from a
   different source.
    ```
    wget --no-check-certificate https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars-fastqs.md5
    ls
    md5sum -c sars-fastqs.md5
    ```
3. Next, we will unzip the file using `tar` with the `-xf` (extract, file;
   respectively) tags, which tells `tar` extract the given file.
    ```
    tar -xf sars-fastqs.tar.gz
    ls
    ```
4.  Download SARS-CoV-2 reference genome and the genome annotation file.

    We will retrieve SARS-CoV-2 reference genome and the annotation from [NCBI](https://www.ncbi.nlm.nih.gov/).
    1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
    2. Type 'SARS-CoV-2' on the search box and select 'Genome' database.
    3. Select the [Genbank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/) hyperlink.
    4. Select the genome version [GCA_009858895.3_ASM985889v3](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/).
    5. Right click on the genome [FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz) and select 'copy link'.
    6. Change into the ```genome``` directory using the command
    ```cd ../genome```.
    7. Use ```wget``` to fetch the file.
    8. Retrieve the feature annotation file [GFF](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz) using ```wget``` command.
    9. Dowload the [md5checksum](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/md5checksums.txt)  using `wget` command and check for integrity of your reference genome (FASTA) and annotation (GFF) files.

        ```
        echo "$(grep *GCA_009858895.3_ASM985889v3_genomic.fna.gz* md5checksums.txt | cut -f1 -d' ') GCA_009858895.3_ASM985889v3_genomic.fna.gz" | md5sum -c -
        echo "$(grep *GCA_009858895.3_ASM985889v3_genomic.gff.gz* md5checksums.txt | cut -f1 -d' ') GCA_009858895.3_ASM985889v3_genomic.gff.gz" | md5sum -c -
        ```

    11. If integrity check of the files has passed (`OK`), Uncompress the ```.gz``` files
       ```
       gunzip *.gz
       ```
    12. Rename the `FASTA` and `GFF` files
        ```
        mv GCA_009858895.3_ASM985889v3_genomic.fna nCoV-2019.fasta
        mv GCA_009858895.3_ASM985889v3_genomic.gff nCoV-2019.gff
        ```