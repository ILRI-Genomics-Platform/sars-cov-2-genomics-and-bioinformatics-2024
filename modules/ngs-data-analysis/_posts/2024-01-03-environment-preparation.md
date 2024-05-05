---
title: Environment preparation
---

## Environment preparation

### Loading modules
1. Clear the environment.
    ```
    module purge
    ```
2. Load modules using the `module load <tool-name>`command.
    ```
    module load fastqc/0.11.7
    module load fastp/0.22.0
    module load kraken/2.0.8-beta
    module load bowtie2/2.3.4.1
    module load samtools/1.11
    module load ivar/1.3.1
    module load bedtools/2.29.0
    module load R/4.2
    module load bcftools/1.11
    module load snpeff/4.1g
    module load multiqc/1.12
    module load nextclade/1.11.0
    module load python/3.9
    ```

    **Optional**
    The above modules can also be loaded using a single command
    ```
    module load fastqc/0.11.7 fastp/0.22.0 \
    kraken/2.0.8-beta bowtie2/2.3.4.1 samtools/1.11 ivar/1.3.1 \
    bedtools/2.29.0 R/4.2 bcftools/1.11 snpeff/4.1g multiqc/1.12 \
    nextclade/1.11.0 python/3.9
    ```
3. To list the loaded modules, type the below command.
    ```
    module list
    ```

## Prepare the reference genome


1. While still in the `genome` directory, we will index the reference sequence using samtools' `faidx`. Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname, seqlength, first-base offset, seqlinewidth` without `\n` (newline character) and `seqlinewidth` with`\n`. This is essential for samtools' operations.

    ```
    samtools faidx nCoV-2019.fasta
    ```
    The above command generates the index for reference genome with the name `nCoV-2019.fasta.fai`.
2. We can take a sneak-view of the generated file and manipulate it for fun, say, to extract the genome size of reference fasta. This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest.
    ```
    cut -f 1,2 nCoV-2019.fasta.fai > nCoV-2019.fasta.sizes
    ```
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bowtie2-build``` command.

    ```
    mkdir /var/scratch/$USER/AfricaCDC_training/genome/bowtie2
    ```

    ```
    bowtie2-build \
      --threads 1 \
      /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
      /var/scratch/$USER/AfricaCDC_training/genome/bowtie2/nCoV-2019
    ```
    The above command generates index files with the suffix `.bt2` for the reference genome with the prefix `nCoV-2019.`
4. Build SnpEff database for the reference genome

    [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/), a variant annotation and predictor needs a database to perform genomic annotations. There are pre-built databases for thousands of genomes, so chances are that your organism of choice already has a SnpEff database available.

    >**Note** We will use pre-built SARS-CoV-2 SnpEff database


    ***Optional***
    In the (unlikely?) event that you need to build one yourself, you can build one using the commands found [here](http://pcingola.github.io/SnpEff/snpeff/build_db/)