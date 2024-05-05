---
title: Sequence data quality control and decontamination
---

## Quality assessment

[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be crirical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the results ```fastqc``` directory
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/fastqc/
    ```
2. Run ```fastqc```
    ```
    fastqc \
        -t 1 \
        -o . \
        /var/scratch/$USER/AfricaCDC_training/data/sars1_R1.fastq.gz \
        /var/scratch/$USER/AfricaCDC_training/data/sars1_R2.fastq.gz
    ```
    ***Optional***
        Run step 3. above for the other 2 samples.

<br>

## Quality and adapter filtering

The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```fastp``` directory.
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/fastp/
    ```

2. Run ```fastp```. `i,I` (input(s)) are for read1, read2; respectively. `o,O` (output(s)) are for the respective read1, read2; respectively. The `2>` construct redirects the standard error channel for saving as a log file.

    ```
    fastp \
        -w 1 \
        -i /var/scratch/$USER/AfricaCDC_training/data/sars1_R1.fastq.gz \
        -I /var/scratch/$USER/AfricaCDC_training/data/sars1_R2.fastq.gz \
        -o sars1_R1.trim.fastq.gz \
        -O sars1_R2.trim.fastq.gz \
        -h sars1.fastp.html \
        -j sars1.fastp.json \
        2> sars1.fastp.log
    ```

    ***Optional***
        Run steps 3 and 4 above for the other 2 samples.
<br>

## Decontamination
At times, sequencing experients will pick up non-target nucleic acids: for instance, host genetic material in SARS-CoV-2  sequencing. Such may obscure our signal of interest in the data; therefore, it is important to minimise or remove such sources of noise (unwanted background).
There are several bioinformatics tools and databases which can be used in querying the reads data in order to remove such noise. A commonly used tool is [Kraken2](https://github.com/DerrickWood/kraken2/wiki).
Kraken2 is a fast and memory efficient tool for taxonomic assignment of metagenomics sequencing reads. ```Kraken2``` can allow us to query the composition of our samples by searching for sequence reads against a pre-formatted database ("contaminant").

>**Note**
>Preformatted Kraken 2 and Bracken indexes can be found here: https://benlangmead.github.io/aws-indexes/k2 and downloaded without need of building new ones from scractch.


In this tutorial, we will use pre-formatted kraken2 ```human``` database to identify human-derived reads in our samples. This may give us an indication of contamination from host reads.

**Quiz:** *What type of contaminants would you think of in a SARS-CoV-2 sequencing experiment?*

---
<details close>
  <summary>Answer</summary>
  Host DNA,
  Host RNA and
  Internal control (PhiX)
</details>

---

1. Human database search
    - Change into the `kraken` directory results
        ```
        cd /var/scratch/$USER/AfricaCDC_training/results/kraken
        ```
    - Run `kraken2`
        ```
        kraken2 \
              -db /var/scratch/$USER/AfricaCDC_training/databases/kraken2-human-db \
              --threads 1 \
              --unclassified-out sars1.unclassified#.fastq \
              --classified-out sars1.classified#.fastq \
              --report sars1.kraken2.report.txt \
              --output sars1.kraken2.out \
              --gzip-compressed \
              --report-zero-counts \
              --paired /var/scratch/$USER/AfricaCDC_training/results/fastp/sars1_R1.trim.fastq.gz \
              /var/scratch/$USER/AfricaCDC_training/results/fastp/sars1_R2.trim.fastq.gz
        ```
    **Quiz:** *How many reads have hit the human genome as targets in the sample(s)?*

    ---
    <details close>
      <summary>Answer</summary>
      sars1  -  60529 reads    (13.00%)<br>
    </details>

    ---

    **Quiz:** *What percent of the sequencing reads are classified as non-human?*

    More information on output formats can be found [here](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats).

    ***Optional***
        Run steps 1 above for the other 2 samples.


    **Quiz:** *Based on the results, which sample do you think will give us good results in terms of genome coverage?*

<br>