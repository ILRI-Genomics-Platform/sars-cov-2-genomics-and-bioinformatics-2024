---
title: Sequence Alignment
---

## Sequence alignment
Aligning sequence reads to a reference genome is the first step in many
comparative genomics pipelines, including pipelines for variant calling,
isoform quantitation and differential gene expression. In many cases, the
alignment step is the slowest. This is because for each read the aligner must
solve a difficult computational problem: determining the read's likely point of
origin with respect to a reference genome. This is always non trivial for
several reasons:
- The reference genome is often very big. Searching big things is harder than
  searching small things.
- You aren’t always looking for exact matches in the reference genome–or, at
least, probably not.

Here we use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), an
ultrafast and memory-efficient tool for aligning sequencing reads to long
reference sequences.


1. Change to the ```bowtie2``` directory.
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/bowtie2
    ```

2. Run the ```bowtie2``` command to align reads to the reference genome.

    ```
    bowtie2 \
          -x /var/scratch/$USER/AfricaCDC_training/genome/bowtie2/nCoV-2019 \
          -1 /var/scratch/$USER/AfricaCDC_training/results/kraken/sars1.unclassified_1.fastq \
          -2 /var/scratch/$USER/AfricaCDC_training/results/kraken/sars1.unclassified_2.fastq \
          --threads 1 \
          --un-conc-gz sars1.unmapped.fastq.gz \
          --local \
          --very-sensitive-local \
          2> sars1.bowtie2.log \
          | samtools view -@ 1 -F4 -bhS -o sars1.trim.dec.bam -
    ```
    ***Optional***
        Run steps 1 and 2. above for the other 2 samples.
<br>

## Sort and Index alignment map
Alignments can often be manipulated using
[```samtools```](http://www.htslib.org/) using several sub commands

1. Sort the converted binary alignment (```.bam```)
    ```
    samtools sort -@ 1 -o sars1.sorted.bam -T sars1 sars1.trim.dec.bam
    ```

2. Index the sorted alignment
    ```
    samtools index -@ 1 sars1.sorted.bam
    ```
    ***Optional***
        Run steps 1 and 2. above for the other 2 samples.

<br>