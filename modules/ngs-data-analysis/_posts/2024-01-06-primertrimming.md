---
title: Primer trimming of amplicons
---

## Primer trimming
Trim amplicon primers using
[```iVar```](https://andersen-lab.github.io/ivar/html/manualpage.html) iVar uses
primer positions supplied in a BED file to soft clip primer sequences from an
aligned and sorted BAM file. Following this, the reads are trimmed based on a
quality threshold (Default: 20). To do the quality trimming, iVar uses a sliding
window approach (Default: 4). The windows slides from the 5' end to the 3' end
and if at any point the average base quality in the window falls below the
threshold, the remaining read is soft clipped. If after trimming, the length of
the read is greater than the minimum length specified (Default: 30), the read is
written to the new trimmed BAM file.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Run the command to trim primers

    ```
    ivar trim \
        -i /var/scratch/$USER/AfricaCDC_training/results/bowtie2/sars1.sorted.bam \
        -b /var/scratch/$USER/AfricaCDC_training/primer-schemes/V3/nCoV-2019.primer.bed \
        -p sars1.primertrimmed \
        -m 30 \
        -q 20 > sars1.ivar.log
    ```
3. Sort the primer trimmed alignment
    ```
    samtools sort \
          -@ 1 \
          -o sars1.primertrimmed.sorted.bam \
          -T sars1 sars1.primertrimmed.bam
    ```
4. Index the sorted primer trimmed alignment
    ```
    samtools index -@ 1 sars1.primertrimmed.sorted.bam
    ```


## Compute coverage 
Here we will use [bedtools](https://github.com/arq5x/bedtools2), your swiss-army
knife for genomic arithmetic and interval manipulation.

1. Change to the output directory ```bedtools```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/bedtools/
    ```

2. Compute coverage
    ```
    bedtools \
        genomecov \
        -d \
        -ibam \
        /var/scratch/$USER/AfricaCDC_training/results/ivar/sars1.primertrimmed.sorted.bam \
        > sars1.coverage
    ```
3. Plot to visualize

    ```
    Rscript /var/scratch/$USER/AfricaCDC_training/scripts/plotGenomecov.R sars1.coverage
    ```
<br>