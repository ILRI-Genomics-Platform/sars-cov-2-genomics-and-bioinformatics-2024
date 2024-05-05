---
title: Variant calling and the VCF file format
---

## Variant calling
`iVar` uses the output of the ```samtools mpileup``` command to call variants -
single nucleotide variants(SNVs) and indels. Pileup format consists of
TAB-separated lines, with each line representing the pileup of reads at a single
genomic position.

Several columns contain numeric quality values encoded as individual ASCII
characters. Each character can range from "!" to "~" and is decoded by taking
its ASCII value and subtracting 33; e.g., "A" encodes the numeric value 32.

The first three columns give the position and reference:

1. Chromosome name.
2. 1-based position on the chromosome.
3. Reference base at this position (this will be "N" on all lines if ```-f``` or
   ```--fasta-ref``` has not been used)

In generating the mpileup, we will use the flags:
```--count-orphans```: Do not skip anomalous read pairs in variant calling.
Anomalous read pairs are those marked in the FLAG field as paired in sequencing
but without the properly-paired flag set.
```--ignore-overlaps```: Disable read-pair overlap detection
```--no-BAQ```: Disable base alignment quality (BAQ) computation

The ```tee``` command, used with a pipe, reads standard input from ```samtools
mpileup```, then writes the output of the program to standard output and
simultaneously copies it into the specified file ```.mpileup```


In order to call variants correctly, the reference file used for alignment must
be passed to `iVar` using the ```-r``` flag. The output of samtools pileup is
piped into ivar variants to generate a ```.tsv``` file with the variants. There
are two parameters that can be set for variant calling using `iVar` - minimum
quality (Default: 20) and minimum frequency (Default: 0.03). Minimum quality is
the minimum quality for a base to be counted towards the ungapped depth to
calculate iSNV frequency at a given position. For insertions, the quality metric
is discarded and the mpileup depth is used directly. Minimum frequency is the
minimum frequency required for a SNV or indel to be reported.

`iVar` can identify codons and translate variants into amino acids using a
[GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
format containing the required coding regions (CDS). In absence of a GFF file,
iVar will not perform the translation and "NA" will be added to the output file
in place of the reference and alternate codons and amino acids.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Call variants

    ```
    samtools mpileup \
            --ignore-overlaps \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            --reference /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
            /var/scratch/$USER/AfricaCDC_training/results/ivar/sars1.primertrimmed.sorted.bam \
            | tee sars1.mpileup \
            | ivar \
                variants \
                -t 0.25 \
                -q 20 \
                -m 10 \
                -g /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.gff \
                -r /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
                -p sars1.variants
    ```
3. Convert the variants from ```.tsv``` to ```.vcf``` (Variant Call Format)

    ```
    python /var/scratch/$USER/AfricaCDC_training/scripts/ivar_variants_to_vcf.py \
      sars1.variants.tsv \
      sars1.vcf \
      --pass_only \
      --allele_freq_thresh 0.75 > sars1.variant.counts.log
    ```
    ## VCF file format

    The header begins the file and provides metadata describing the body     of the file.
    Header lines are denoted as starting with `#`.
    Special keywords in the header are denoted with `##`.
    Recommended keywords   include fileformat, fileDate and reference.

    The header contains keywords that optionally semantically and syntactically describe the fields used in the body of the file, notably INFO, FILTER, and FORMAT.


    |   |      Name    |  Brief description (see the specification for details)[VCF](https://samtools.github.io/hts-specs/VCFv4.1.pdf).  |
    |---|:-------------|:---------------------------------------------------------|
    | 1 |  CHROM       |The name of the sequence (typically a chromosome) on which the variation is being called.                                                           |
    | 2 |  POS         |The 1-based position of the variation on the given sequence.                                                          |
    | 3 |  ID          |The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.                                                           |
    | 4 |  REF         |The reference base (or bases in the case of an indel) at the given position on the given reference sequence.                                                          |
    | 5 |  ALT         |The list of alternative alleles at this position.                                                           |
    | 6 |  QUAL        |A quality score associated with the inference of the given alleles.                                                          |
    | 7 |  FILTER      |A flag indicating which of a given set of filters the variation has failed or PASS if all the filters were passed successfully.                                                          |
    | 8 |  INFO        |An extensible list of key-value pairs (fields) describing the variation.                                                          |
    | 9 |  FORMAT      |An (optional) extensible list of fields for describing the samples                                                          |
    | + |  SAMPLES     |For each (optional) sample described in the file, values are given for the fields listed in FORMAT                                                           |

4. Compress vcf file
    ```
    bgzip -c sars1.vcf > sars1.vcf.gz
    ```

5. Create tabix index from a sorted bgzip tab-delimited genome file
    ```
    tabix -p vcf -f sars1.vcf.gz
    ```
6. Generate stats from VCF file
    ```
    bcftools stats sars1.vcf.gz > sars1.stats.txt
    ```
<br>