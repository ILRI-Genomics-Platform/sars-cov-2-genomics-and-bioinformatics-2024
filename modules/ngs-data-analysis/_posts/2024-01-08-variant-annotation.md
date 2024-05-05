---
title: Variant annotation
---

## Variant annotation

We will use [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/). It
annotates and predicts the effects of genetic variants on genes and proteins
(such as amino acid changes). It requires a configured SnpEff database with the
annotation or features of the genome.

1. Change to the output directory ```snpeff```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/snpeff/
    ```

2. Annotate and predict variants

    ```
    java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
        nCoV-2019 \
        -c /var/scratch/$USER/AfricaCDC_training/databases/snpeff_db/snpeff.config \
        -dataDir /var/scratch/$USER/AfricaCDC_training/databases/snpeff_db/data \
        /var/scratch/$USER/AfricaCDC_training/results/ivar/sars1.vcf.gz \
        > sars1.ann.vcf
    ```
3. Compress vcf file
    ```
    bgzip -c sars1.ann.vcf > sars1.ann.vcf.gz
    ```
4. Rename the ```summary.html``` and ```genes.txt``` file
    ```
    mv snpEff_summary.html sars1.summary.html
    mv snpEff_genes.txt sars1.genes.txt
    ```

5. Create tabix index from a sorted bgzip tab-delimited genome file
    ```
    tabix -p vcf -f sars1.ann.vcf.gz
    ```

6. Generate stats from VCF file
    ```
    bcftools stats sars1.ann.vcf.gz > sars1.stats.txt
    ```

7. Filter variants
    [SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) annotates
    genomic variants using databases, filters, and manipulates genomic annotated
    variants. Once you annotated your files using SnpEff, you can use SnpSift to
    help you filter large genomic datasets in order to find the most significant
    variants for your experiment.
    ```
      java -Xmx4g -jar /export/apps/snpeff/4.1g/SnpSift.jar \
            extractFields \
            -s "," \
            -e "." \
            sars1.ann.vcf.gz \
            CHROM POS REF ALT \
            "ANN[*].GENE" "ANN[*].GENEID" \
            "ANN[*].IMPACT" "ANN[*].EFFECT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" \
            "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
            "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \
            "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
            "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \
            "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \
            > sars1.snpsift.txt
    ```
<br>