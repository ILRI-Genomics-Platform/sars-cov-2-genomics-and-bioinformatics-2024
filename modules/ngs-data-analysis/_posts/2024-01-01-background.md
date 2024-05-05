---
title: Background
---

## Background on sequence data
We will use a dataset comprising of raw sequence reads of SARS-CoV-2 samples
obtained from a sequencing run on NextSeq 550 platorm at [ILRI](www.ilri.org).
NextSeq 550 flowcell uses 4 lanes; and so, 4 reads of data per sequenced sample
corresponding to the

- `https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars1_R1.fastq.gz` 
- `https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars1_R2.fastq.gz` 
- `https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars2_R1.fastq.gz` 
- `https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars2_R2.fastq.gz` 
- `https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars3_R1.fastq.gz` 
- `https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/sars3_R2.fastq.gz`

4 lanes are generated with suffixes L001, L002, L003 and L004. The dataset we
are using in this tutorial comprises of already concatenated sequences. These
reads can be combined/concatenated into a single file bearing in mind the type
of library sequencing either ```single``` or ```paired-end```.

<br>