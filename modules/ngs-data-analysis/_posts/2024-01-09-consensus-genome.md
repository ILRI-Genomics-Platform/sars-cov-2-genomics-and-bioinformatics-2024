---
title: Consensus genome reconstruction
---

## Consensus genome reconstruction
To generate a consensus sequence iVar uses the output of samtools mpileup
command. The mpileup output must be piped into ivar consensus. There are five
parameters that can be set:
- minimum quality ```-q``` (Default: 20).
- minimum frequency threshold ```-t``` (Default: 0).
- minimum depth to call a consensus ```-m``` (Default: 10).
- a flag ```-n``` to exclude nucleotides from regions with depth less than the
  minimum depth and a character to call in regions with coverage lower than the
  speicifed minimum depth (Default: 'N').

Minimum quality is the minimum quality of a base to be considered in
calculations of variant frequencies at a given position. Minimum frequency
threshold is the minimum frequency that a base must match to be called as the
consensus base at a position. If one base is not enough to match a given
frequency, then an ambigious nucleotide is called at that position. Minimum
depth is the minimum required depth to call a consensus. If ```-k``` flag is set
then these regions are not included in the consensus sequence. If ```-k``` is
not set then by default, a 'N' is called in these regions. You can also specfy
which character you want to add to the consensus to cover regions with depth
less than the minimum depth. This can be done using ```-n``` option. It takes
one of two values: ```-``` or ```N```.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Generate pileup and consensus genome sequences

    ```
    samtools \
            mpileup \
            --reference /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            -aa \
            /var/scratch/$USER/AfricaCDC_training/results/ivar/sars1.primertrimmed.sorted.bam \
            | tee sars1.mpileup \
            | ivar \
                consensus \
                -t 0.75 \
                -q 20 \
                -m 10 \
                -n N \
                -p sars1.cons
    ```
The ```tee``` command reads from the standard input and writes to both standard
output and one or more files at the same time. ```tee``` is mostly used in
combination with other commands through piping.

3. Rename the consensus genome header

    ```
    sed -i '/^>/s/Consensus_\(.*\)_threshold.*/\1/' sars1.cons.fa
    ```
<br>