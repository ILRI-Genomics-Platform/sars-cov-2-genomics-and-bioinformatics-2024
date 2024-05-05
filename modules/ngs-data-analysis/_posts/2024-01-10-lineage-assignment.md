---
title: Clade and lineage assignment
---

## Clade assignment
[**Nextclade**](https://docs.nextstrain.org/projects/nextclade/en/stable/) is a
tool within the [**Nextrain**](https://nextstrain.org/) collection that uses
sequence differences for their assignment to
[clades](https://clades.nextstrain.org/). It also reports suspect quality issues
with such sequences. There are both [web-](https://clades.nextstrain.org/) and
[command-line-interfaces](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html)
for *nextclade*. To run it in the command-line, we need some reference files:
genome, feature map, origin tree, primers and quality configurations. Luckily,
for SARS-CoV-2, these can be easily retrieved using the same tool, otherwise,
you will have to create/retrieve accordingly.

1. Get the reference dataset
    ```
    nextclade dataset get --name 'sars-cov-2' --reference 'MN908947' --output-dir /var/scratch/$USER/AfricaCDC_training/nextclade_db
    ```
2. Perform clade assignment
    ```
    nextclade \
       --input-fasta /var/scratch/$USER/AfricaCDC_training/results/ivar/sars1.cons.fa \
       --input-dataset /var/scratch/$USER/AfricaCDC_training/nextclade_db \
       --input-root-seq /var/scratch/$USER/AfricaCDC_training/nextclade_db/reference.fasta \
       --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
       --input-gene-map /var/scratch/$USER/AfricaCDC_training/nextclade_db/genemap.gff \
       --input-tree /var/scratch/$USER/AfricaCDC_training/nextclade_db/tree.json \
       --input-qc-config /var/scratch/$USER/AfricaCDC_training/nextclade_db/qc.json \
       --input-pcr-primers /var/scratch/$USER/AfricaCDC_training/nextclade_db/primers.csv \
       --output-csv /var/scratch/$USER/AfricaCDC_training/results/nextclade/sars1.csv \
       --output-tree /var/scratch/$USER/AfricaCDC_training/results/nextclade/sars1.auspice.json \
       --output-dir /var/scratch/$USER/AfricaCDC_training/results/nextclade/ \
       --output-basename sars1.cons \
       2> /var/scratch/$USER/AfricaCDC_training/results/nextclade/sars1.nextclade.log
    ```

3. Visualization: The output of Nextclade includes a phylogenetic tree in
   `.json` format. This tree can be visualized in
   [Auspice](https://auspice.us/). First let us download the the `.json` file:
```
cp /var/scratch/$USER/AfricaCDC_training/results/nextclade/sars1.auspice.json ~/
```
In your local computer use `scp` to copy the file to any desired destination:
```
scp <user_name>@hpc.ilri.cgiar.org:/home/<user_name>/*.json <destination_folder>
```
Open [Auspice](https://auspice.us/) and drag and drop the `.json` file in the [Auspice](https://auspice.us/). Now edit the tree.
  - In `Dataset` click the drop down arrow and select `ncov`, below it select `open` and below it select `global`.
  - In `Color By` click the drop down arrow and select `clade`.
  - Do any other adjustments as you wish.


## Lineage Assignment

[Phylogenetic Assignment of Named Global Outbreak Lineages
(Pangolin)](https://cov-lineages.org/resources/pangolin.html) implements a
dynamic nomenclature of SARS-CoV-2 lineages, known as the Pango nomenclature. To
assign [Pangolin Lineages](https://cov-lineages.org/lineage_list.html), we will
use [the web version of pangolin](https://pangolin.cog-uk.io/). It also has a
robust [command-line version](https://github.com/cov-lineages/pangolin) that we
will look into later. With the web version we must retrieve out consensus genome
from the analysis server (HPC) to our local computer. Follow the following steps
to assign your query sequences pangolin lineages. Also here is a
[***tutorial***](https://cov-lineages.org/resources/pangolin/tutorial.html).

1. In your bowser open [Pangolin Web Application](https://pangolin.cog-uk.io/).
   This is the online version of
   [Pangolin](https://github.com/cov-lineages/pangolin).
2. Now copy-paste/drag-drop your consensus file to the site and click `Start Analysis`
3. Once done, download the results and if you need to, copy the names of
   sequences that failed the analysis to a file. The download is called
   `results.csv`.

Here is how pangolin performs the analysis:
![alt text](https://cov-lineages.org/assets/images/pangolin_pipelines.svg "Pangolin Analysis Workflow")

Alternatively to perform the commandline analysis for Pangolin, let us proceed
as follows. We will need to use a singularity image (think of a singularity
image as a ready-to-use container within which we have packaged all the software
needed to do a certain task) in this case packaging pangolin softwares*.
1. Let us create a directory to store our image:
```
mkdir /var/scratch/$USER/AfricaCDC_training/singularity
```
2. Download the image:
```
singularity pull --dir /var/scratch/$USER/AfricaCDC_training/singularity/ \
                    --force docker://staphb/pangolin:latest
```
3. Download Pangolin's referemce data: Downloads/updates to latest release of
   pangoLEARN and constellations
```
singularity run /var/scratch/$USER/AfricaCDC_training/singularity/pangolin_latest.sif \
          pangolin --update-data \
          --datadir /var/scratch/$USER/AfricaCDC_training/pangolin_db
```
4. Conduct Pangolin Lineage assignment:
```
singularity run /var/scratch/$USER/AfricaCDC_training/singularity/pangolin_latest.sif
          pangolin /var/scratch/$USER/AfricaCDC_training/results/ivar/sars1.cons.fa \
          --alignment \
          --usher \
          --max-ambig 0.3 \
          --min-length 25000 \
          --outdir /var/scratch/$USER/AfricaCDC_training/results/pangolin/ \
          --outfile sars1.pangolin.usher.csv \
          --datadir /var/scratch/$USER/AfricaCDC_training/pangolin_db/ \
          2> /var/scratch/$USER/AfricaCDC_training/results/pangolin/sars1.pangolin.usher.log
```
<br>