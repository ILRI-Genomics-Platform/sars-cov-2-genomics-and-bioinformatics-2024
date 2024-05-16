---
title: Phylodynamics inference
published: true
---

#### The transmission dynamics of SARS-CoV-2 alpha (B.1.1.7) in Kenya

<br>

#### Summary
This exercise offers a guide on how to conduct maximum likelihood phylodynamic
inference to reconstruct the evolutionary dynamics based on a
set of virus sequences that have been isolated at different points in time
(‘heterochronous’ data) using [TreeTime](https://github.com/neherlab/treetime). 
The focus is on SARS-CoV-2 virus evolution, in particular on the emergence and
viral exchanges of the alpha (B.1.1.7) variant of SARS-CoV-2 in Kenya.

The alpha data set has NOT been officially analyzed and ongoing analyses is
underway. The aim is to obtain estimates of several evolution
parameters including the time of the most recent common ancestor and viral
exchanges between different locations in Kenya.



#### Introduction
The first step is to obtain datasets from [GISAID](https://gisaid.org/),
specifically the EpiCoV database. One can download the dataset interactively on
the website or programatically access the GISAID API using the wrapper R 
package [GISAIDR](https://github.com/Wytamma/GISAIDR). We have already retrieved
the data for this tutorial.



1. Login to the hpc and begin an interactive session on `compute06` with `1`
   core. Replace the `user` with your actual username

    ```bash
    ssh -y <user> hpc.ilri.cgiar.org
    interactive -w compute06 -c 1
    ```


2. Change into the global temporary directory, `/var/scratch` and create a
   directory named `AfricaCDC_training`. Once in the newly created directory, 
   initiate a symbolic link of the parent directory on your temporary directory

    ```bash
    cd /var/scratch

    mkdir -p $USER/AfricaCDC_training

    cd $USER/AfricaCDC_training

    ln -s /var/scratch/global/AfricaCDC_training/phylogenetic-datasets .
    ```

3. Create intermediate directories/sub-directories for different outputs

    ```bash
    mkdir -p phylodynamics/{deduplicated,alignment,models,iqtree,treetime}
    ```

4. Deduplicate sequences on composition

    ```
    module load seqkit/0.11.0

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/deduplicated

    touch B.1.1.7.duplicates.txt

    seqkit rmdup -s \
      -D  B.1.1.7.duplicates.txt < /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/Kenya/sarscov-kenya-B.1.1.7.fasta > \
      B.1.1.7.dedup.fasta
    ```

5. Perform multiple sequence alignment

    ```
    module load mafft/7.475

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment

    mafft --thread 1 \
      /var/scratch/$USER/AfricaCDC_training/phylodynamics/deduplicated/B.1.1.7.dedup.fasta > \
      B.1.1.7.align.fasta
    ```

6. Test evolutionary models

    ```
    module load modeltest/0.1.7

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/models

    modeltest-ng \
    --input /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.align.fasta \
    --datatype nt \
    -p 1 \
    --models HKY,GTR \
    -t ml \
    -o B.1.1.7.model
    ```


7. Maximum likelihood tree inference

    ```
    module load iqtree/2.2.0

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/iqtree

    iqtree2 \
    -s /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.align.fasta \
    -m TEST \
    -T 1 \
    -redo \
    -bb 1000 \
    --prefix B.1.1.7
    ```


8. Move your results in the `phylodynamics` to the `global` temporary directory.
    Dowload the generated `B.1.1.7.treefile` to your local computer and inspect
    the temporal signal in the sequences using the program `TempEST`. The first
    command is simply tranferring the directory `phylodynamics` from `compute06` to
    the global temporary directory in the head node of `hpc`. 

    The second command should be typed on your local computer terminal emulator e.g
    `WSL` or `MobaXterm`
    The `<DESTINATION>` in the second command is location where you want the directory
    `phylodynamics` to be transferred into. The second command will prompt for a
    `password`. Please use the password that you were provided with at the start of
    the course.

    ```bash
    rsync -avP --partial /var/scratch/$USER/AfricaCDC_training/phylodynamics /var/scratch/global/$USER/

    scp <USER>@hpc.ilri.cgiar.org:/var/scratch/global/$USER/phylodynamics/iqtree/B.1.1.7.treefile <DESTINATION>
    ```

- Start the TempEST application. TempEST by default will initiate the import menu.
Search for the file `B.1.1.7.treefile` that you downloaded and select it then
click `Open`

{% figure [caption:"Starting TempEST"] [class:"caption"] %}
![Powering TempEST for temporal signal assessment](/img/tempest_prompt.png)
{% endfigure %}


- Since we want to assess the temporal signal in our heterochronous sequences, we
will provide the sample collection dates which are embedded in the taxa names.
We will therefore use the option `Parse Dates` and select 
`Defined by a prefix and its order`. The taxa names are given as GISAID
accession|country|pangolin lineage|date separated by `|`. Date is the `last`
item, so we use `last`.

- We will select the option `Parse as a calendar date` with the `Date format` as
`yyyy-MM-dd` then click `OK`.

{% figure [caption:"Parsing dates"] [class:"caption"] %}
![Parsing dates](/img/parse_dates.png)
{% endfigure %}

- Select `Best-fitting root` using the default `Heuristic-residual mean squared`
  Function.

{% figure [caption:"Fitting the regression model"] [class:"caption"] %}
![Fit the regression model](/img/fit_model.png)
{% endfigure %}

- Select `Root-to-tip` tab.

{% figure [caption:"Inspect root-to-tip divergence vs Time"] [class:"caption"] %}
![Inspect root-to-tip divergence vs Time](/img/inspect_root_to_tip.png)
{% endfigure %}

>**Quiz**
>1. When was the probable date that B.1.1.7 variant occur in Kenya? 
>2. What is the substitution rate ob B.1.1.7 variant in this ML inference? 
>3. Are there potential outliers in the dataset?

- To remove potential outliers, we can inspect the `Residuals` tab and identify
  any sequences that deviate more than 0.0001 from the residual mean.

{% figure [caption:"Select potential outliers"] [class:"caption"] %}
![Select potential outlier](/img/highlight_outliers.png)
{% endfigure %}

- Potential outliers have been identified and written to a file named `b.1.1.7.potential.outliers.txt`. 

9. We can remove these potential outliers from the tree and repeat the ML tree
  inference step

    ```bash
    grep ">" /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.align.fasta | sed 's/>//g' > /var/scratch/$USER/AfricaCDC_training/phylodynamics/deduplicated/B.1.1.7.seqids.txt

    awk 'NR==FNR{a[$0]=1;next}!a[$0]' /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/b.1.1.7.outliers.txt /var/scratch/$USER/AfricaCDC_training/phylodynamics/deduplicated/B.1.1.7.seqids.txt > /var/scratch/$USER/AfricaCDC_training/phylodynamics/deduplicated/B.1.1.7.final.ids.txt

    module load seqtk/1.3

    seqtk subseq \
    /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.align.fasta \
    /var/scratch/$USER/AfricaCDC_training/phylodynamics/deduplicated/B.1.1.7.final.ids.txt > \
    /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.filtered.fasta
    ```

10. Re-run the  multiple sequence alignment 

    ```bash
    mafft \
    --thread 1 \
    /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.filtered.fasta > \
    /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.filtered.align.fasta
    ```

11. Re-run the ML tree building

    ```bash
    module load iqtree/2.2.0

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/iqtree

    iqtree2 \
    -s /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.filtered.align.fasta \
    -m TEST \
    -T 1 \
    -redo \
    -bb 1000 \
    --prefix B.1.1.7.filtered
    ```

12. Copy the ouput to your global temporary directory and download the new tree
    to your local computer

    ```bash
    rsync -avP --partial /var/scratch/$USER/AfricaCDC_training/phylodynamics /var/scratch/global/$USER/

    scp <USER>@hpc.ilri.cgiar.org:/var/scratch/global/$USER/phylodynamics/iqtree/B.1.1.7.treefile <DESTINATION>
    ```

13. Inspect the new tree in TempEST.

>**Quiz**
>1. Is there any improvement to the temporal signal for molecular clock
>   analysis? 

<details markdown="1">
<summary>Hint</summary>
Look at the results of the Heuristic residual mean squared function
</details>

14. Time-scaled phylogenetic inference

    ```bash

    module load treetime/0.11.3

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/treetime

    treetime \
    --tree /var/scratch/$USER/AfricaCDC_training/phylodynamics/iqtree/B.1.1.7.filtered.treefile \
    --aln /var/scratch/$USER/AfricaCDC_training/phylodynamics/alignment/B.1.1.7.filtered.align.fasta \
    --dates /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/Kenya/sarscov-kenya-B.1.1.7.csv \
    --outdir .
    ```

15. Basic viral dispersal analysis

    A migration model can be fitted to the time calibrated tree topologies
    in TreeTime, mapping the location of sampled sequences to the external tips of
    the tree. 
    The mugration model of TreeTime also infer the most likely location for internal
    nodes in the tree. 
    This can be achieved by counting the number of state changes by iterating over each
    phylogeny from the root to the external tips. State changes are counted when an
    internal node transitions from one location to a different location in the
    resulting child-node or tip(s). The timing of transition events is then recorded
    which serve as the estimated import or export event.

    ```bash
    module load R/4.2

    cd /var/scratch/$USER/AfricaCDC_training/phylodynamics/treetime

    treetime mugration \
      --tree /var/scratch/$USER/AfricaCDC_training/phylodynamics/iqtree/B.1.1.7.filtered.treefile \
      --states /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/Kenya/sarscov-kenya-B.1.1.7.csv \
      --attribute location \
      --confidence \
      --outdir .


    python /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/ancestral_locations.py \
        --tree /var/scratch/$USER/AfricaCDC_training/phylodynamics/treetime/annotated_tree.nexus \
        --mrsd 2021.893150 \
        --prefix B.1.1.7 \
        --outDir .

    Rscript \
        /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/plotDiscreteLocations.R \
        --input b.1.1.7-kenya.annotated_tree_events.csv \
        --prefix b.1.1.7


    Rscript \
        /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/plotTree.R \
        --tree /var/scratch/$USER/AfricaCDC_training/phylodynamics/treetime/timetree.nexus \
        --metadata /var/scratch/$USER/AfricaCDC_training/phylogenetic-datasets/Kenya/sarscov-kenya-B.1.1.7.csv \
        --prefix b.1.1.7
    ```

16. Transfer the `b.1.1.7.viral-exchanges.pdf` and the `b.1.1.7.time-calibrated-tree.pdf` to your global temporary
    directory and then download to your local computer

    ```bash
    rsync -avP --partial /var/scratch/$USER/AfricaCDC_training/phylodynamics /var/scratch/global/$USER/

    scp <USER>@hpc.ilri.cgiar.org:/var/scratch/global/$USER/phylodynamics/treetime/*.pdf <DESTINATION>
    ```
>**Quiz**
>1. Which cities/towns had the highest exportation events of B.1.1.7 variant?
>2. Which town/city was the biggest exporter of B.1.1.7 variant in Kenya? 