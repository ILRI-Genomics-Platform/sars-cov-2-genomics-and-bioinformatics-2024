---
title: Data and metadata management and sharing
---

## Data Retrieval and Review
Having sequenced our samples in both MiSeq (Illumina) and MinION (ONT) we can now transfer the sequence output to the HPC where we will conduct the bioinformatics analysis.

### Transfer of data: Illumina
MiSeq is based on Windows and data transfer will be done by copy-pasting the data to HPC through the network. Ensure the transfer is completed successfully without errors.

### Transfer of data: MinION
The MinION sequencer stores its sequencing output in a Linux based computer. To transfer the data we logged into computer and transferred the data on the command line as follows.
```
rsync -avP <path-to-the-directory-with_sequencing-ouput>/ <username>:<path-to-the-directory-to-store-sequencing-ouput>/
```
Replaced `<path-to-the-directory-with_sequencing-ouput>` with the path to the directory storing the sequencing output. Replaced `<HPC-login-username>` with your hpc login username (i.e user##@hpc.ilri.cgiar.org) and `<path-to-the-directory-to-store-sequencing-output>` with path to the directory you want to store the data in the HPC.
Example:`
rsync -avP /media/SeqData_LTS/20220405_1121_MN2816_FAH91436_2b8e9827/ <username>@hpc.ilri.cgiar.org:/var/scratch/global/$USER/20220405_1121_MN2816_FAH91436_2b8e9827`

### Reviewing data: Illumina
Change working directory into the directory that stores the Illumina dataset
```
cd /var/scratch/global/miseq/
```
Now let us view what is the sequencing output. The output is a FASTA format. Note: in the second command replace `<one-of-the-fastq.gz>`
 with the name of the fastq.gz file
 ```
ls fastq/
less -S fastq/<one-of-the-fastq.gz>
```

### Reviewing data: ONT
Change working directory into directory that stores the ONT data. Note: replace `<name-of-run-folder>` with the name of the run folder
```
cd /var/scratch/global/ONT/
ls <name-of-run-folder>
```

## Working with metadata
Metadata is the data associated with your main data; it describes your data in a manner that can allow drawing information upon analysing the data. Often, the importance of `metadata` is ignored; things like how to capture, store, encode and organise metadata, especially, having the downstream analyses and interpretation in mind.

> **IMPORTANT**
> Data without metadata is, mostly, garbage.

We will take a look at an example metadata to highlight some concerns with metadata, and the reason(s) why they are important in data analyses workflows.

test_lab|case_id|lab_id|loc|age|gender|occup|samp_type|symp|vacc_state|coll_dt|confir_dt|recep_dt
|---|---|---|---|---|---|---|---|---|---|---|---|---|
<strong>TESTING_LAB</strong>|<strong>CASE-ID</strong>|<strong>SampleNumber</strong>|<strong>LOCATION</strong>|<strong>AGE</strong>|<strong>GENDER (M/F)</strong>|<strong>OCCUPATION</strong>|<strong>SAMPLE TYPE</strong>|<strong>SYMPTOMS SHOWN (COUGH;FEVER;ETC)</strong>|<strong>VACCINATION STATUS</strong>|<strong>DATE OF SAMPLE COLLECTION</strong>|<strong>LAB CONFIRMATION DATE</strong>|<strong>DATE SAMPLE RECEIVED IN THE LAB</strong>
LABA|place/id/date|COVD0308|Some Place|80|F|Food handler|NP & OP Swab|Asymptomatic|Yes|5th June 2020|08-Jun-20|08/06/20
LABB|COM/SARS001/2022|COVD360|Comoros|38|F|None|NP-OP Swab|FC;CO;H|Yes|5th June 2020|08-Jun-20|08/06/20
LABC|DJI/SARS001/2023|COVD 273|Djibouti|30|M|Business|NP-OP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABD|SWZ/SARS001/2024|COVD154|Eswatini|23|Female|Food seller|OP-NP Swap|No symptoms|Yes|5th June 2020|08-Jun-20|08/06/20
LABE|ETH/SARS001/2025|COVD0875|Ethiopia|34|M|Targeted testing|NP-OP Swab|Asymptomatic|Yes|5th June 2020|08-Jun-20|08/06/20
LABF|LBY/SARS001/2027|COVD00672|Libya|18|Female|Food handler|NP-OP Swab|Fever/Chills, Cough, Headache|Yes|5th June 2020|08-Jun-20|08/06/20
LABG|MDG/SARS001/2028|COVD499|Madagascar|25|F|Targeted testing|NP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABH|MUS/SARS001/2029|COVD078|Mauritius|25|F|Not indicated|NP-OP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABI|SYC/SARS001/2030|COVD579|Seychelles|2 years|Male|Targeted testing|NP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABJ|SOM/SARS001/2031|300|Somalia|26|F|Food handler|NP-OP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABK|SSD/SARS001/2031|COVD00381|South Sudan|45|M|NA|NP-OP Swab|Asymptomatic|Yes|5th June 2020|08-Jun-20|08/06/20

What are some of the issues you notice with the above metadata?

---
<details close>
  <summary>Answer</summary>

  - Inconsistent header naming
  - Future dates
  - Mixed data types within column
  - Inconsistent date formats
  - Inconsistent sample type capturing
  - Inconsistent symptoms
  - Spaces in headers
  - Use of commas
  - ...
</details>


---

***Computation-wise, whichever way (poor/good) you choose to organise your data, ensure consistency***.

<br>