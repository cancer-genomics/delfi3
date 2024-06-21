# User Guide

## Overview
This pipeline generates DELFI `features` from paired-end `fastq` files. delfi3 contains scripts in `align` that align the fastq files using Bowtie2 (saved in `cram`) after trimming and filtering by fastp. After alignment,
the duplicates are marked and removed and stored in `bed`. GC counts, CpG types and counts are computed are added to `bed` and saved in `starch` and fastq, alignment and fragment QC are stored in `qc`. Fragment counts per
bin for the selected bin size are computed after GC-correction and removal of black-listed regions in `bins`. Finally, `features` like z-scores, frament ratios and binned coverages are saved.


## Obtaining Delfi3
This pipeline is meant to be run on a cluster. The scripts make use of SLURM for submitting jobs.
* Log in to the cluster
* Start an interactive session on a compute node 
* Create an empty project directory. Navigate to that directory.
* From within your project directory, run the following line of code: `git clone git@github.com:cancer-genomics/delfi3.git`. This will download the delfi3 pipeline from GitHub into a folder called `delfi3/`
* Create an empty folder called `fastq` in the project directory. Populate this folder with fastq files (or symbolic links to fastq files) you wish to process. IMPORTANT: The Delfi pipeline makes stringent assumptions about how fastq files are named, so this is the time to format them. For each sample ID, the symbolic links should follow the following pattern: `${sample_id}_R1.fastq.gz  ${sample_id}_R2.fastq.gz`
* Find the file: `delfi3/delfi/samples.txt`. Populate this file with all the sample IDs corresponding to the fastq files. For each pair of fastq files, there should be a single sample ID. There should be one sample ID per line, and there should be no duplicates in the file.
* For clarity, the table below shows the relationship between the contents of the `fastq` directory and the contents of the `samples.txt` file.

| `fastq` Directory     | `samples.txt` file |
| :---:        |    :----:   |
| CGPLH1000P_R1.fastq.gz<br/>CGPLH1000P_R2.fastq.gz      | CGPLH1000P       |
| CGPLH1001P_R1.fastq.gz<br/>CGPLH1001P_R2.fastq.gz   | CGPLH1001P        |
| CGPLH1002P_R1.fastq.gz<br/>CGPLH1002P_R2.fastq.gz |  CGPLH1002P |
| CGPLH1003P_R1.fastq.gz<br/>CGPLH1003P_R2.fastq.gz |  CGPLH1003P |

At this point, your project directory should follow this structure:

```
project_directory_2023-01-31/
├── delfi3
│   ├── alignment
│   ├── delfi
│   ├── LICENSE
│   ├── README.md
│   ├── references
│   └── VERSION
└── fastq
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
```

Software to be installed to run this pipeline
* fastp
* Bowtie2
* Samtools 
* Bedtools
* Bedops
* R 4.3
Note that, in the scripts, Bowtie2, Samtools, Bedtools and R are loaded into the environment using `module load` command.

You will need to update `/delfi/RUN_DELFI.sh` with your local paths to `Bedops` and `fastp`.

## Reference Files

The reference files are not included in the repo because of large file size. Save these files in the `/references/hg19/` directory.

Download the hg19 reference genome from UCSC Downloads Feb. 2009 (GRCh37/hg19) sequence files. Use bowtie2-build to create index of downloaded fasta file for alignment.

Under `references/genome/`, it is necessary to add `cytosine_ref.hg19.starch`, a starch format file of the location of all cytosine nucleotides in the hg19 reference genome.


## Running the Delfi pipeline

* Within your project directory, navigate to `delfi3/delfi/`. 
* Double-check that the parameters defined in RUN_DELFI.sh are correct (genome build, sequencer, bin size, etc.)
* `./RUN_DELFI.sh`


## Publications

Cristiano S, Leal A, Phallen J, Fiksel J, Adleff V, Bruhm DC, Jensen SØ, Medina JE, Hruban C, White JR, Palsgrove DN, Niknafs N, Anagnostou V, Forde P, Naidoo J, Marrone K, Brahmer J, Woodward BD, Husain H, van Rooijen KL, Ørntoft MW, Madsen AH, van de Velde CJH, Verheij M, Cats A, Punt CJA, Vink GR, van Grieken NCT, Koopman M, Fijneman RJA, Johansen JS, Nielsen HJ, Meijer GA, Andersen CL, Scharpf RB, Velculescu VE. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019 Jun;570(7761):385-389. doi: 10.1038/s41586-019-1272-6. Epub 2019 May 29. PMID: 31142840

Mathios D, Johansen JS, Cristiano S, Medina JE, Phallen J, Larsen KR, Bruhm DC, Niknafs N, Ferreira L, Adleff V, Chiao JY, Leal A, Noe M, White JR, Arun AS, Hruban C, Annapragada AV, Jensen SØ, Ørntoft MW, Madsen AH, Carvalho B, de Wit M, Carey J, Dracopoli NC, Maddala T, Fang KC, Hartman AR, Forde PM, Anagnostou V, Brahmer JR, Fijneman RJA, Nielsen HJ, Meijer GA, Andersen CL, Mellemgaard A, Bojesen SE, Scharpf RB, Velculescu VE. Detection and characterization of lung cancer using cell-free DNA fragmentomes. Nat Commun. 2021 Aug 20;12(1):5060. doi: 10.1038/s41467-021-24994-w. PMID: 34417454
