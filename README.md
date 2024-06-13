# User Guide

## Before you start

* Set up SSH key between GitHub and your JHPCE account. (https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
* Install PlasmaTools to your R environment:
```
#!/bin/bash
module load conda_R/4.3
git clone git@github.com:cancer-genomics/PlasmaTools.git
R CMD build PlasmaTools
R CMD INSTALL PlasmaTools_0.0.11.tar.gz
```

## Setting up your environment

* Log in to JHPCE
* Start a qrsh session. 
* Load git module (`module load git/2.28.0`)
* Create an empty directory with a descriptive name relating to the set of fastq files you wish to process (e.g. `danish-batch01-wgs-plasma-2021-06-02`). This will be your analysis directory. Navigate to that directory.
* From within your analysis directory, run the following line of code: `git clone git@github.com:cancer-genomics/delfi_pipeline.git`. This will download the delfi pipeline from GitHub into a folder called `delfi_pipeline/`
* Create an empty folder called `fastq`. Populate this folder with symbolic links to all the fastq files you wish to process. IMPORTANT: The Delfi pipeline makes stringent assumptions about how fastq files are named, so this is the time to format them. For each sample ID, the symbolic links should follow the following pattern: `${sample_id}_R1.fastq.gz  ${sample_id}_R2.fastq.gz`
* Find the file: `delfi_pipeline/delfi/samples.txt`. Populate this file with all the sample IDs corresponding to the fastq files. For each pair of fastq files, there should be a single sample ID. There should be one sample ID per line, and there should be no duplicates in the file.
* For clarity, the table below shows the relationship between the contents of the `fastq` directory and the contents of the `samples.txt` file.

| `fastq` Directory     | `samples.txt` file |
| :---:        |    :----:   |
| CGPLH1000P_R1.fastq.gz<br/>CGPLH1000P_R2.fastq.gz      | CGPLH1000P       |
| CGPLH1001P_R1.fastq.gz<br/>CGPLH1001P_R2.fastq.gz   | CGPLH1001P        |
| CGPLH1002P_R1.fastq.gz<br/>CGPLH1002P_R2.fastq.gz |  CGPLH1002P |
| CGPLH1003P_R1.fastq.gz<br/>CGPLH1003P_R2.fastq.gz |  CGPLH1003P |

At this point, your analysis directory should follow this structure:

```
analysis_directory_2023-01-31/
│
├── delfi_pipeline/
│   ├── alignment/
│   ├── aneuploidy/
│   ├── delfi/
│   ├── prediction/
│   ├── preprocess/
│   └── references/
│       └── hg19/
└── fastq/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
```

## Running the Delfi pipeline

* Within your analysis directory, navigate to `delfi_pipeline/delfi/`. 
* Double-check that the parameters defined in RUN_DELFI.sh are correct (genome build, sequencer, bin size, etc.)
* `./RUN_DELFI.sh`


# References

Cristiano S, Leal A, Phallen J, Fiksel J, Adleff V, Bruhm DC, Jensen SØ, Medina JE, Hruban C, White JR, Palsgrove DN, Niknafs N, Anagnostou V, Forde P, Naidoo J, Marrone K, Brahmer J, Woodward BD, Husain H, van Rooijen KL, Ørntoft MW, Madsen AH, van de Velde CJH, Verheij M, Cats A, Punt CJA, Vink GR, van Grieken NCT, Koopman M, Fijneman RJA, Johansen JS, Nielsen HJ, Meijer GA, Andersen CL, Scharpf RB, Velculescu VE. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019 Jun;570(7761):385-389. doi: 10.1038/s41586-019-1272-6. Epub 2019 May 29. PMID: 31142840

Mathios D, Johansen JS, Cristiano S, Medina JE, Phallen J, Larsen KR, Bruhm DC, Niknafs N, Ferreira L, Adleff V, Chiao JY, Leal A, Noe M, White JR, Arun AS, Hruban C, Annapragada AV, Jensen SØ, Ørntoft MW, Madsen AH, Carvalho B, de Wit M, Carey J, Dracopoli NC, Maddala T, Fang KC, Hartman AR, Forde PM, Anagnostou V, Brahmer JR, Fijneman RJA, Nielsen HJ, Meijer GA, Andersen CL, Mellemgaard A, Bojesen SE, Scharpf RB, Velculescu VE. Detection and characterization of lung cancer using cell-free DNA fragmentomes. Nat Commun. 2021 Aug 20;12(1):5060. doi: 10.1038/s41467-021-24994-w. PMID: 34417454
