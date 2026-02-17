# Phylogenetic tree pipeline
```
ramisha-azhar-phylogenetic-tree/
├── README.md
├── code.py
│
├──fastqc/
|
├──multiqc/
│   └── multiqc_data/
│       ├── fastqc-status-check-heatmap.txt
│       ├── fastqc_adapter_content_plot.txt
│       ├── fastqc_per_base_n_content_plot.txt
│       ├── fastqc_per_base_sequence_quality_plot.txt
│       ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│       ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│       ├── fastqc_per_sequence_quality_scores_plot.txt
│       ├── fastqc_sequence_counts_plot.txt
│       ├── fastqc_sequence_duplication_levels_plot.txt
│       ├── fastqc_sequence_length_distribution_plot.txt
│       ├── multiqc_citations.txt
│       ├── multiqc_fastqc.txt
│       ├── multiqc_general_stats.txt
│       ├── multiqc_software_versions.txt
│       └── multiqc_sources.txt
│
├──trimmed_reads/
│   ├── ERR4079133_1.fastq.gz_trimming_report.txt
│   ├── ERR4079133_1_val_1.fq.gz
│   ├── ERR4079133_2.fastq.gz_trimming_report.txt
│   ├── ERR4079133_2_val_2.fq.gz
│   ├── ERR4079285_1.fastq.gz_trimming_report.txt
│   ├── ERR4079285_1_val_1.fq.gz
│   ├── ERR4079285_2.fastq.gz_trimming_report.txt
│   └── ERR4079285_2_val_2.fq.gz
│
├── mash_results/
│       ├── ERR4079133_1_val_1.fq.gz.msh
│       ├── ERR4079133_1_val_1.fq.gz_vs_ref.txt
│       ├── ERR4079133_2_val_2.fq.gz.msh
│       ├── ERR4079133_2_val_2.fq.gz_vs_ref.txt
│       ├── ERR4079285_1_val_1.fq.gz.msh
│       ├── ERR4079285_1_val_1.fq.gz_vs_ref.txt
│       ├── ERR4079285_2_val_2.fq.gz.msh
│       ├── ERR4079285_2_val_2.fq.gz_vs_ref.txt
│       └── Fv10027Complete.fasta.msh
│
├──alignment/
│   ├── ERR4079133.sam
│   ├── ERR4079133.bam
│   ├── ERR4079133.sorted.bam
│   ├── ERR4079133_dedup.bam
│   ├── ERR4079285_sorted.bam
│   ├── ERR4079285_dedup.bam
│   ├── Fv10027Complete.fasta.amb
│   ├── Fv10027Complete.fasta.ann
│   └── alignment_stats/
│       ├── ERR4079133_dedup_metrics.txt
│       ├── ERR4079133_flagstat.txt
│       ├── ERR4079285_dedup_metrics.txt
│       └── ERR4079285_flagstat.txt
│
├── variant_calling/
│   ├── ERR4079133.g.vcf.gz
│   ├── ERR4079133.g.vcf.gz.csi
│   ├── ERR4079133.vcf.gz
│   ├── ERR4079133.vcf.gz.csi
│   ├── ERR4079285.g.vcf.gz
│   ├── ERR4079285.g.vcf.gz.csi
│   ├── ERR4079285.vcf.gz
│   ├── ERR4079285.vcf.gz.csi
│   ├── merged.vcf.gz
│   └── merged.vcf.gz.csi
│
├── 07_phylogenetics/
│   ├── distance_matrix.txt
│   ├── pairwise_stats.txt
│   └── two_taxon_tree.nwk
│
└── raw_data/
    └── fasta_files/
        └── data/
            ├── ERR4079133_1.fastq.gz
            ├── ERR4079133_2.fastq.gz
            ├── ERR4079285_1.fastq.gz
            └── ERR4079285_2.fastq.gz
```
This repository contains a complete bioinformatics workflow built using a Conda environment on Ubuntu. It includes all major steps of a variant‑analysis and phylogenomics pipeline, with tools primarily installed from the Bioconda channel. The workflow covers quality control, read trimming, homology checking, read alignment, variant calling, and a custom SNP‑to‑Tree pipeline for generating pairwise SNP distances and a 2‑taxon phylogenetic tree.

Tools and stages included in this repository

**FastQC & MultiQC** – quality assessment of raw and trimmed FASTQ files

**Trimmed reads** – preprocessing and adapter removal

**Homology check** – verifying sample identity and contamination

**Read alignment** – mapping reads to a reference genome

**variant calling** - it is used for variant calling 

**SNP‑to‑Tree pipeline** – merging VCFs, generating SNP matrices, computing pairwise distances, and building a 2‑taxon Newick tree

This repository provides scripts, commands, and environment setup instructions to reproduce the entire workflow from FASTQ files to a final phylogenetic tree.

### Full Bioinformatics Pipeline Diagram
```
                 ┌──────────────────────────────┐
                 │        Raw FASTQ Files        │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │   FastQC + MultiQC (QC)       │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │       Trimmed Reads           │
                 │     (Adapter removal)         │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │     Homology / Identity       │
                 │          Checking             │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │      Reference Genome         │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │        BWA index              │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │         BWA mem               │
                 │     (Read Alignment)          │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │     SAM → BAM (samtools)      │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │   Sorted + Indexed BAM        │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │         Variant VCFs          │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │   bcftools merge + index      │
                 │     (merged.vcf.gz)           │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │       vcf2phylip              │
                 │  (FASTA + PHYLIP SNP matrix)  │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │  Python SNP Distance Script   │
                 │  - comparable SNPs            │
                 │  - differences                │
                 │  - percent identity           │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │   Distance Matrix (2×2)       │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │   2‑Taxon Newick Tree         │
                 └──────────────────────────────┘
                               │
                               ▼
                 ┌──────────────────────────────┐
                 │       iTOL Visualization      │
                 └──────────────────────────────┘

```


## the project
I want to create the phylogenetic tree of the **Fusarium oxysporum** it is a filamentous fungus that lives in soil and is known for causing Fusarium wilt, a serious plant disease affecting many crops worldwide.It enters the plant through the **roots**, colonizes the **xylem vessels**, and blocks water transport, causing:
 Wilting,
 Yellowing of leaves,
 Stunted growth and 
 Plant death
 
## conda installation and activation 
I have used **conda** as paskage and environmental manager and among both distribution is used anaconda to install the tools from the channels that is important for my pipeline below are the steps to install the conda 

```
sudo apt update && sudo apt upgrade -y #Update system packages
sudo apt install wget bzip2 -y #Install required tools
wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh #Download the latest Anaconda installer
bash Anaconda3-2024.10-1-Linux-x86_64.sh #Run the installer
source ~/.bashrc  #Activate Anaconda
conda --version #Confirm installation
conda 24.11.3
```
conda environmnet activation
```
conda activate #to see if we have conda
conda list ##to see all the default packages in the conda
conda env list # to see all our environment that we have created
base                  *  /home/ramisha_azhar/anaconda3
the * star means we are basically in this environmnet and base is the default environmnet
```

to create the environmnet
```
conda create --name mynewenv   #the name of my environmnet is mynewenv
conda activate mynewenv   #activate the environment
conda deactivate     #deactivate the environmnet and return to the base
conda env list  #to see the list of all the environmnet created
# conda environments:
#
base                  *  /home/ramisha_azhar/anaconda3
mynewenv                 /home/ramisha_azhar/anaconda3/envs/mynewenv
```

## fastqc files and multi qc files
When we add Bioconda to Conda, we gain access to thousands of bioinformatics packages, but we still need to install each tool
Check if Bioconda channel is configured
```
conda config list

channels:
  - bioconda
  - conda-forge
  - defaults
```
If we  don’t see bioconda, add it:
```
conda config --add channels bioconda

conda config --add channels conda-forge

conda config --add channels defaults
```

This downloads and installs FastQC and Multi QC into our Conda environment.
```
conda install -c bioconda fastqc multiqc
```

we can see if our tools are installed by the below commands

```
fastqc --version   #FastQC v0.12.1

multiqc --version  #multiqc, version 1.32

conda list | egrep 'fastqc|multiqc' 
#fastqc                    0.12.1               hdfd78af_0    bioconda
#multiqc                   1.32               pyhdfd78af_1    bioconda

conda list fastqc
# Name                    Version                   Build  Channel
#fastqc                    0.12.1               hdfd78af_0    bioconda

which fastqc #/home/ramisha_azhar/anaconda3/bin/fastqc #it tells us the exact path

which multiqc #/home/ramisha_azhar/anaconda3/bin/multiqc
```
the output of files for fastqc and multi qc

<img width="262" height="331" alt="image" src="https://github.com/user-attachments/assets/333aaec9-3256-451e-9de3-73d9b7d140af" />


<img width="172" height="97" alt="image" src="https://github.com/user-attachments/assets/da9512eb-04c5-4b03-8790-7dbf161c596c" />

## trimmed reads

Trim Galore is a preprocessing tool for cleaning raw FASTQ files before downstream analysis. It bundles two other tools—Cutadapt and FastQC—so we can trim adapters, remove low‑quality bases, and generate QC reports in a single command.

how to install trim-galore

```
conda install trim-galore

trim_galore --version  #version 0.6.10

 conda list | grep trim-galore
#trim-galore               0.6.10               hdfd78af_1    bioconda

which trim_galore
#/home/ramisha_azhar/anaconda3/bin/trim_galore
```

For Trim Galore , channel order matters because it depends on other tools—mainly Cutadapt and FastQC—and those dependencies must come from the correct channels to avoid version conflict

Make sure our channels are set correctly:
```
conda config --show channels

channels:
  - defaults
  - bioconda
  - conda-forge
```
trim galore command
```
trim_galore --paired -q 20 --fastqc --length 30 -o ERR4079133_1.fastq.gz ERR4079133_2.fastq.gz
```
### output files of the trim.galore

**1. Trimmed FASTQ files**
   
- for Paired‑end reads:
  
sample_R1_val_1.fq.gz    


sample_R2_val_2.fq.gz

**2. FastQC reports (automatically generated)**
   
- For each trimmed FASTQ file, you get:

*_fastqc.html  


*_fastqc.zip

These contain quality plots, adapter content, per‑base quality, etc.

**3. Trimming report**
   
A text summary describing:

how many reads were trimmed

how many bases removed

adapter detection

quality trimming

read length distribution

**- File name example:**
  
sample_R1_val_1_report.txt  


sample_R2_val_2_report.txt

<img width="302" height="382" alt="image" src="https://github.com/user-attachments/assets/8347dc1c-5124-4b35-a9b3-32a1f7648e7a" />

<img width="255" height="340" alt="image" src="https://github.com/user-attachments/assets/3997aaee-dec6-47f8-9e18-c055775fe7ea" />

## Homology check

Mash is a fast genome and metagenome distance estimation tool.
It uses MinHash sketches to quickly compare large sequences.
it gives us a fast, lightweight way to estimate how similar our samples are before we invest time in heavy downstream analysis.


It quickly checks sample similarity using MinHash sketches.

It helps us confirm that your samples belong to the expected species or strain.

It can detect contamination or mislabeled samples early.

It lets you cluster genomes or reads to see overall relationships.

It’s extremely fast, so you get a quality check without running full alignment or assembly.

how to install

```
conda install -c bioconda mash
```

to see if it is installed properly

```
mash --version  #version 2.3

conda list | grep mash
#mash                      2.3                  hc74b729_7    bioconda

which mash   #/home/ramisha_azhar/anaconda3/bin/mash
```
- we have two types of files:

**.msh** → Mash sketch files


**_vs_ref** → Mash comparison results (distance to reference)

**1 .msh files**

These are Mash sketch files — compressed MinHash representations of our FASTQ trimmed reads and our reference genome.

- In your folder:
  
ERR4079133_1_val_1.fq.gz.msh


ERR4079133_2_val_2.fq.gz.msh


ERR4079285_1_val_1.fq.gz.msh


ERR4079285_2_val_2.fq.gz.msh


Fv10027Complete.fasta.msh


Each .msh file is a sketch of one read
The reference genome also has its own sketch (Fv10027Complete.fasta.msh).
These sketches are what Mash uses to compute distances.

**2. _vs_ref files**

These are Mash distance results comparing of each read sketch to our reference sketch.

- In your folder:
  
ERR4079133_1_val_1.fq.gz_vs_ref

 
ERR4079133_2_val_2.fq.gz_vs_ref


ERR4079285_1_val_1.fq.gz_vs_ref


ERR4079285_2_val_2.fq.gz_vs_ref

What they contain:
Each file is a tab‑delimited text output showing:

- sample name

- reference name

- Mash distance

- p‑value

- number of shared hashes

This tells us  how similar each read file is to the reference genome.

we now have:

- Sketches of all our reads

- Sketch of the reference genome

- Distance results showing how close each sample is to the referenc

  <img width="301" height="383" alt="image" src="https://github.com/user-attachments/assets/fdb6e369-8104-4677-8c31-2e14d9ef79cf" />

## Alignment of reads

### Alignment Workflow (BWA + SAMtools)
```  
  FASTQ (trimmed reads)
        │
        ▼
Reference FASTA ──► [BWA index] ──► Indexed reference (.bwt, .pac, .sa, etc.)
        │
        ▼
[BWA mem] ──► SAM file (text alignments)
        │
        ▼
[SAMtools view] ──► BAM file (binary alignments)
        │
        ▼
[SAMtools sort] ──► Sorted BAM file
        │
        ▼
[SAMtools index] ──► BAM index (.bai)
        │
        ▼
Post‑processing ──► Deduplication + Stats (flagstat, markdup)
        │
        ▼
Final outputs: Clean BAM + QC reports
```

  tools we used
  1. bwa
  2. samtools
    
- **Genome indexing (BWA)** → prepares the reference for alignment.
- **BAM indexing (SAM tools)** → prepares alignment files for visualization and downstream analysis
     

Genome mapping using BWA (Burrows-Wheeler Aligner) refers to the process of aligning short DNA sequencing reads to a reference genome using the BWA software tool. 

```
conda install -c bioconda bwa
```

```
conda list bwa
# Name                    Version                   Build      Channel
# bwa                       0.7.18               he4a0461_1    bioconda
```

## **Steps in Genome Mapping with BWA**

**- Prepare Reference Genome**
   
     - Download FASTA file of the reference genome.
    
    ```
    Fv10027Complete.fasta  #this is my reference genome file
    
    to see the file on the terminal
    
    less Fv10027Complete.fasta
    head Fv10027Complete.fasta
    ```
- index the reference genome

```
  bwa index reference.fasta
```

Indexing the genome with BWA is like building a “search engine” for our reference FASTA or fai or fna, so millions of reads can be aligned quickly and efficiently.
these are the reference genome indexing files


.fasta

- .fasta.amb
- .fasta.ann
- .fasta.bwt
- .fasta.pac
- .fasta.sa

Tools like BWA MEM require the reference to be indexed before they can run. Without the index files (.bwt, .pac, .sa, etc.), BWA simply won’t work.

 **- align short-reads**
 
aligning reads to the reference genome using bwa mem command
we have paired end read samples so we use both the reads of the sample  and we will use trimmed samples as part of the pipeline this will generate sam file for us 

 ```
bwa mem reference.fasta read_val_1.fq read_val_2.fq > output.sam
```

Input:

- **Trimmed FASTQ files**
- Indexed reference genome

Output:

- **SAM file** (plain text alignment file)

The SAM file contains:

- Mapping positions
- CIGAR strings
- Flags
- Alignment quality scores
- Metadata

to see sam files

```
head ERR4079133.sam

less ERR4079133.sam
```

 **- Convert SAM → BAM (samtools view)**

 to use samtools commands we need to install samtool first

```
conda install -c bioconda samtools
```

to check if it is installed in our conda environment

```
conda list samtools
# Name                    Version                   Build  Channel
samtools                  1.21                 h50ea8bc_0    bioconda

which samtools
#/home/ramisha_azhar/anaconda3/bin/samtools
```
This gives you:

samtools view

samtools sort

samtools index

samtools flagstat

samtools stats

and many other useful commands

All in one installation.

SAM files are large and slow to process, so they are converted to compressed BAM format.
bam file is same sam file no changes in the coordinate or anything 

```
samtools view -S -b ERR4079133.sam > ERR4079133.bam

#-S is text file
#-b binary file

#to view the bam file
samtools view  ERR4079133.bam | less

```

- Smaller, binary BAM file
- This prints the BAM contents in **SAM-like text format**.
- You can pipe it to `less` or `head` to scroll or preview.

**- Sort the BAM File (samtools sort)**

Sorting arranges reads based on their genomic coordinate.

```
samtools sort alignment.bam -o alignment_sorted.bam
#-o it is output
```

to view the sorted bam file

```
samtools view  ERR4079133.sorted.bam | head
```
We need a sorted BAM file because most downstream tools (indexing, variant calling, visualization) require reads to be ordered by genomic coordinates. Sorting makes access faster, indexing possible, and ensures correct analysis results  

This prints just the header section (the @SQ, @HD, @PG lines):
```
samtools view -H  ERR4079133.bam
samtools view -H  ERR4079133.sorted.bam
```
It shows exactly the same header we saw in our SAM file, but pulled from the BAM.

**- index bam sorted file**
  
```
samtools index sample.sorted.bam
```

it generate statistics of the sorted bam file

```
samtools flagstats  ERR4079133.sorted.bam
```

**- use picard**

```
conda install -c bioconda picard
```

Picard is a toolkit that cleans and improves your BAM files so they are ready for high‑quality downstream analysis, especially variant calling.

we are using this picard tool :

**MarkDuplicates**

This tool:

scans our BAM file which is sorted.bam file

finds PCR duplicates (reads that came from the same original DNA fragment)

removes them (because you set REMOVE_DUPLICATES=true)

writes a clean BAM file

produces a metrics report

our command

```
MarkDuplicates I=sorted.bam O=dedup.bam M=metrics.txt REMOVE_DUPLICATES=true

```
This creates:

sample_dedup.bam → the cleaned BAM

sample_dedup.bam.bai → index

sample_dedup_metrics.txt → duplicate statistics

This dedup BAM is the one you use for variant calling.

Picard cleans your BAM file by removing PCR duplicates so your variant calling is accurate and trustworthy.

the files we get in the end

- Alignment files
  
sample_sorted.bam  

sample_sorted.bam.bai

- alignmnet stat + Dedup files
  
sample_flagstat.txt  

sample_dedup.bam  

sample_dedup.bam.bai  

sample_dedup_metrics.txt

sample_dedup.bam #is the file you use for variant calling.

## variant calling
tools to download 

```
conda install -c bioconda bcftools
bcftools --version
```
there is common error that we see while varient calling that is **(libgsl.so.25)**
it says 
<img width="697" height="41" alt="image" src="https://github.com/user-attachments/assets/6444f879-4b29-4f65-be43-88532f11e155" />

<img width="706" height="60" alt="image" src="https://github.com/user-attachments/assets/714d1870-e978-408b-a427-7accc4a8fa92" />

```
libgsl.so.25 => not found
```
This confirms that bcftools cannot find the GSL library, so the error you saw earlier is because the library is missing or not linked correctly.


to resolve this we should follow certain steps 

- create a new conda environment we can call it variants
  
```
conda create -n variants python=3.10
conda activate variants
conda deactivate             #to return to the base environment

```

- now install all the tools necessary for variant calling in this environment
  
```
 conda install -c bioconda bcftools samtools htslib
 conda install -c conda-forge gsl #in the terminal it get the gsl from the default channel make sure it is from conda forge
```

- conda channel order in our environment
  
```
this is the default setting
 conda config --show channels
 channels:
  - bioconda
  - defaults
  - conda-forge
  - https://repo.anaconda.com/pkgs/main
  - https://repo.anaconda.com/pkgs/r
  ```
we need to download gsl from the conda-forge not the default channel so we have to do some changes for that either we have to change channel order or we use the command that overdo and download the library from the right channel

```
conda config --show channelsconda install -n variant --override-channels -c conda-forge gsl=2.8
```
if the above command does not work use the below one

```
# Install aligned versions

conda install -y -c conda-forge -c bioconda --strict-channel-priority \

bcftools=1.22 htslib=1.22 samtools=1.22 gsl openssl
```
- **bcftools 1.2** and **samtools 1.6** are the main issue. Upgrade both to **1.22** (or 1.21) and keep them aligned with **htslib**.
- Your **GSL 2.8** and **OpenSSL 3.6.0** from **conda‑forge** are good.
- Use **strict channel priority**

to see all the tools installed

```
conda list --show-channel-urls | grep -E 'bcftools|htslib|samtools|openssl|gsl'
```
<img width="691" height="159" alt="image" src="https://github.com/user-attachments/assets/ba72acbb-c9e0-4c5c-9b7d-145ea305f697" />

```
ldd $(which bcftools) | grep gsl

#libgsl.so.25 => /home/ramisha_azhar/anaconda3/envs/variant/bin/../lib/libgsl.so.25 (0x00007f1731215000)
```
now we can see the library is found

these are the step we should take before running the pipeline of variant calling

command for variant calling to generate vcf files of our samples
- bcftools mpileup
- bcftools call
- bcftools view

my input is ERR4079133_dedup.bam which I get at the end of the alignment process 
```
bcftools mpileup -O b -o ERR4079133.bcf -f Fv10027Complete.fasta 
  --threads 8 -q 20 -Q 30 ERR4079133_dedup.bam
```
“Summarize all bases at every position in the genome so I can call SNPs later.”
  
```
bcftools call --ploidy 2 -m -v -o ERR4079133vcf ERR4079133.bcf
```
This step produces your final variant calls.

we can see our vcf files with the help of **bcftools view** command

Shows metadata, sample names, contigs, filters, FORMAT fields

```
bcftools view -h file.vcf.gz
bcftools view  ERR4079133.vcf.gz | less -S # to see whole file while scrolling
```
Count all variant records
```
grep -v -c '^#'  ERR4079133.vcf.gz  #740804
```
Count only SNPs
```
bcftools view -v snps ERR4079133.vcf.gz | grep -v -c '^#' #271009
```
## SNP‑to‑Tree Pipeline

Script should:

Merge the VCFs

Index the merged VCF

Run vcf2phylip.py

Produce:

merged.min2.phy

merged.min4.fasta

```
Multiple VCF.gz files (per-sample variants)
        │
        ▼
[bcftools merge] ──► merged.vcf.gz (multi‑sample VCF)
        │
        ▼
[bcftools index] ──► merged.vcf.gz.csi (indexed VCF)
        │
        ▼
[vcf2phylip] ──► SNP alignment: #Convert merged VCF to FASTA
                 • merged.min2.fasta (FASTA) #Use merged.min4.fasta as our SNP matrix(FASTA SNP matrix)
                 • merged.min2.phy   (PHYLIP)
        │
        ▼
[Python SNP distance script] ──►
                 • SNP difference count
                 • Percent identity
                 • Distance matrix
        │
        ▼
[2‑taxon tree builder] ──► two_taxon_tree.nwk (Newick tree)
        │
        ▼
[iTOL visualization] ──► Final phylogenetic figure
```


these are core bioinformatics tools, and we must have all of them installed for the pipeline to run smoothly 

- bcftools
Used for:

merging VCFs

indexing VCFs

querying SNPs

```
conda install -c bioconda bcftools
```
SNP counts in merged.vcf.gz

```
bcftools view -H merged.vcf.gz | wc -l
#291076
```
- Install vcf2phylip
### Install dependencies (vcf2phylip needs Biopython)

```
conda install -c bioconda biopython
```

### Download the script directly into your conda environment

```
cd $CONDA_PREFIX/bin
wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
chmod +x vcf2phylip.py

#  Test it
vcf2phylip.py -h
```
If we see the help menu, it’s installed correctly.

## Build IQ‑TREE 

It requires at least 3 taxa to infer a phylogenetic tree.
With two samples, IQ‑TREE will always stop 

```
conda install -y -c bioconda iqtree

#Basic IQ‑TREE command
iqtree2 -s merged.min2.fasta -m GTR+G -nt AUTO

```
IQ‑TREE requires **at least 3 sequences or samples** to infer a branching phylogeny.
This is why you get:

Code

**`Alignment must have at least 3 sequences`**

This is not a bug — it’s a mathematical requirement.
we can not bulid the iq-tree because of the lack of the samples

## compute_snp_distance,build_distance_matrix and build_two_taxon_tree

All three are scientifically valid for two samples, and they give you a complete comparative analysis even without IQ‑TREE.

### compute_snp_distance:

icompute_snp_distance() uses bcftools query to extract genotypes from a merged VCF file and compares the two samples SNP‑by‑SNP. It counts the number of comparable SNPs, the number of differences, and calculates the percent identity between the genomes. The function outputs these statistics to a text file. The only external tool required is bcftools, while the rest of the processing is performed in Python.

Query the VCF using bcftools
Used for extracting genotype information from the VCF.
```
bcftools query -f "%CHROM %POS [%GT ]\n" merged.vcf.gz
```
- result
```
Total comparable SNPs: 246932
Different SNPs: 7873
Percent identity: 96.8117%
```
## build_distance_matrix:

build_distance_matrix() takes the number of SNP differences and the total comparable SNPs, calculates the pairwise genetic distance (p‑distance), and formats it into a simple 2×2 distance matrix for the two samples. The matrix is printed and saved to a text file.

That’s the whole idea — it converts your raw SNP counts into a standard distance‑matrix format used in phylogenetics.

```
p = diff / total

                ERR4079133     ERR4079285
ERR4079133        0.0000       0.031883
ERR4079285        0.031883       0.0000
```
Distance from a sample to itself = 0

Distance between the two samples = p

This is the standard format used in phylogenetics.

## build_two_taxon_tree:

build_two_taxon_tree() takes the pairwise genetic distance (p‑distance) between two samples and converts it into a valid Newick‑formatted phylogenetic tree. Because there are only two samples, each branch is assigned half of the total distance, and the function writes the resulting 2‑taxon tree to a .nwk file

```
(ERR4079133:0.015942,ERR4079285:0.015942);
```
## Visualize our tree in iTOL 

https://itol.embl.de/
<img width="1296" height="108" alt="image" src="https://github.com/user-attachments/assets/df622557-b69a-4df2-bf98-00c484a931f1" />
