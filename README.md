# 🦠 SARS-Covid Genomic Analysis Pipeline

This project implements a complete workflow for analyzing blood samples from SARS-CoV-2 patients in **Bash** language. The pipeline processes paired samples to determine the viral genome sequence, assesses the coverage of the S protein gene, identifies genomic variations, determines viral strains via the Pangolin server, and constructs a phylogenetic tree.

---

## 📋 Table of Contents

- [Project Description](#project-description)
- [Technologies & Tools](#technologies--tools)
- [Installation & Setup](#installation--setup)
- [Example Data & Results](#example-data--results)
- [Important Notes](#important-notes)

---

## 📝 Project Description

The project performs the following tasks:
1. **Read Quality Control:**  
   - Assess read quality using 'fasqc'.
   - Identify issues such as high N-content, low quality in short reads, and contamination.

2. **Trimming:**  
   - Remove low-quality bases and adapter sequences using `fastp`.
   - Generate a report with `fasqc` and `multiqc`.

3. **Reference Indexing and Mapping:**  
   - Index the reference genome using Bowtie2.
   - Map reads to the reference and add Read Group (@RG) information.

4. **Post-Mapping Processing:**  
   - Sort and index BAM files using `samtools`.
   - Examine mapping statistics and evaluate close organisms.

5. **S-Gene Extraction and Coverage Analysis:**  
   - Extract coordinates of the S-gene from the annotation.
   - Use `bedtools` to assess the coverage of the S gene.

6. **Variant Calling and Consensus Generation:**  
   - Identify genomic variations using `freebayes`.
   - Merge VCF files and generate consensus sequences with `bcftools`.

7. **Protein Analysis:**  
   - Translate nucleotide sequences into amino acids.
   - Extract the S-gene and combine it with known S-protein sequences.
   - Run BLAST searches and perform multiple sequence alignment.
   - Identify amino acid substitutions and build a phylogenetic tree with `FastTree`.

8. **Contamination Check:**  
   - Use `kraken2` and the Pavian report for contamination analysis.

9. **Detection of Amino Acid Substitutions in the S Protein:**  
   - Extract S-protein sequence and record genomic substitutions.

---

## 🔧 Technologies & Tools

The pipeline is implemented using a variety of command-line tools and bioinformatics software, including:
- **Quality Control & Trimming:** fasqc, fastp, multiqc
- **Mapping & Post-Processing:** Bowtie2, samtools, bedtools
- **Variant Calling & Consensus:** freebayes, bcftools, bgzip, tabix, vcf-merge
- **Protein Analysis:** BLAST (blastp), makeblastdb, linsi (MAFFT), alv, FastTree
- **Contamination Analysis:** kraken2

---

## 🚀 Installation & Setup

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/Elskene-Sashina/NGS_covid.git
   cd your-project-repo
 
2. **Install Dependencies:**
   Ensure you have the required command-line tools installed (e.g., fastp, bowtie2, samtools, etc.).
   Install any necessary bioinformatics libraries and databases (e.g., Kraken2 database, BLAST database).

3. **Set Up File Paths:**
   Update file paths and reference locations in the scripts according to your project directory structure.

---

## 📊 Example Data

Trimming Report:
[View Report](https://drive.google.com/file/d/1e7lG9g4twnREYO8nXEShPGszEyXnsLCR/view)

S-gene Coordinates Extraction:
[View Coordinates](https://drive.google.com/file/d/1Emeike-Y3-FfoEX2W16cWKUAF-8Gat2S/view)

Strain Analysis (Pangolin):
[View Report](https://drive.google.com/file/d/1dpG9EZ1bDP4R569iMsQSoTMkoiND_vWX/view)

Phylogenetic Tree:
[View Tree](https://drive.google.com/file/d/1dpG9EZ1bDP4R569iMsQSoTMkoiND_vWX/view)

---

## 🔍 Important Notes

1. **Parameter Adjustments:**
   Make sure to modify the file paths, sample identifiers, and parameters (e.g., quality thresholds) based on your experimental setup.
 
2. **Tool Versions & Databases:**
   Ensure you have the required command-line tools installed (e.g., fastp, bowtie2, samtools, etc.).
   Install any necessary bioinformatics libraries and databases (e.g., Kraken2 database, BLAST database).

3. **In-line Comments:**
   Each command in the pipeline contains comments for clarification. Please review them to ensure the steps are executed as expected.
