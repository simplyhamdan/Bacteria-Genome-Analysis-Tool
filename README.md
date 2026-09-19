# Bacterial Genome Analysis Tool

A Python-based bioinformatics tool for analysing and comparing bacterial genomes using **Biopython**.

The program loads bacterial genome sequences from FASTA files and performs genome statistics, comparative analysis, k-mer analysis, visualisation, and basic biological interpretation.

## Features

### Genome Loading

* Automatically loads bacterial genomes from FASTA files
* Uses Biopython for sequence parsing
* Handles genome loading errors

### Genome Selection

* Displays available genomes
* Allows multiple genomes to be selected for analysis
* Validates genome selections and prevents duplicate selections

### Individual Genome Analysis

For each selected genome, the tool calculates:

* Genome length
* A, C, G, and T base counts
* GC content
* AT content

### Genome Comparison

Compares selected bacterial genomes and identifies:

* Highest and lowest GC content
* GC content differences
* Largest and smallest genomes
* Genome size differences

### K-mer Analysis

* Counts k-mers within genome sequences
* Identifies the most frequent k-mers
* Compares k-mer frequency patterns between genomes

### Visualizations

Generates plots for:

* Genome size comparison
* GC content comparison
* AT content comparison
* Base composition
* K-mer frequency

### Biological Interpretation

Provides basic interpretations of:

* Low, moderate, and high GC content
* Differences in genome size
* Possible biological factors contributing to genome-size variation

## Project Structure

```text
Bacterial-Genome-Analysis-Tool/
├── data/
│   ├── e_coli.fasta
│   ├── b_subtilis.fasta
│   └── s_aureus.fasta
├── main.py
├── README.md
└── .gitignore
```

## Dataset

The project uses complete bacterial genome sequences in FASTA format obtained from the **NCBI RefSeq** database.

The current dataset contains:

* *Escherichia coli* K-12 MG1655
* *Bacillus subtilis* subsp. subtilis 168
* *Staphylococcus aureus* subsp. aureus NCTC 8325

These genomes provide different genome sizes, GC contents, and sequence composition for comparative analysis.

## Technologies

* Python
* Biopython
* Matplotlib
* Pathlib

## Installation

Clone the repository and navigate into the project directory:

```bash
git clone <repository-url>
cd Bacterial-Genome-Analysis-Tool
```

Install the required Python packages:

```bash
pip install biopython matplotlib
```

## How to Run

Run the program with:

```bash
python main.py
```

The program displays the available genomes and allows you to select which genomes to analyse.

Example:

```text
Available genomes:
[1] NC_000964.3
[2] NC_000913.3
[3] NC_007795.1

Select genomes to analyze (e.g. 1,3): 1,2,3
```

You can then choose from the analysis menu:

```text
Analysis Menu:
[1] Individual Genome Analysis
[2] Compare Genomes
[3] K-mer Analysis
[4] Visualizations
[5] Biological Interpretation
[6] Exit
```

## Example Analysis

For the three genomes in the dataset, the tool can identify differences such as:

* *E. coli* having the largest genome
* *S. aureus* having the smallest genome
* Differences in GC content between organisms
* Differences in k-mer frequency patterns

These comparisons demonstrate how bacterial genomes can vary in size and nucleotide composition.

## What I Learned

This project built on the sequence-analysis concepts from the previous projects and introduced larger-scale biological data analysis.

Key concepts included:

* Working with complete bacterial genome sequences
* FASTA file parsing with Biopython
* Genome statistics and base composition
* Comparative genomic analysis
* K-mer analysis
* Biological data visualization
* Interpreting basic genomic characteristics
* Organizing biological data for computational analysis

## Purpose

This project was developed as a practical introduction to **computational genomics and bioinformatics programming**.

It demonstrates how Python can be used to process biological sequence data, perform comparative genomic analysis, and visualize genomic characteristics.

## Future Improvements

Possible future extensions include:

* Gene prediction and annotation
* ORF analysis across genomes
* Codon usage analysis
* GC content analysis across genomic regions
* Sequence alignment
* SNP and mutation detection
* Genome similarity analysis
* Automated report generation

## Project Series

Part of a series of progressively more advanced bioinformatics projects.

**Project 1:** DNA Sequencing Toolkit — Fundamental DNA sequence operations
**Project 2:** Protein Sequence Analyzer — Protein sequence and physicochemical analysis
**Project 3:** DNA Sequence Analyzer — Comparative DNA sequence analysis
**Project 4:** Bacterial Genome Analysis Tool — Genome-scale sequence analysis and visualization

The projects progress from individual sequence processing to comparative analysis and genome-scale biological data analysis.
