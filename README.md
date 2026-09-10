# Bacterial Genome Analysis Tool

A Python-based bioinformatics tool for analyzing and comparing bacterial genomes using **Biopython**.

The program loads bacterial genome sequences from FASTA files and provides genome statistics, comparative analysis, k-mer analysis, visualizations, and basic biological interpretation.

## Features

### Genome Loading

* Automatically loads bacterial genomes from FASTA files
* Uses Biopython for sequence parsing
* Handles genome loading errors

### Genome Selection

* Displays all available genomes
* Allows the user to select multiple genomes
* Validates genome selections and prevents duplicate selections

### Individual Genome Analysis

For each selected genome, the tool calculates:

* Genome length
* A, C, G, and T base counts
* GC content
* AT content

### Genome Comparison

Compares multiple bacterial genomes and identifies:

* Highest GC content
* Lowest GC content
* GC content difference
* Largest genome
* Smallest genome
* Genome size difference

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

## Technologies Used

* **Python**
* **Biopython**
* **Matplotlib**
* **Pathlib**

## Project Structure

```text
Bacterial-Genome-Analysis-Tool/
│
├── data/
│   ├── e_coli.fasta
│   ├── b_subtilis.fasta
│   └── s_aureus.fasta
│
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

These genomes provide differences in genome size, GC content, and sequence composition that can be explored using the tool.

## Installation

Clone the repository and navigate into the project directory:

```bash
git clone https://github.com/simplyhamdan/Bacteria-Genome-Analysis-Tool.git
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

The program will display the available genomes and allow you to select which genomes to analyze.

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
* Differences in GC content between the organisms
* Differences in k-mer frequency patterns

These comparisons demonstrate how bacterial genomes can vary in size and nucleotide composition.

## Purpose

This project was developed as a practical introduction to **computational genomics and bioinformatics programming**.

It demonstrates the use of Python to process biological sequence data, perform comparative genomic analysis, and generate visual representations of genomic characteristics.

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

## Author

**Hamdan Sajith**

Biotechnology student interested in **bioinformatics, computational biology, and genomic data analysis**.
