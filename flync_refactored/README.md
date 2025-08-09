# FLYNC Refactored

This repository contains a refactored version of the FLYNC pipeline for lncRNA discovery and classification in *Drosophila melanogaster*. The original pipeline was a collection of bash, R, and Python scripts, orchestrated by a master bash script. This version refactors the pipeline into a modern, reproducible, and scalable workflow using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

## Overview

The pipeline performs the following steps:
1.  Downloads RNA-seq data from SRA.
2.  Maps the reads to the reference genome using HISAT2.
3.  Assembles transcripts using StringTie.
4.  Merges the assemblies and compares them to the reference annotation.
5.  Predicts the coding potential of new transcripts using CPAT.
6.  Classifies new transcripts as coding or non-coding.
7.  Performs differential expression analysis using Ballgown (if metadata is provided).
8.  Extracts features for the machine learning model.
9.  Predicts lncRNAs using a pre-trained machine learning model.
10. Generates a final table with the results.

## Installation

1.  **Clone this repository:**
    ```bash
    git clone <repository-url>
    cd flync_refactored
    ```

2.  **Create the conda environment:**
    All the required dependencies are listed in the `environment.yaml` file. You can create the conda environment using the following command:
    ```bash
    conda env create -f environment.yaml
    ```

3.  **Activate the conda environment:**
    ```bash
    conda activate flync_refactored
    ```

## Usage

To run the pipeline, you need to provide a `samples.txt` file with the SRA accession numbers, one per line. A `metadata.csv` file is also required for the differential expression analysis.

Once you have the input files, you can run the pipeline using the following command:
```bash
snakemake --cores <number-of-cores>
```

This will run the entire pipeline and generate the final results in the `results/` directory.

## Workflow

The pipeline is defined in the `Snakefile`. It consists of a series of rules, where each rule corresponds to a step in the pipeline. Snakemake automatically determines the order of execution of the rules based on their dependencies.

## Output

The main output files are:
-   `results/final-results.csv`: A table with the final results, including the machine learning classification, differential expression analysis, and transcript expression levels.
-   `results/new-non-coding.gtf`: A GTF file with the new non-coding transcripts.
-   `results/new-coding.gtf`: A GTF file with the new coding transcripts.
