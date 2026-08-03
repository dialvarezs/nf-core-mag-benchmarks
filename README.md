# nf-core/mag benchmarking

This repository contains benchmarking data for the nf-core/mag pipeline. The data is used to evaluate the performance of the pipeline on different datasets and configurations and provide a reference for users to estimate the computational resources required for their own analyses.

## Datasets

Both datasets used for benchmarking contains paired-end short reads (SR) and long reads (LR) from the same sample, to enable all possible assembly strategies (short-read only, long-read only, hybrid assembly). The datasets are summarized in the tables below.

### Dataset 1: ZymoBIOMICS Fecal Reference

| Sample     | SR_Gbp | SR_Size (pb) | LR_Gbp | LR_N50 (kb) | LR_reads (M) |
| ---------- | ------ | ------------ | ------ | ----------- | ------------ |
| zymo_fecal | 15.60  | 2x150        | 11.08  | 10.0        | 2.00         |

#### Sources

- Long reads: a subset of the ZymoBIOMICS Fecal Reference published by Oxford Nanopore Technologies (ONT) in 2025 (https://epi2me.nanoporetech.com/zymo_fecal_2025.05/). Flow Cell PAU85136, basecalling SUP.
- Short reads: Illumina dataset uploaded to The ZymoBIOMICS Fecal Reference database (https://fecalreferencedb.com/) by the The BioCollective.

### Dataset 2: Maghini et al. (2025) human gut metagenomes

| Sample | SR_Gbp | SR_Size (pb) | LR_Gbp | LR_N50 (kb) | LR_reads (M) |
| ------ | ------ | ------------ | ------ | ----------- | ------------ |
| D01    | 9.44   | 2x150        | 6.98   | 6.3         | 3.43         |
| D02    | 6.89   | 2x150        | 10.50  | 11.8        | 3.73         |
| D03    | 15.81  | 2x150        | 3.82   | 3.9         | 2.14         |
| D04    | 18.96  | 2x150        | 9.68   | 12.1        | 2.78         |
| D05    | 12.01  | 2x150        | 13.31  | 7.4         | 4.11         |
| D06    | 14.26  | 2x150        | 13.67  | 10.4        | 3.98         |
| D07    | 31.69  | 2x150        | 6.71   | 3.4         | 4.29         |
| D08    | 10.15  | 2x150        | 12.31  | 13.3        | 2.75         |
| D09    | 9.26   | 2x150        | 9.24   | 6.6         | 5.07         |
| D10    | 9.76   | 2x150        | 4.12   | 5.2         | 2.36         |

#### Sources

Data published by Maghini et al. (2025) in "Illumina complete long read assay yields contiguous bacterial genomes from human gut metagenomes" (https://doi.org/10.1128/msystems.01531-24). The dataset is available at the NCBI SRA under BioProject [PRJNA940499](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA940499).

## Experimental setup

### Setup

- Software: nf-core/mag v5.5.0
- Hardware: a HPC cluster based on AMD EPYC Zen 2 processors, each node with 48 cores and 512 GB RAM. The cluster uses a Slurm workload manager.
