# nf-core/mag benchmarking

This repository contains benchmarking data for the [nf-core/mag](https://nf-co.re/mag/) pipeline. The benchmarks measure how the pipeline performs across datasets and configurations, and give users a reference for estimating the computational resources required for their analyses.

## Datasets

Both datasets contain paired-end short reads (SR) and long reads (LR) from the same sample, so that every assembly strategy is possible: short-read only, long-read only, and hybrid assembly.

### Dataset 1: ZymoBIOMICS Fecal Reference

| Sample     | SR_Gbp | SR_Size (pb) | SR_read_pairs (M) | LR_Gbp | LR_N50 (kb) | LR_reads (M) |
| ---------- | ------ | ------------ | ----------------- | ------ | ----------- | ------------ |
| zymo_fecal | 15.60  | 2x150        | 52.01             | 11.08  | 10.0        | 2.00         |

#### Sources

- Long reads: a subset of the ZymoBIOMICS Fecal Reference published by Oxford Nanopore Technologies (ONT) in 2025 (https://epi2me.nanoporetech.com/zymo_fecal_2025.05/). Flow Cell PAU85136, basecalling SUP.
- Short reads: Illumina dataset uploaded to The ZymoBIOMICS Fecal Reference database (https://fecalreferencedb.com/) by The BioCollective.

### Dataset 2: Maghini et al. (2025) human gut metagenomes

| Sample | SR_Gbp | SR_Size (pb) | SR_read_pairs (M) | LR_Gbp | LR_N50 (kb) | LR_reads (M) |
| ------ | ------ | ------------ | ----------------- | ------ | ----------- | ------------ |
| D01    | 9.44   | 2x150        | 31.91             | 6.98   | 6.3         | 3.43         |
| D02    | 6.89   | 2x150        | 23.26             | 10.50  | 11.8        | 3.73         |
| D03    | 15.81  | 2x150        | 53.59             | 3.82   | 3.9         | 2.14         |
| D04    | 18.96  | 2x150        | 63.98             | 9.68   | 12.1        | 2.78         |
| D05    | 12.01  | 2x150        | 40.53             | 13.31  | 7.4         | 4.11         |
| D06    | 14.26  | 2x150        | 48.28             | 13.67  | 10.4        | 3.98         |
| D07    | 31.69  | 2x150        | 107.27            | 6.71   | 3.4         | 4.29         |
| D08    | 10.15  | 2x150        | 34.32             | 12.31  | 13.3        | 2.75         |
| D09    | 9.26   | 2x150        | 31.30             | 9.24   | 6.6         | 5.07         |
| D10    | 9.76   | 2x150        | 32.80             | 4.12   | 5.2         | 2.36         |

#### Sources

Data published by Maghini et al. (2025) in "Illumina complete long read assay yields contiguous bacterial genomes from human gut metagenomes" (https://doi.org/10.1128/msystems.01531-24). The dataset is available at the NCBI SRA under BioProject [PRJNA940499](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA940499).

## Experimental setup

### Software and hardware

- Software: nf-core/mag v5.5.0, Nextflow v26.04.6, Apptainer v1.5.0
- Hardware: HPC cluster based on AMD EPYC Zen 2 processors, each node with 48 cores and 512 GB RAM. Slurm as workload manager.

### Pipeline parameters

Each dataset was processed in a single run, with minimal changes over pipeline defaults. All assemblers (2 short-read, 2 long-read, 1 hybrid) and all binners (6) were enabled, plus polishing of the long-read assemblies, assembly QC (QUAST, ALE), bin refinement with DAS Tool, bin QC by BUSCO, CheckM2, QUAST and GUNC, annotation with Prodigal and Prokka, and taxonomic classification by GTDB-Tk and CATpack (ZymoBIOMICS only).

Because `postbinning_input` is set to `both`, every post-binning step runs on the bins of each of the 6 binners and on the DAS Tool refined set: 7 bin sets per assembly, so 35 per sample across the 5 assemblies.

Both runs share the parameters below, with one exception: `cat_db` was set only for the ZymoBIOMICS dataset, because of CATpack's heavy computational requirements.

```yaml
input: input/samplesheet.csv
outdir: results

run_pypolca: true

refine_bins_dastool: true
postbinning_input: both

exclude_unbins_from_postbinning: true

busco_db_lineage: auto_prok
busco_clean: true
run_checkm2: true
run_gunc: true

busco_db: input/databases/busco/20260522/
checkm2_db: input/databases/checkm2/v3/uniref100.KO.1.dmnd
gtdb_db: input/databases/gtdb/release232/
gunc_db: input/databases/gunc/gunc_db_progenomes2.1.dmnd
cat_db: input/databases/catpack/20231120_CAT_gtdb/ # ZymoBIOMICS only
```

### Resource overrides

Assemblers and binners all get 8 CPUs on a single socket, so that their resource use stays comparable. metaSPAdes gets more memory than the pipeline default to avoid out-of-memory failures, and COMEBin a longer time limit. Prokka errors are ignored, so that a failure on an individual bin does not stop the run.

```groovy
process {
    withName: 'PROKKA' {
        errorStrategy = 'ignore'
    }
    withName: 'METASPADES|METASPADESHYBRID' {
        memory = 125.GB
    }

    withName: 'MEGAHIT|METASPADES|METASPADESHYBRID|FLYE|METAMDBG_ASM' {
        cpus           = 8
        clusterOptions = '--cores-per-socket=8'
    }
    withName: 'METABAT2_METABAT2|MAXBIN2|CONCOCT_CONCOCT|SEMIBIN_SINGLEEASYBIN|METABINNER_METABINNER|COMEBIN_RUNCOMEBIN' {
        cpus           = 8
        clusterOptions = '--cores-per-socket=8'
    }

    withName: 'COMEBIN_RUNCOMEBIN' {
        time = 48.h
    }
}
```

## Analysis

The notebooks read the pipeline outputs under `data/` and share the code in `src/mag_bench/`.

### [notebooks/01_basic_metrics.ipynb](notebooks/01_basic_metrics.ipynb)

Assembly and bin metrics per dataset and assembler, from the QUAST summaries and the bin summary tables: mean assembly length, contig count and N50, then the number of bins recovered and their mean length. Describes what the pipeline produced, not what it cost.

### [notebooks/02_benchmarking_metrics.ipynb](notebooks/02_benchmarking_metrics.ipynb)

The six supplementary figures of the manuscript, written to `figures/`. Compute use and storage per analysis stage, resource footprint per assembler and binner, the most resource-intensive processes, and how assemblers and binners scale with sequencing input and assembly size. Built from the extended execution traces.

### [notebooks/03_configuration_estimates.ipynb](notebooks/03_configuration_estimates.ipynb)

What two minimal configurations would have cost, one short-read and one long-read, obtained by selecting from the same traces only the tasks each configuration would have launched. Reports jobs, CPU core-hours, work-directory size and peak memory per sample, both in total and broken down by stage. Runtime is left out on purpose: wall-clock depends on the critical path and on cluster contention, not on the sum of task times.

### [notebooks/04_process_resources.ipynb](notebooks/04_process_resources.ipynb)

A lookup table of what every process costs per task: CPU time, peak memory, runtime and work-directory size, as a median and an observed range. The numbers are what tasks consumed, not what the pipeline reserved for them, and the ranges pool both datasets and every assembly a process ran on. The same table lives in [data/process_resources.csv](data/process_resources.csv), for looking things up without running the notebook.
