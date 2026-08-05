# nf-core/mag benchmarking

This repository contains benchmarking data for the [nf-core/mag](https://nf-co.re/mag/) pipeline. The data is used to evaluate the performance of the pipeline on different datasets and configurations and provide a reference for users to estimate the computational resources required for their own analyses.

## Datasets

Both datasets contain paired-end short reads (SR) and long reads (LR) from the same sample, to enable all possible assembly strategies (short-read only, long-read only, hybrid assembly). The datasets are summarized in the tables below.

### Dataset 1: ZymoBIOMICS Fecal Reference

| Sample     | SR_Gbp | SR_Size (pb) | LR_Gbp | LR_N50 (kb) | LR_reads (M) |
| ---------- | ------ | ------------ | ------ | ----------- | ------------ |
| zymo_fecal | 15.60  | 2x150        | 11.08  | 10.0        | 2.00         |

#### Sources

- Long reads: a subset of the ZymoBIOMICS Fecal Reference published by Oxford Nanopore Technologies (ONT) in 2025 (https://epi2me.nanoporetech.com/zymo_fecal_2025.05/). Flow Cell PAU85136, basecalling SUP.
- Short reads: Illumina dataset uploaded to The ZymoBIOMICS Fecal Reference database (https://fecalreferencedb.com/) by The BioCollective.

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

### Software and hardware

- Software: nf-core/mag v5.5.0, Nextflow v26.04.6, Apptainer v1.5.0
- Hardware: HPC cluster based on AMD EPYC Zen 2 processors, each node with 48 cores and 512 GB RAM. Slurm as workload manager.

### Pipeline parameters

Each dataset was processed in a single run, with minimal changes over pipeline defaults. All assemblers (2 short-read, 2 long-read, 1 hybrid) and all binners (6) were enabled, plus polishing of the long-read assemblies, assembly QC (QUAST, ALE), bin refinement with DAS Tool, bin QC by BUSCO, CheckM2, QUAST and GUNC, annotation with Prodigal and Prokka, and taxonomic classification by GTDB-Tk and CATpack (ZymoBIOMICS only).

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

Assemblers and binners all get 8 CPUs on a single socket, so that their resource use is comparable to each other. metaSPAdes gets more memory than the pipeline default to avoid out-of-memory failures, and COMEBin a longer time limit. Prokka errors are ignored, so that a failure on an individual bin does not stop the run.

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
