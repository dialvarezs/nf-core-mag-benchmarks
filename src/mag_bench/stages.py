"""Group Nextflow processes into the pipeline stages used in the figures.

A trace row is one *task*, and the number of tasks per stage varies by
orders of magnitude: an assembler runs once per sample, while Prokka runs
once per bin (>10 000 tasks here). Aggregating by stage is what makes the
runs comparable, so :data:`PROCESSES` is the backbone of every figure.

Unlike the CPU allocations, which are recorded in the run outputs, this
mapping is editorial: nothing in the pipeline config knows that QUAST-on-
assemblies and QUAST-on-bins belong at opposite ends of the workflow, or
that seqkit is a size filter rather than a QC step. It is maintained by
hand on purpose, and an unmapped process raises rather than falling into
a catch-all, so a pipeline bump surfaces as a review moment.
"""

from __future__ import annotations

import polars as pl

# Display order, coarse-to-fine along the pipeline.
STAGE_ORDER = [
    "Short-read preprocessing",
    "Long-read preprocessing",
    "Short-read assembly",
    "Hybrid assembly",
    "Long-read assembly",
    "Assembly polishing (pypolca)",
    "Assembly QC (QUAST)",
    "Assembly QC (ALE)",
    "Read mapping + depth",
    "Binning",
    "Bin refinement (DAS Tool)",
    "Bin QC (BUSCO)",
    "Bin QC (CheckM2)",
    "Bin QC (QUAST)",
    "Bin QC (GUNC)",
    "Annotation (Prodigal/Prokka)",
    "Taxonomic classification (GTDB-Tk)",
    "Reporting",
]

#: Granularity = what one task of a process corresponds to. It drives the
#: point labels ("points = individual assemblies") and explains why stage
#: totals are dominated by processes that run once per bin.
GRANULARITY_ORDER = [
    "per sample",
    "per assembly",
    "per assembly x binner",
    "per bin",
    "once per run",
]

#: process -> (stage, granularity). Verified against the `include` statements
#: and call sites of nf-core/mag; the non-obvious calls are commented.
PROCESSES: dict[str, tuple[str, str]] = {
    "FASTQC_RAW": ("Short-read preprocessing", "per sample"),
    "FASTP": ("Short-read preprocessing", "per sample"),
    "FASTQC_TRIMMED": ("Short-read preprocessing", "per sample"),
    "BOWTIE2_PHIX_REMOVAL_BUILD": ("Short-read preprocessing", "once per run"),
    "BOWTIE2_PHIX_REMOVAL_ALIGN": ("Short-read preprocessing", "per sample"),
    "NANOPLOT_RAW": ("Long-read preprocessing", "per sample"),
    "PORECHOP_ABI": ("Long-read preprocessing", "per sample"),
    "CHOPPER": ("Long-read preprocessing", "per sample"),
    "NANOPLOT_FILTERED": ("Long-read preprocessing", "per sample"),
    "MEGAHIT": ("Short-read assembly", "per sample"),
    "METASPADES": ("Short-read assembly", "per sample"),
    "GUNZIP_SHORTREAD_ASSEMBLIES": ("Short-read assembly", "per sample"),
    "METASPADESHYBRID": ("Hybrid assembly", "per sample"),
    # CAT_FASTQ alias: pools the long reads that feed metaSPAdes-hybrid
    "POOL_LONG_READS": ("Hybrid assembly", "per sample"),
    "FLYE": ("Long-read assembly", "per sample"),
    "METAMDBG_ASM": ("Long-read assembly", "per sample"),
    "GUNZIP_LONGREAD_ASSEMBLIES": ("Long-read assembly", "per sample"),
    "PYPOLCA_RUN": ("Assembly polishing (pypolca)", "per sample"),
    # QUAST is called twice in the pipeline, on assemblies and on bins; the
    # two calls sit at opposite ends of the workflow and differ by ~7x in
    # task count, so they are kept apart.
    "QUAST": ("Assembly QC (QUAST)", "per assembly"),
    "ALE": ("Assembly QC (ALE)", "per assembly"),
    "BOWTIE2_ASSEMBLY_BUILD": ("Read mapping + depth", "per assembly"),
    "BOWTIE2_ASSEMBLY_ALIGN": ("Read mapping + depth", "per assembly"),
    "MINIMAP2_ASSEMBLY_INDEX": ("Read mapping + depth", "per assembly"),
    "MINIMAP2_ASSEMBLY_ALIGN": ("Read mapping + depth", "per assembly"),
    "METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_SHORTREAD": ("Read mapping + depth", "per assembly"),
    "METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_LONGREAD": ("Read mapping + depth", "per assembly"),
    "CONVERT_DEPTHS": ("Read mapping + depth", "per assembly"),
    "MAG_DEPTHS": ("Read mapping + depth", "per assembly x binner"),
    "MAG_DEPTHS_SUMMARY": ("Reporting", "once per run"),
    "METABAT2_METABAT2": ("Binning", "per assembly x binner"),
    "MAXBIN2": ("Binning", "per assembly x binner"),
    "ADJUST_MAXBIN2_EXT": ("Binning", "per assembly x binner"),
    "CONCOCT_CUTUPFASTA": ("Binning", "per assembly x binner"),
    "CONCOCT_CONCOCTCOVERAGETABLE": ("Binning", "per assembly x binner"),
    "CONCOCT_CONCOCT": ("Binning", "per assembly x binner"),
    "CONCOCT_MERGECUTUPCLUSTERING": ("Binning", "per assembly x binner"),
    "CONCOCT_EXTRACTFASTABINS": ("Binning", "per assembly x binner"),
    "SEMIBIN_SINGLEEASYBIN": ("Binning", "per assembly x binner"),
    "METABINNER_TOOSHORT": ("Binning", "per assembly x binner"),
    "METABINNER_KMER": ("Binning", "per assembly x binner"),
    "METABINNER_METABINNER": ("Binning", "per assembly x binner"),
    "METABINNER_BINS": ("Binning", "per assembly x binner"),
    "COMEBIN_RUNCOMEBIN": ("Binning", "per assembly x binner"),
    "SPLIT_FASTA": ("Binning", "per assembly x binner"),
    # not a QC step: seqkit only measures each bin so the --bin_min_size /
    # --bin_max_size filter can be applied
    "SEQKIT_STATS": ("Binning", "per assembly x binner"),
    "RENAME_PREDASTOOL": ("Bin refinement (DAS Tool)", "per assembly x binner"),
    "DASTOOL_FASTATOCONTIG2BIN": ("Bin refinement (DAS Tool)", "per bin"),
    "DASTOOL_DASTOOL": ("Bin refinement (DAS Tool)", "per assembly"),
    "RENAME_POSTDASTOOL": ("Bin refinement (DAS Tool)", "per assembly"),
    "BUSCO_BUSCO": ("Bin QC (BUSCO)", "per assembly x binner"),
    "CONCAT_BUSCO_TSV": ("Bin QC (BUSCO)", "once per run"),
    "CHECKM2_PREDICT": ("Bin QC (CheckM2)", "per assembly x binner"),
    "CONCAT_CHECKM2_TSV": ("Bin QC (CheckM2)", "once per run"),
    "QUAST_BINS": ("Bin QC (QUAST)", "per assembly x binner"),
    "CONCAT_QUAST_SUMMARY": ("Bin QC (QUAST)", "once per run"),
    "GUNC_RUN": ("Bin QC (GUNC)", "per assembly x binner"),
    "CONCAT_GUNC_TSV": ("Bin QC (GUNC)", "once per run"),
    "PRODIGAL": ("Annotation (Prodigal/Prokka)", "per assembly"),
    "PROKKA": ("Annotation (Prodigal/Prokka)", "per bin"),
    "GTDBTK_CLASSIFYWF": ("Taxonomic classification (GTDB-Tk)", "per assembly x binner"),
    "GTDBTK_SUMMARY": ("Taxonomic classification (GTDB-Tk)", "once per run"),
    "MULTIQC": ("Reporting", "once per run"),
    "BIN_SUMMARY": ("Reporting", "once per run"),
}

STAGE_OF = {process: stage for process, (stage, _) in PROCESSES.items()}
GRANULARITY_OF = {process: grain for process, (_, grain) in PROCESSES.items()}

#: process -> (assembler token as it appears in Nextflow task tags, display
#: name). One entry per tool, so the roster, the tag vocabulary and the
#: figure labels cannot drift apart.
ASSEMBLER_TOOLS = {
    "MEGAHIT": ("MEGAHIT", "MEGAHIT\n(short-read)"),
    "METASPADES": ("SPAdes", "metaSPAdes\n(short-read)"),
    "METASPADESHYBRID": ("SPAdesHybrid", "metaSPAdes\n(hybrid)"),
    "FLYE": ("FLYE", "metaFlye\n(long-read)"),
    "METAMDBG_ASM": ("METAMDBG", "metaMDBG\n(long-read)"),
}
BINNER_TOOLS = {
    "METABAT2_METABAT2": ("MetaBAT2", "MetaBAT2"),
    "MAXBIN2": ("MaxBin2", "MaxBin2"),
    "CONCOCT_CONCOCT": ("CONCOCT", "CONCOCT"),
    "SEMIBIN_SINGLEEASYBIN": ("SemiBin2", "SemiBin2"),
    "METABINNER_METABINNER": ("MetaBinner", "MetaBinner"),
    "COMEBIN_RUNCOMEBIN": ("COMEBin", "COMEBin"),
}

TOOL_NAMES = {p: name for p, (_, name) in (ASSEMBLER_TOOLS | BINNER_TOOLS).items()}
ASSEMBLERS = [name for _, name in ASSEMBLER_TOOLS.values()]
BINNERS = [name for _, name in BINNER_TOOLS.values()]
#: process -> tag token, used to line assembler tasks up with QUAST stats.
TAG_OF_ASSEMBLER = {p: tag for p, (tag, _) in ASSEMBLER_TOOLS.items()}


def annotate(trace: pl.DataFrame) -> pl.DataFrame:
    """Add stage / granularity / display-name columns to a parsed trace.

    An unmapped process raises here rather than landing in a catch-all
    bucket, so a new pipeline release fails on load instead of quietly
    reshaping a figure.
    """
    unmapped = set(trace["process"].unique()) - set(PROCESSES)
    if unmapped:
        raise KeyError(f"processes missing from stages.PROCESSES: {sorted(unmapped)}")
    return trace.with_columns(
        stage=pl.col("process").replace_strict(STAGE_OF),
        granularity=pl.col("process").replace_strict(GRANULARITY_OF),
        tool=pl.col("process").replace_strict(TOOL_NAMES, default=None, return_dtype=pl.String),
    )


def parse_tag(trace: pl.DataFrame) -> pl.DataFrame:
    """Split the Nextflow task tag into assembler and sample.

    Tags are positional and vary by process, e.g. ``D01`` (assembler),
    ``FLYE-D01`` (binner), ``FLYE-MetaBAT2-unclassified-unrefined-D01``
    (bin QC). Matching known tokens is more robust than splitting on "-",
    because sample names themselves contain hyphens in other datasets.
    Longest token first, so ``SPAdesHybrid`` is not shadowed by ``SPAdes``.
    """
    tokens = sorted(set(TAG_OF_ASSEMBLER.values()), key=len, reverse=True)
    return trace.with_columns(
        tag_assembler=pl.col("tag").str.extract(rf"^({'|'.join(tokens)})-"),
        # The sample is always the last hyphen-separated token, but it picks
        # up decorations: ".<bin number>" on per-bin tasks and "_run<n>_raw"
        # on the short-read preprocessing tasks (one row per sequencing run).
        sample=(
            pl.col("tag")
            .str.split("-")
            .list.last()
            .str.replace(r"\..*$", "")
            .str.replace(r"_run\d+.*$", "")
        ),
    )
