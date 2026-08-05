"""Map Nextflow processes to analysis stages and tools."""

from __future__ import annotations

import polars as pl

# Pipeline display order.
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
    "Taxonomic classification (CAT)",
    "Taxonomic classification (GTDB-Tk)",
    "Reporting",
]

#: Units represented by individual tasks.
GRANULARITY_ORDER = [
    "per sample",
    "per assembly",
    "per assembly x binner",
    "per bin",
    "once per run",
]

#: Process-to-stage and task-granularity mapping.
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
    # CAT_FASTQ alias for reads passed to metaSPAdes-hybrid.
    "POOL_LONG_READS": ("Hybrid assembly", "per sample"),
    "FLYE": ("Long-read assembly", "per sample"),
    "METAMDBG_ASM": ("Long-read assembly", "per sample"),
    "GUNZIP_LONGREAD_ASSEMBLIES": ("Long-read assembly", "per sample"),
    "PYPOLCA_RUN": ("Assembly polishing (pypolca)", "per sample"),
    # Assembly-level QUAST; bin-level QUAST is QUAST_BINS.
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
    "METABAT2_METABAT2": ("Binning", "per assembly"),
    "MAXBIN2": ("Binning", "per assembly"),
    "ADJUST_MAXBIN2_EXT": ("Binning", "per assembly"),
    "CONCOCT_CUTUPFASTA": ("Binning", "per assembly"),
    "CONCOCT_CONCOCTCOVERAGETABLE": ("Binning", "per assembly"),
    "CONCOCT_CONCOCT": ("Binning", "per assembly"),
    "CONCOCT_MERGECUTUPCLUSTERING": ("Binning", "per assembly"),
    "CONCOCT_EXTRACTFASTABINS": ("Binning", "per assembly"),
    "SEMIBIN_SINGLEEASYBIN": ("Binning", "per assembly"),
    "METABINNER_TOOSHORT": ("Binning", "per assembly"),
    "METABINNER_KMER": ("Binning", "per assembly"),
    "METABINNER_METABINNER": ("Binning", "per assembly"),
    "METABINNER_BINS": ("Binning", "per assembly"),
    "COMEBIN_RUNCOMEBIN": ("Binning", "per assembly"),
    "SPLIT_FASTA": ("Binning", "per assembly x binner"),
    # Used for bin-size filtering, not QC.
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
    "CATPACK_BINS": ("Taxonomic classification (CAT)", "per assembly x binner"),
    "CATPACK_ADDNAMES_BINS": ("Taxonomic classification (CAT)", "per assembly x binner"),
    "CATPACK_SUMMARISE_BINS": ("Taxonomic classification (CAT)", "per assembly x binner"),
    "GTDBTK_CLASSIFYWF": ("Taxonomic classification (GTDB-Tk)", "per assembly x binner"),
    "GTDBTK_SUMMARY": ("Taxonomic classification (GTDB-Tk)", "once per run"),
    "MULTIQC": ("Reporting", "once per run"),
    "BIN_SUMMARY": ("Reporting", "once per run"),
}

STAGE_OF = {process: stage for process, (stage, _) in PROCESSES.items()}
GRANULARITY_OF = {process: grain for process, (_, grain) in PROCESSES.items()}

#: Process-to-tag-token and display-name mapping.
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
#: Process-to-tag token for matching QUAST statistics.
TAG_OF_ASSEMBLER = {p: tag for p, (tag, _) in ASSEMBLER_TOOLS.items()}


def annotate(trace: pl.DataFrame) -> pl.DataFrame:
    """Add stage, granularity, and tool labels to a trace."""
    unmapped = set(trace["process"].unique()) - set(PROCESSES)
    if unmapped:
        raise KeyError(f"processes missing from stages.PROCESSES: {sorted(unmapped)}")
    return trace.with_columns(
        stage=pl.col("process").replace_strict(STAGE_OF),
        granularity=pl.col("process").replace_strict(GRANULARITY_OF),
        tool=pl.col("process").replace_strict(TOOL_NAMES, default=None, return_dtype=pl.String),
    )


def parse_tag(trace: pl.DataFrame) -> pl.DataFrame:
    """Extract assembler and sample identifiers from task tags."""
    tokens = sorted(set(TAG_OF_ASSEMBLER.values()), key=len, reverse=True)
    return trace.with_columns(
        tag_assembler=pl.col("tag").str.extract(rf"^({'|'.join(tokens)})-"),
        # Remove per-bin and sequencing-run suffixes.
        sample=(
            pl.col("tag")
            .str.split("-")
            .list.last()
            .str.replace(r"\..*$", "")
            # CONCOCT and SemiBin2 number their bins with an underscore instead of a dot.
            .str.replace(r"_\d+$", "")
            .str.replace(r"_run\d+.*$", "")
        ),
    )
