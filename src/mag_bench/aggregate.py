"""Turn the per-task trace into the tables the figures are drawn from.

Two things have to happen before any number is comparable across the four
runs:

**Run attribution.** Each dataset was processed twice -- once assembling
short/hybrid/long reads (``*_hybrid``) and once re-assembling the polished
long reads (``*_polish``). The second run repeats most of the workflow, so
summing all four runs would double-count. Every stage is therefore taken
from the hybrid run, except those that exist *only* in the polish run.

**Per-sample normalisation.** The two datasets differ in size (maghini has
10 samples, zymo has 1), so raw totals say more about the dataset than
about the pipeline. Stage totals are divided by the number of samples:
"what one sample costs to push through this stage".
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .stages import ASSEMBLERS, BINNERS, STAGE_ORDER, TAG_OF_ASSEMBLER


def polish_only_stages(trace: pl.DataFrame) -> list[str]:
    """Stages that appear in the polish run and nowhere else.

    Derived rather than listed: the rule *is* "whatever the hybrid run did
    not do", so computing it means a future polish-only stage cannot be
    silently dropped from the totals.
    """
    seen = trace.group_by("assembly_mode").agg(stages=pl.col("stage").unique())
    by_mode = dict(zip(seen["assembly_mode"], seen["stages"]))
    return sorted(set(by_mode.get("polish", [])) - set(by_mode.get("hybrid", [])))


def canonical(trace: pl.DataFrame) -> pl.DataFrame:
    """Drop the tasks that would be double-counted across the two runs."""
    polish_only = polish_only_stages(trace)
    from_hybrid = (pl.col("assembly_mode") == "hybrid") & ~pl.col("stage").is_in(polish_only)
    from_polish = (pl.col("assembly_mode") == "polish") & pl.col("stage").is_in(polish_only)
    return trace.filter(from_hybrid | from_polish)


def stage_totals(canonical_trace: pl.DataFrame, samples: pl.DataFrame) -> pl.DataFrame:
    """Per-sample compute of each stage.

    Takes the already-attributed frame from :func:`canonical`.
    """
    totals = ["cpu_hours", "wall_hours", "written_gb", "workdir_gb"]
    return (
        canonical_trace.group_by("dataset", "stage")
        .agg(
            n_tasks=pl.len(),
            cpu_hours_total=pl.col("cpu_hours").sum(),
            wall_hours_total=pl.col("realtime_h").sum(),
            written_gb_total=pl.col("written_gb").sum(),
            workdir_gb_total=pl.col("workdir_gb").sum(),
        )
        .join(samples, on="dataset")
        .with_columns(
            *[
                (pl.col(f"{name}_total") / pl.col("n_samples")).alias(f"{name}_per_sample")
                for name in totals
            ],
        )
        .with_columns(stage=pl.col("stage").cast(pl.Enum(STAGE_ORDER)))
        .sort("stage", "dataset")
    )


#: Stages measured from one dataset only, because the other's number does
#: not describe the pipeline a user would run today.
STAGE_SOURCE = {
    # maghini never ran Prokka, so its annotation cost is Prodigal alone
    "Annotation (Prodigal/Prokka)": "zymo",
    # the zymo runs predate the removal of `.transpose()` in
    # binning_refinement, which submitted one DAS Tool contig2bin job per bin
    # instead of one per binner; maghini ran the batched version now in dev
    "Bin refinement (DAS Tool)": "maghini",
}


def stage_budget(totals: pl.DataFrame) -> pl.DataFrame:
    """One figure per stage per sample, pooled over the datasets that ran it.

    Derived from :func:`stage_totals` rather than re-aggregating the trace,
    so "sum core-hours per stage, divide by sample count" is defined once.

    Pooling happens on the totals, not by averaging the two datasets'
    per-sample figures: maghini contributes 10 samples and zymo 1, so an
    average of the two would weight the mock community tenfold. Stages that
    only one dataset ran are divided by that dataset's sample count.

    ``low``/``high`` carry the per-sample figure each dataset gave on its
    own, so the single number can be reported without hiding that Binning
    differs threefold between them.
    """
    unknown = set(STAGE_SOURCE) - set(totals["stage"].cast(pl.String))
    if unknown:
        raise KeyError(f"STAGE_SOURCE names stages that do not exist: {sorted(unknown)}")

    chosen = totals.filter(
        pl.col("dataset")
        == pl.col("stage")
        .cast(pl.String)
        .replace_strict(STAGE_SOURCE, default=pl.col("dataset"))
    ).with_columns(jobs_one_dataset=pl.col("n_tasks") / pl.col("n_samples"))
    return (
        chosen.group_by("stage")
        .agg(
            datasets=pl.col("dataset").sort().str.join("+"),
            n_datasets=pl.len(),
            n_samples=pl.col("n_samples").sum(),
            cpu_hours=pl.col("cpu_hours_total").sum(),
            n_tasks=pl.col("n_tasks").sum(),
            low=pl.col("cpu_hours_per_sample").min(),
            high=pl.col("cpu_hours_per_sample").max(),
            jobs_low=pl.col("jobs_one_dataset").min(),
            jobs_high=pl.col("jobs_one_dataset").max(),
        )
        .with_columns(
            core_hours_per_sample=pl.col("cpu_hours") / pl.col("n_samples"),
            jobs_per_sample=pl.col("n_tasks") / pl.col("n_samples"),
            # read off the result, not off STAGE_SOURCE: GUNC and GTDB-Tk are
            # also single-dataset, simply because only zymo ran them
            single_dataset=pl.col("n_datasets") == 1,
        )
        .with_columns(
            pct=100 * pl.col("core_hours_per_sample") / pl.col("core_hours_per_sample").sum()
        )
        .sort("core_hours_per_sample", descending=True)
    )


def tool_tasks(canonical_trace: pl.DataFrame) -> pl.DataFrame:
    """One row per task of an assembler or binner, for the footprint figures.

    Each row is one assembly (or one assembly x binner combination), since
    the caller has already applied the run attribution.
    """
    tools = ASSEMBLERS + BINNERS
    return (
        canonical_trace.filter(pl.col("tool").is_in(tools))
        .with_columns(tool=pl.col("tool").cast(pl.Enum(tools)))
        .sort("tool")
    )


def load_assembly_stats(data_dir: str | Path = "../data") -> pl.DataFrame:
    """QUAST metrics per assembly, used to relate resource use to assembly size."""
    data_dir = Path(data_dir)
    stats = pl.concat(
        [
            pl.read_csv(p / "quast_assembly.csv", infer_schema_length=None, null_values=[""])
            for p in sorted(data_dir.iterdir())
            if p.is_dir()
        ],
        how="diagonal_relaxed",
    )
    return stats.select(
        dataset="dataset",
        assembly_mode="experiment",
        sample="sample",
        tag_assembler="assembler",
        n_contigs=pl.col("# contigs (>= 0 bp)"),
    )


def with_assembly_stats(trace: pl.DataFrame, stats: pl.DataFrame) -> pl.DataFrame:
    """Attach assembly size to every task that can be traced back to one assembly.

    Assembler tasks are tagged with the sample only, so the assembler has
    to come from the process name; binner tasks carry both in the tag.
    Tasks that belong to no single assembly (read preprocessing, reporting)
    keep null stats and drop out of the scaling figures.
    """
    keyed = trace.with_columns(
        tag_assembler=pl.coalesce(
            pl.col("tag_assembler"),
            pl.col("process").replace_strict(
                TAG_OF_ASSEMBLER, default=None, return_dtype=pl.String
            ),
        )
    )
    joined = keyed.join(
        stats, on=["dataset", "assembly_mode", "sample", "tag_assembler"], how="left"
    )

    # The scaling captions quote a per-panel n, so a tool silently losing
    # its stats must not pass unnoticed.
    orphaned = (
        joined.filter(pl.col("tool").is_not_null() & pl.col("n_contigs").is_null())["tool"]
        .unique()
        .to_list()
    )
    if orphaned:
        raise ValueError(f"no assembly stats matched for: {sorted(orphaned)}")
    return joined
