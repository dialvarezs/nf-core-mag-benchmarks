"""Estimate the cost of simpler pipeline configurations from the benchmark traces."""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from .stages import BINNER_TOOLS, GRANULARITY_OF, STAGE_ORDER, TAG_OF_ASSEMBLER

#: Every process a binner needs, including its preparation and extraction steps.
BINNER_PROCESSES = {
    "MetaBAT2": ["METABAT2_METABAT2"],
    "MaxBin2": ["MAXBIN2", "ADJUST_MAXBIN2_EXT"],
    "CONCOCT": [
        "CONCOCT_CUTUPFASTA",
        "CONCOCT_CONCOCTCOVERAGETABLE",
        "CONCOCT_CONCOCT",
        "CONCOCT_MERGECUTUPCLUSTERING",
        "CONCOCT_EXTRACTFASTABINS",
    ],
    "SemiBin2": ["SEMIBIN_SINGLEEASYBIN"],
    "MetaBinner": [
        "METABINNER_TOOSHORT",
        "METABINNER_KMER",
        "METABINNER_METABINNER",
        "METABINNER_BINS",
    ],
    "COMEBin": ["COMEBIN_RUNCOMEBIN"],
}

#: Task tags name the binner that produced the bin set.
TAG_OF_BINNER = {tag for tag, _ in BINNER_TOOLS.values()}

#: Granularities whose tasks belong to one binner.
_PER_BINNER = ("per assembly x binner", "per bin")


@dataclass(frozen=True)
class Config:
    """A pipeline configuration to cost out, as a subset of the benchmarked run."""

    name: str
    assembler: str
    binners: tuple[str, ...]
    stages: tuple[str, ...]

    def __post_init__(self):
        unknown = set(self.stages) - set(STAGE_ORDER)
        if unknown:
            raise KeyError(f"{self.name}: unknown stages {sorted(unknown)}")
        if self.assembler not in TAG_OF_ASSEMBLER:
            raise KeyError(f"{self.name}: unknown assembler {self.assembler}")
        unknown = set(self.binners) - set(BINNER_PROCESSES)
        if unknown:
            raise KeyError(f"{self.name}: unknown binners {sorted(unknown)}")


def _tag_binner() -> pl.Expr:
    return pl.col("tag").str.extract(rf"-({'|'.join(sorted(TAG_OF_BINNER))})-")


def select_tasks(trace: pl.DataFrame, config: Config) -> pl.DataFrame:
    """Return the tasks a configuration would have launched."""
    other_assemblers = [p for p in TAG_OF_ASSEMBLER if p != config.assembler]
    other_binners = [p for b, ps in BINNER_PROCESSES.items() if b not in config.binners for p in ps]
    per_binner = [p for p, grain in GRANULARITY_OF.items() if grain in _PER_BINNER]

    # Tasks either carry the chosen assembler in their tag, are its own process, or are
    # shared by the whole run (preprocessing and reporting).
    right_assembler = (pl.col("tag_assembler") == TAG_OF_ASSEMBLER[config.assembler]) | (
        pl.col("process") == config.assembler
    )
    return trace.filter(
        pl.col("stage").is_in(config.stages),
        ~pl.col("process").is_in(other_assemblers + other_binners),
        right_assembler | pl.col("tag_assembler").is_null(),
        ~pl.col("process").is_in(per_binner) | _tag_binner().is_in(config.binners),
    )


def _totals(tasks: pl.DataFrame, samples: pl.DataFrame, group: list[str]) -> pl.DataFrame:
    return (
        tasks.group_by(group)
        .agg(
            n_tasks=pl.len(),
            cpu_hours=pl.col("cpu_hours").sum(),
            workdir_gb=pl.col("workdir_gb").sum(),
            peak_rss_gb=pl.col("peak_rss_gb").max(),
        )
        .join(samples, on="dataset")
        .with_columns(
            jobs_per_sample=(pl.col("n_tasks") / pl.col("n_samples")).round(1),
            cpu_hours_per_sample=(pl.col("cpu_hours") / pl.col("n_samples")).round(1),
            workdir_gb_per_sample=(pl.col("workdir_gb") / pl.col("n_samples")).round(1),
            peak_rss_gb=pl.col("peak_rss_gb").round(1),
        )
    )


def estimate(trace: pl.DataFrame, samples: pl.DataFrame, configs: list[Config]) -> pl.DataFrame:
    """Cost each configuration per dataset."""
    return pl.concat(
        _totals(select_tasks(trace, c), samples, ["dataset"])
        .with_columns(config=pl.lit(c.name))
        .select(
            "config",
            "dataset",
            "jobs_per_sample",
            "cpu_hours_per_sample",
            "workdir_gb_per_sample",
            "peak_rss_gb",
        )
        for c in configs
    ).sort("config", "dataset")


def estimate_by_stage(trace: pl.DataFrame, samples: pl.DataFrame, config: Config) -> pl.DataFrame:
    """Break one configuration down by analysis stage, pooling the datasets."""
    tasks = select_tasks(trace, config)
    pooled = _totals(tasks, samples, ["dataset", "stage"])
    return (
        pooled.group_by("stage")
        .agg(
            jobs_per_sample=pl.col("n_tasks").sum() / pl.col("n_samples").sum(),
            cpu_hours_per_sample=pl.col("cpu_hours").sum() / pl.col("n_samples").sum(),
            workdir_gb_per_sample=pl.col("workdir_gb").sum() / pl.col("n_samples").sum(),
            peak_rss_gb=pl.col("peak_rss_gb").max(),
        )
        .with_columns(
            pl.col("jobs_per_sample").round(1),
            pl.col("cpu_hours_per_sample").round(1),
            pl.col("workdir_gb_per_sample").round(1),
            stage=pl.col("stage").cast(pl.Enum(STAGE_ORDER)),
        )
        .sort("cpu_hours_per_sample", descending=True)
    )
