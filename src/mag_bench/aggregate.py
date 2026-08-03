"""Aggregate task traces into figure-ready tables."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .stages import ASSEMBLERS, BINNERS, STAGE_ORDER, TAG_OF_ASSEMBLER


def polish_only_stages(trace: pl.DataFrame) -> list[str]:
    """Return stages unique to polish runs."""
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
    """Aggregate stage metrics per sample."""
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


def storage_per_sample(canonical_trace: pl.DataFrame) -> pl.DataFrame:
    """Aggregate storage by dataset, stage, and biological sample."""
    return (
        canonical_trace.filter(pl.col("granularity") != "once per run")
        .group_by("dataset", "stage", "sample")
        .agg(
            written_gb=pl.col("written_gb").sum(),
            workdir_gb=pl.col("workdir_gb").sum(),
        )
        .with_columns(stage=pl.col("stage").cast(pl.Enum(STAGE_ORDER)))
    )


#: Dataset overrides for stages with non-comparable runs.
STAGE_SOURCE = {
    # maghini did not run Prokka.
    "Annotation (Prodigal/Prokka)": "zymo",
    # zymo used the older per-bin DAS Tool implementation.
    "Bin refinement (DAS Tool)": "maghini",
}


def stage_budget(totals: pl.DataFrame, metric: str = "cpu_hours") -> pl.DataFrame:
    """Pool per-sample stage metrics across eligible datasets."""
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
            metric_total=pl.col(f"{metric}_total").sum(),
            n_tasks=pl.col("n_tasks").sum(),
            low=pl.col(f"{metric}_per_sample").min(),
            high=pl.col(f"{metric}_per_sample").max(),
            jobs_low=pl.col("jobs_one_dataset").min(),
            jobs_high=pl.col("jobs_one_dataset").max(),
        )
        .with_columns(
            value_per_sample=pl.col("metric_total") / pl.col("n_samples"),
            jobs_per_sample=pl.col("n_tasks") / pl.col("n_samples"),
            single_dataset=pl.col("n_datasets") == 1,
        )
        .with_columns(
            pct=100 * pl.col("value_per_sample") / pl.col("value_per_sample").sum()
        )
        .sort("value_per_sample", descending=True)
    )


def tool_tasks(canonical_trace: pl.DataFrame) -> pl.DataFrame:
    """Return assembler and binner tasks for footprint figures."""
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
        assembly_length=pl.col("Total length"),
    )


def with_assembly_stats(trace: pl.DataFrame, stats: pl.DataFrame) -> pl.DataFrame:
    """Attach assembly size to traceable assembler and binner tasks."""
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

    # Fail if a plotted tool cannot be matched to assembly statistics.
    orphaned = (
        joined.filter(pl.col("tool").is_not_null() & pl.col("assembly_length").is_null())["tool"]
        .unique()
        .to_list()
    )
    if orphaned:
        raise ValueError(f"no assembly stats matched for: {sorted(orphaned)}")
    return joined


#: Sequencing yield from the published dataset tables, in Gbp.
READ_YIELD_GBP = {
    ("zymo", "zymo_fecal"): {"sr": 15.60, "lr": 11.08},
    ("maghini", "D01"): {"sr": 9.44, "lr": 6.98},
    ("maghini", "D02"): {"sr": 6.89, "lr": 10.50},
    ("maghini", "D03"): {"sr": 15.81, "lr": 3.82},
    ("maghini", "D04"): {"sr": 18.96, "lr": 9.68},
    ("maghini", "D05"): {"sr": 12.01, "lr": 13.31},
    ("maghini", "D06"): {"sr": 14.26, "lr": 13.67},
    ("maghini", "D07"): {"sr": 31.69, "lr": 6.71},
    ("maghini", "D08"): {"sr": 10.15, "lr": 12.31},
    ("maghini", "D09"): {"sr": 9.26, "lr": 9.24},
    ("maghini", "D10"): {"sr": 9.76, "lr": 4.12},
}


def with_read_yield(tasks_with_stats: pl.DataFrame) -> pl.DataFrame:
    """Attach the sequencing input consumed by each assembler stage."""
    yields = pl.DataFrame(
        [
            {"dataset": dataset, "sample": sample, "sr_gbp": gbp["sr"], "lr_gbp": gbp["lr"]}
            for (dataset, sample), gbp in READ_YIELD_GBP.items()
        ]
    )
    return tasks_with_stats.join(yields, on=["dataset", "sample"], how="left").with_columns(
        input_gbp=pl.when(pl.col("stage") == "Short-read assembly")
        .then(pl.col("sr_gbp"))
        .when(pl.col("stage") == "Long-read assembly")
        .then(pl.col("lr_gbp"))
        .when(pl.col("stage") == "Hybrid assembly")
        .then(pl.col("sr_gbp") + pl.col("lr_gbp"))
        .otherwise(None)
    )
