"""Per-process resource reference table."""

from __future__ import annotations

import polars as pl

from .stages import STAGE_ORDER


def _range(column: str, decimals: int) -> pl.Expr:
    """Format the observed span of a metric as ``min-max``."""
    return (
        pl.col(column).min().round(decimals).cast(pl.String)
        + "-"
        + pl.col(column).max().round(decimals).cast(pl.String)
    ).alias(f"{column}_range")


def process_resources(trace: pl.DataFrame) -> pl.DataFrame:
    """Summarise observed resource use per process, pooling datasets and assemblers.

    Values are what tasks consumed, not what the pipeline reserved for them. Ranges
    span both datasets and every assembly a process ran on, so they describe the
    spread a user should expect rather than the effect of any single input.
    """
    return (
        trace.group_by("stage", "process", "granularity")
        .agg(
            pl.len().alias("n_tasks"),
            pl.col("peak_rss_gb").median().round(2),
            _range("peak_rss_gb", 2),
            pl.col("realtime_h").median().round(2),
            _range("realtime_h", 2),
            pl.col("workdir_gb").median().round(2),
            _range("workdir_gb", 2),
            pl.col("cpu_hours").median().round(2),
            _range("cpu_hours", 2),
            pl.col("cpu_hours").sum().round(1).alias("cpu_hours_total"),
        )
        .with_columns(stage=pl.col("stage").cast(pl.Enum(STAGE_ORDER)))
        .select(
            "stage",
            "process",
            "granularity",
            "n_tasks",
            "cpu_hours",
            "cpu_hours_range",
            "cpu_hours_total",
            "peak_rss_gb",
            "peak_rss_gb_range",
            "realtime_h",
            "realtime_h_range",
            "workdir_gb",
            "workdir_gb_range",
        )
        .sort("stage", "cpu_hours_total", descending=[False, True])
    )
