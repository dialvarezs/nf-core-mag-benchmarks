"""Parse nf-core/mag execution traces and reports."""

from __future__ import annotations

from pathlib import Path

import polars as pl

_SIZE_UNITS = {"B": 1, "KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40, "PB": 2**50}
_MISSING = ["-", "", "NA"]


def _duration_s(column: str) -> pl.Expr:
    """Convert a Nextflow duration column to seconds."""
    value = pl.col(column)

    def part(pattern: str, seconds: float) -> pl.Expr:
        return value.str.extract(pattern).cast(pl.Float64).fill_null(0.0) * seconds

    total = (
        part(r"([\d.]+)d\b", 86400.0)
        + part(r"([\d.]+)h\b", 3600.0)
        + part(r"([\d.]+)m\b", 60.0)
        + part(r"([\d.]+)ms\b", 1e-3)
        + part(r"([\d.]+)s\b", 1.0)
    )
    # Preserve missing durations after summing optional components.
    return pl.when(value.is_null()).then(None).otherwise(total)


def _size_bytes(column: str) -> pl.Expr:
    """Nextflow memory size ("2.1 GB", binary multiples) to bytes."""
    value = pl.col(column)
    return value.str.extract(r"([\d.]+)").cast(pl.Float64) * value.str.extract(
        r"([KMGTP]?B)$"
    ).replace_strict(_SIZE_UNITS, default=None, return_dtype=pl.Float64)


def load_run(run_dir: Path) -> pl.DataFrame:
    """Read one run directory (e.g. ``data/maghini``)."""
    traces = sorted(run_dir.glob("pipeline_info/execution_trace_*_extended.txt"))
    if len(traces) != 1:
        raise FileNotFoundError(f"{run_dir}: expected 1 extended trace, found {len(traces)}")

    raw = pl.read_csv(traces[0], separator="\t", infer_schema_length=None, null_values=_MISSING)

    return raw.with_columns(
        dataset=pl.lit(run_dir.name),
        # Split the qualified process name and trailing task tag.
        process=pl.col("name").str.split(" (").list.get(0).str.split(":").list.last(),
        tag=pl.col("name").str.extract(r"\(([^)]*)\)$"),
        realtime_s=_duration_s("realtime"),
        cpu_pct=pl.col("%cpu").str.strip_chars_end("%").cast(pl.Float64),
        peak_rss_b=_size_bytes("peak_rss"),
        wchar_b=_size_bytes("wchar"),
    )


def load_traces(data_dir: str | Path = "../data") -> pl.DataFrame:
    """Load all runs under ``data_dir`` into one table."""
    data_dir = Path(data_dir)
    runs = [p for p in sorted(data_dir.iterdir()) if p.is_dir()]
    trace = pl.concat([load_run(p) for p in runs], how="diagonal_relaxed")

    return trace.with_columns(
        realtime_h=pl.col("realtime_s") / 3600,
        # CPU time consumed, not reserved capacity.
        cpu_hours=pl.col("cpu_pct") / 100 * pl.col("realtime_s") / 3600,
        peak_rss_gb=pl.col("peak_rss_b") / 2**30,
        written_gb=pl.col("wchar_b") / 2**30,
        workdir_gb=pl.col("disk_bytes") / 2**30,
    )


def sample_counts(trace: pl.DataFrame) -> pl.DataFrame:
    """Count biological samples from preprocessing task tags."""
    return (
        trace.filter(pl.col("process").is_in(["FASTP", "CHOPPER", "FASTQC_RAW", "NANOPLOT_RAW"]))
        .with_columns(sample=pl.col("tag").str.replace(r"_run\d+.*$", ""))
        .group_by("dataset")
        .agg(n_samples=pl.col("sample").n_unique())
        .sort("dataset")
    )
