"""Parse nf-core/mag execution traces and reports."""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

_SIZE_UNITS = {"B": 1, "KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40, "PB": 2**50}
_MISSING = ["-", "", "NA"]

#: Extract the required fields without parsing the report's JavaScript payload.
_REPORT_PAYLOAD = re.compile(r"window\.data\s*=\s*(\{.*?\});", re.S)
_REPORT_FIELDS = {
    name: re.compile(rf'"{name}":"(\d+)"')
    for name in ("task_id", "cpus", "attempt")
}


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


def _report_cpus(run_dir: Path) -> pl.DataFrame | None:
    """Read per-task CPU reservations from an execution report."""
    reports = sorted(run_dir.glob("pipeline_info/execution_report_*.html"))
    if not reports:
        return None
    match = _REPORT_PAYLOAD.search(reports[0].read_text(errors="replace"))
    if match is None:
        return None
    payload = match.group(1)
    columns = {name: rx.findall(payload) for name, rx in _REPORT_FIELDS.items()}
    if not columns["task_id"]:
        return None
    if len({len(values) for values in columns.values()}) != 1:
        raise ValueError(f"{run_dir}: execution report fields do not line up")
    return pl.DataFrame(columns).cast(pl.Int64)


def load_run(run_dir: Path) -> pl.DataFrame:
    """Read one run directory (e.g. ``data/maghini_hybrid``)."""
    traces = sorted(run_dir.glob("pipeline_info/execution_trace_*_extended.txt"))
    if len(traces) != 1:
        raise FileNotFoundError(f"{run_dir}: expected 1 extended trace, found {len(traces)}")

    dataset, _, assembly_mode = run_dir.name.partition("_")
    raw = pl.read_csv(traces[0], separator="\t", infer_schema_length=None, null_values=_MISSING)

    parsed = raw.with_columns(
        dataset=pl.lit(dataset),
        assembly_mode=pl.lit(assembly_mode),
        # Split the qualified process name and trailing task tag.
        process=pl.col("name").str.split(" (").list.get(0).str.split(":").list.last(),
        tag=pl.col("name").str.extract(r"\(([^)]*)\)$"),
        realtime_s=_duration_s("realtime"),
        cpu_pct=pl.col("%cpu").str.strip_chars_end("%").cast(pl.Float64),
        peak_rss_b=_size_bytes("peak_rss"),
        wchar_b=_size_bytes("wchar"),
    )

    reserved = _report_cpus(run_dir)
    if reserved is None:
        # Attempts are unknown when the report has no task payload.
        return parsed.with_columns(
            cpus=pl.lit(None, dtype=pl.Int64),
            attempt=pl.lit(None, dtype=pl.Int64),
        )
    return parsed.join(reserved, on="task_id", how="left")


def load_traces(data_dir: str | Path = "../data") -> pl.DataFrame:
    """Load all runs under ``data_dir`` into one table."""
    data_dir = Path(data_dir)
    runs = [p for p in sorted(data_dir.iterdir()) if p.is_dir()]
    trace = _fill_missing_reservations(
        pl.concat([load_run(p) for p in runs], how="diagonal_relaxed")
    )

    return trace.with_columns(
        realtime_h=pl.col("realtime_s") / 3600,
        # CPU time consumed, not reserved capacity.
        cpu_hours=pl.col("cpu_pct") / 100 * pl.col("realtime_s") / 3600,
        cpu_efficiency=pl.col("cpu_pct") / (100 * pl.col("cpus")),
        peak_rss_gb=pl.col("peak_rss_b") / 2**30,
        written_gb=pl.col("wchar_b") / 2**30,
        workdir_gb=pl.col("disk_bytes") / 2**30,
    )


def _fill_missing_reservations(trace: pl.DataFrame) -> pl.DataFrame:
    """Impute missing reservations from recorded first attempts."""
    known = (
        trace.filter(pl.col("cpus").is_not_null() & (pl.col("attempt") == 1))
        .group_by("process")
        .agg(default_cpus=pl.col("cpus").mode().first())
    )
    filled = trace.join(known, on="process", how="left").with_columns(
        cpus=pl.coalesce("cpus", "default_cpus")
    )
    unresolved = filled.filter(pl.col("cpus").is_null())["process"].unique().to_list()
    if unresolved:
        raise ValueError(f"no CPU allocation recorded for: {sorted(unresolved)}")
    return filled.drop("default_cpus")
def sample_counts(trace: pl.DataFrame) -> pl.DataFrame:
    """Count biological samples from preprocessing task tags."""
    return (
        trace.filter(pl.col("process").is_in(["FASTP", "CHOPPER", "FASTQC_RAW", "NANOPLOT_RAW"]))
        .with_columns(sample=pl.col("tag").str.replace(r"_run\d+.*$", ""))
        .group_by("dataset")
        .agg(n_samples=pl.col("sample").n_unique())
        .sort("dataset")
    )
