"""Parse nf-core/mag Nextflow execution traces into a tidy table.

The pipeline was run four times (two datasets x two assembly modes). Each
run directory under ``data/`` holds a ``pipeline_info/`` folder with an
``execution_trace_*_extended.txt``: the stock Nextflow trace plus
work-directory size, file counts and symlink counts.

Nextflow prints human-readable durations ("4m 3s") and memory sizes
("2.1 GB", binary multiples), so the raw columns need parsing before any
arithmetic. Tasks that failed carry no metrics at all -- Nextflow writes
"-" for %cpu, peak_rss and wchar -- so they are dropped by default.

**CPU allocations.** The trace has no ``cpus`` column, but the HTML
execution report next to it embeds the same task list *with* the cores
Nextflow reserved, per task and per retry attempt. That matters: the
pipeline scales resources with ``task.attempt``, so a retried metaSPAdes
job holds 20 or 30 cores rather than 10, and reading the allocation off
the process name alone would misreport both the footprint labels and the
CPU-efficiency figure.
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

_SIZE_UNITS = {"B": 1, "KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40, "PB": 2**50}
_MISSING = ["-", "", "NA"]

#: The report embeds its task list as a JavaScript object literal, which is
#: not valid JSON (it escapes forward slashes and single quotes) and runs to
#: 17 MB. Only three integer fields are needed, so they are scanned out
#: directly -- both cheaper and free of the JSON-vs-JS escaping problem.
#: Verified against a full parse of all four reports: identical values.
_REPORT_PAYLOAD = re.compile(r"window\.data\s*=\s*(\{.*?\});", re.S)
_REPORT_FIELDS = {
    name: re.compile(rf'"{name}":"(\d+)"')
    for name in ("task_id", "cpus", "attempt")
}


def _duration_s(column: str) -> pl.Expr:
    """Nextflow duration ("1h 2m 3s", "187ms") to seconds.

    The word boundaries are load-bearing: without them the minute and
    second terms would also match inside "187ms".
    """
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
    # fill_null(0) above would otherwise turn a missing duration into 0 s
    return pl.when(value.is_null()).then(None).otherwise(total)


def _size_bytes(column: str) -> pl.Expr:
    """Nextflow memory size ("2.1 GB", binary multiples) to bytes."""
    value = pl.col(column)
    return value.str.extract(r"([\d.]+)").cast(pl.Float64) * value.str.extract(
        r"([KMGTP]?B)$"
    ).replace_strict(_SIZE_UNITS, default=None, return_dtype=pl.Float64)


def _report_cpus(run_dir: Path) -> pl.DataFrame | None:
    """Cores reserved per task, read from the HTML execution report.

    Returns ``None`` when the report carries no per-task payload: Nextflow
    drops it for runs above ~10 000 tasks, which is the case for one of
    the four runs here.
    """
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
        # "NFCORE_MAG:MAG:BINNING:METABAT2_METABAT2 (D01)" -> process + tag
        process=pl.col("name").str.split(" (").list.get(0).str.split(":").list.last(),
        tag=pl.col("name").str.extract(r"\(([^)]*)\)$"),
        realtime_s=_duration_s("realtime"),
        cpu_pct=pl.col("%cpu").str.strip_chars_end("%").cast(pl.Float64),
        peak_rss_b=_size_bytes("peak_rss"),
        wchar_b=_size_bytes("wchar"),
    )

    reserved = _report_cpus(run_dir)
    if reserved is None:
        # attempt stays null rather than defaulting to 1: this run's retries
        # are indistinguishable from first attempts here, and calling them
        # first attempts would pair a retry's peak with the base reservation
        return parsed.with_columns(
            cpus=pl.lit(None, dtype=pl.Int64),
            attempt=pl.lit(None, dtype=pl.Int64),
        )
    return parsed.join(reserved, on="task_id", how="left")


def load_traces(data_dir: str | Path = "../data") -> pl.DataFrame:
    """Load every run under ``data_dir`` into one tidy table.

    Failed tasks are kept: Nextflow records no metrics for them, so they
    contribute only nulls, but the caller needs them to count samples.
    """
    data_dir = Path(data_dir)
    runs = [p for p in sorted(data_dir.iterdir()) if p.is_dir()]
    trace = _fill_missing_reservations(
        pl.concat([load_run(p) for p in runs], how="diagonal_relaxed")
    )

    return trace.with_columns(
        realtime_h=pl.col("realtime_s") / 3600,
        # actual CPU time consumed, not the allocation
        cpu_hours=pl.col("cpu_pct") / 100 * pl.col("realtime_s") / 3600,
        cpu_efficiency=pl.col("cpu_pct") / (100 * pl.col("cpus")),
        peak_rss_gb=pl.col("peak_rss_b") / 2**30,
        written_gb=pl.col("wchar_b") / 2**30,
        workdir_gb=pl.col("disk_bytes") / 2**30,
    )


def _fill_missing_reservations(trace: pl.DataFrame) -> pl.DataFrame:
    """Give the run whose report has no payload the usual reservation.

    Every process in that run also ran in a run that does have a payload,
    so the first-attempt reservation can be read off those instead of being
    transcribed from the pipeline config by hand.
    """
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
    """Number of biological samples per dataset, read off the preprocessing tags.

    Short-read tasks are tagged ``<sample>_run<n>`` and long-read tasks
    ``<sample>``, so the run suffix has to be stripped before counting.
    """
    return (
        trace.filter(pl.col("process").is_in(["FASTP", "CHOPPER", "FASTQC_RAW", "NANOPLOT_RAW"]))
        .with_columns(sample=pl.col("tag").str.replace(r"_run\d+.*$", ""))
        .group_by("dataset")
        .agg(n_samples=pl.col("sample").n_unique())
        .sort("dataset")
    )
