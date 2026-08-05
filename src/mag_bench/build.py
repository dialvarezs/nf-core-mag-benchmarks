"""Build and save the supplementary figures."""

from __future__ import annotations

import time
from pathlib import Path

import polars as pl

from . import figures as fg
from .aggregate import (
    load_assembly_stats,
    stage_budget,
    stage_totals,
    storage_per_sample,
    tool_tasks,
    with_assembly_stats,
    with_read_yield,
)
from .stages import annotate, parse_tag
from .trace import load_traces, sample_counts


def prepare(data_dir: str | Path = "data") -> dict[str, pl.DataFrame]:
    """Prepare all data frames used by the figures."""
    # Failed preprocessing tasks are still needed to count samples.
    raw = parse_tag(annotate(load_traces(data_dir)))
    samples = sample_counts(raw)
    kept = raw.filter(pl.col("status") != "FAILED")
    tasks = tool_tasks(kept)
    totals = stage_totals(kept, samples)
    return {
        "trace": kept,
        "samples": samples,
        "budget": stage_budget(totals),
        "storage_budget": stage_budget(totals, metric="workdir_gb"),
        "storage": storage_per_sample(kept),
        "tasks": tasks,
        "scaling": with_read_yield(with_assembly_stats(tasks, load_assembly_stats(data_dir))),
    }


def build_all(data_dir: str | Path = "data") -> dict[str, object]:
    frames = prepare(data_dir)
    return {
        "S1_stage_compute": fg.fig_stage_compute(frames["budget"], frames["trace"]),
        "S2_storage": fg.fig_storage(frames["storage"], frames["storage_budget"]),
        "S3_tool_footprint": fg.fig_tool_footprint(frames["tasks"]),
        "S4_top_processes": fg.fig_top_processes(frames["trace"]),
        "S5_scaling_assemblers": fg.fig_scaling_assemblers(frames["scaling"]),
        "S6_scaling_binners": fg.fig_scaling_binners(frames["scaling"]),
    }


def save_all(
    figs: dict[str, object], out_dir: str | Path = "figures", format: str = "png", dpi: int = 300
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, plot in figs.items():
        path = out_dir / f"{name}.{format}"
        started = time.perf_counter()
        plot.save(path, dpi=dpi, verbose=False)
        print(f"{path}  ({time.perf_counter() - started:.1f}s)", flush=True)
        written.append(path)
    return written


if __name__ == "__main__":
    # Use a non-interactive backend for command-line rendering.
    import matplotlib

    matplotlib.use("Agg")

    save_all(build_all())
