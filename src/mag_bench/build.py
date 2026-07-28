"""Build every supplementary figure and write it to disk.

Run as ``uv run python -m mag_bench.build`` from the project root, or call
:func:`prepare` from a notebook to get the same frames the figures use.
"""

from __future__ import annotations

import time
from pathlib import Path

import polars as pl

from . import figures as fg
from .aggregate import (
    canonical,
    load_assembly_stats,
    stage_budget,
    stage_totals,
    tool_tasks,
    with_assembly_stats,
)
from .stages import annotate, parse_tag
from .trace import load_traces, sample_counts


def prepare(data_dir: str | Path = "data") -> dict[str, pl.DataFrame]:
    """Parse the traces once and derive every frame the figures need.

    The notebook and the command line both go through here, so the two
    cannot drift apart.
    """
    # the sample count is read off the preprocessing tags, so it is taken
    # before the failures are dropped
    raw = parse_tag(annotate(load_traces(data_dir)))
    samples = sample_counts(raw)
    succeeded = raw.filter(pl.col("status") != "FAILED")
    kept = canonical(succeeded)
    tasks = tool_tasks(kept)
    totals = stage_totals(kept, samples)
    return {
        "trace": kept,
        "totals": totals,
        "budget": stage_budget(totals),
        "tasks": tasks,
        "scaling": with_assembly_stats(tasks, load_assembly_stats(data_dir)),
    }


def build_all(data_dir: str | Path = "data") -> dict[str, object]:
    frames = prepare(data_dir)
    return {
        "S1_stage_compute": fg.fig_stage_compute(frames["budget"], frames["trace"]),
        "S2_assembler_footprint": fg.fig_assembler_footprint(frames["tasks"]),
        "S3_binner_footprint": fg.fig_binner_footprint(frames["tasks"]),
        "S4_storage": fg.fig_storage(frames["totals"]),
        "S5_scaling_assemblers": fg.fig_scaling_assemblers(frames["scaling"]),
        "S6_scaling_binners": fg.fig_scaling_binners(frames["scaling"]),
        "S7_top_processes": fg.fig_top_processes(frames["trace"]),
    }


def save_all(
    figs: dict[str, object], out_dir: str | Path = "figures", dpi: int = 300
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, plot in figs.items():
        path = out_dir / f"{name}.png"
        started = time.perf_counter()
        plot.save(path, dpi=dpi, verbose=False)
        print(f"{path}  ({time.perf_counter() - started:.1f}s)", flush=True)
        written.append(path)
    return written


if __name__ == "__main__":
    # Composed figures ask matplotlib for a canvas as they are assembled,
    # which picks up an interactive backend when a display is present.
    import matplotlib

    matplotlib.use("Agg")

    save_all(build_all())
