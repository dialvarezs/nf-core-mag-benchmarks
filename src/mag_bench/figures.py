"""Supplementary resource-use figures for the nf-core/mag v5 manuscript."""

from __future__ import annotations

from dataclasses import dataclass

import plotnine as p9
import polars as pl
from plotnine.composition import Beside

from .stages import ASSEMBLERS, BINNERS, GRANULARITY_ORDER, STAGE_ORDER

#: Minimum observations required for a boxplot.
MIN_N_FOR_BOX = 8
MAX_JITTER_POINTS_PER_GROUP = 500
JITTER_SEED = 0
STAGE_JITTER_SIZE = 0.9
STAGE_JITTER_ALPHA = 0.35


@dataclass(repr=False)
class _WeightedBeside(Beside):
    """Place plots side by side with explicit relative widths."""

    width_ratios: tuple[float, ...]

    def _create_gridspec(self, figure, nest_into):
        from plotnine._mpl.gridspec import p9GridSpec

        self.gridspec = p9GridSpec(
            self.nrow,
            self.ncol,
            figure,
            width_ratios=self.width_ratios,
            nest_into=nest_into,
        )

    def draw(self, *, show: bool = False):
        from contextlib import nullcontext
        from warnings import warn

        from plotnine._mpl.layout_manager import PlotnineCompositionLayoutEngine
        from plotnine._mpl.layout_manager._layout_tree import LayoutTree
        from plotnine._mpl.layout_manager._spaces import LayoutSpaces
        from plotnine.exceptions import PlotnineWarning

        class WeightedLayoutEngine(PlotnineCompositionLayoutEngine):
            def execute(self, figure):
                renderer = figure._get_renderer()
                lookup_spaces = {}
                with getattr(renderer, "_draw_disabled", nullcontext)():
                    for plotspec in self.composition.plotspecs:
                        lookup_spaces[plotspec.plot] = LayoutSpaces(plotspec.plot)

                tree = LayoutTree.create(self.composition, lookup_spaces)
                tree.align_axis_titles()
                tree.align()
                plot_widths = tree.plot_widths
                non_panel_widths = [
                    plot_width - panel_width
                    for plot_width, panel_width in zip(plot_widths, tree.panel_widths)
                ]
                panel_unit = (sum(plot_widths) - sum(non_panel_widths)) / sum(self.composition.width_ratios)
                widths = [
                    margin + panel_unit * ratio
                    for margin, ratio in zip(non_panel_widths, self.composition.width_ratios)
                ]
                self.composition.gridspec.set_width_ratios(widths)

                for plot, spaces in lookup_spaces.items():
                    params = spaces.get_gridspec_params()
                    if not params.valid:
                        warn(
                            "The figure is too small to contain all plots.",
                            PlotnineWarning,
                            stacklevel=2,
                        )
                        break
                    plot.facet._panels_gridspec.layout(params)
                    spaces.items._adjust_positions(spaces)

        figure = super().draw(show=False)
        figure.set_layout_engine(WeightedLayoutEngine(self))
        figure.canvas.draw()
        if show:
            figure.show()
        return figure


DATASET_COLOURS = {"maghini": "#1F78B4", "zymo": "#E8A33D"}
DATASET_LABELS = {"maghini": "maghini (human gut, 10 samples)", "zymo": "zymo (mock, 1 sample)"}
GRANULARITY_COLOURS = {
    "per sample": "#8FBFE0",
    "per assembly": "#2E9E7C",
    "per assembly x binner": "#E8A33D",
    "per bin": "#C4553B",
    "once per run": "#B0AFAC",
}

#: Shared metric labels.
METRIC_LABELS = {
    "peak_rss_gb": "Peak memory (GB)",
    "realtime_h": "Wall-clock time (h)",
    "cpu_hours": "CPU core-hours",
    "workdir_gb": "Work-directory size (GB)",
}


def _labels(*columns: str) -> dict[str, str]:
    """Return labels in panel order."""
    return {column: METRIC_LABELS[column] for column in columns}


#: Shared dimensions and typography.
FIGURE_WIDTH = 13
BASE_SIZE = 11


def theme_mag(height: float = 7) -> p9.theme:
    """Return the shared figure theme."""
    return p9.theme_bw(base_size=BASE_SIZE) + p9.theme(
        figure_size=(FIGURE_WIDTH, height),
        legend_position="bottom",
        panel_grid_minor=p9.element_blank(),
        panel_grid_major_y=p9.element_line(colour="#E6E6E6", size=0.4),
        panel_grid_major_x=p9.element_line(colour="#E6E6E6", size=0.4),
        panel_border=p9.element_rect(colour="#4D4D4D", size=0.5),
        strip_background=p9.element_rect(fill="#F0F0F0", colour="#4D4D4D", size=0.5),
        strip_text=p9.element_text(weight="bold", size=BASE_SIZE),
        plot_title=p9.element_text(weight="bold", size=BASE_SIZE + 3, ha="left"),
        plot_subtitle=p9.element_text(size=BASE_SIZE, colour="#4D4D4D", ha="left"),
        plot_caption=p9.element_text(size=BASE_SIZE - 2, colour="#4D4D4D", ha="right"),
        # Prevent long categorical labels from shifting titles off-canvas.
        plot_title_position="plot",
        plot_caption_position="plot",
        legend_key=p9.element_blank(),
        legend_title=p9.element_text(weight="bold"),
        axis_title=p9.element_text(size=BASE_SIZE),
    )


def _dataset_scale():
    return p9.scale_colour_manual(values=DATASET_COLOURS, labels=DATASET_LABELS, name="Dataset")


def si_labels(breaks) -> list[str]:
    """Format axis breaks with SI suffixes."""
    labels = []
    for value in breaks:
        if value is None:
            labels.append("")
            continue
        magnitude, suffix = next(
            (m, s)
            for m, s in ((1e12, "T"), (1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""))
            if abs(value) >= m or m == 1
        )
        scaled = value / magnitude
        text = f"{scaled:g}" if scaled == int(scaled) else f"{scaled:.10g}"
        labels.append(f"{text}{suffix}")
    return labels


def gb_labels(breaks) -> list[str]:
    """Format GB-valued breaks as byte-size units."""
    return si_labels([None if value is None else value * 1e9 for value in breaks])


def _log_y():
    # plotnine mutates scales when composing plots.
    return p9.scale_y_log10(labels=si_labels)


def _big_groups(df: pl.DataFrame, group_cols: list[str]) -> pl.DataFrame:
    """Return groups large enough for a boxplot."""
    return df.group_by(group_cols).len().filter(pl.col("len") >= MIN_N_FOR_BOX).select(group_cols)


def _sample_for_jitter(
    df: pl.DataFrame,
    group_cols: list[str],
    max_points: int = MAX_JITTER_POINTS_PER_GROUP,
    seed: int = JITTER_SEED,
) -> pl.DataFrame:
    """Sample points within visual groups for compact vector output."""
    return df.group_by(group_cols, maintain_order=True).map_groups(
        lambda group: group.sample(n=min(group.height, max_points), seed=seed)
    )


def to_long(df: pl.DataFrame, metrics: dict[str, str], id_vars: list[str]) -> pl.DataFrame:
    """Melt metrics, preserve panel order, and remove invalid log values."""
    return (
        df.select(*id_vars, *metrics)
        .unpivot(index=id_vars, variable_name="metric", value_name="value")
        .with_columns(metric=pl.col("metric").replace_strict(metrics).cast(pl.Enum(list(metrics.values()))))
        .drop_nulls("value")
        .filter(pl.col("value") > 0)
    )


#: Fixed baseline for bars on logarithmic axes.
BAR_BASELINE = 0.01

PER_JOB_PANEL = "CPU core-hours per job"
JOBS_PANEL = "Jobs launched"
TOTAL_PANEL = "Total core-hours"
TOTAL_DISK_PANEL = "Total left on disk"


def _summary_bar_panel(rows: pl.DataFrame, stages: list[str], panel_label: str, labels) -> p9.ggplot:
    """Plot one pooled summary bar per stage."""
    summary = rows.with_columns(
        panel=pl.lit(panel_label),
        label_x=pl.max_horizontal("value", "high"),
    ).with_columns(
        text=pl.when(pl.col("value") >= 10)
        .then(pl.col("value").round(0).cast(pl.Int64).cast(pl.String))
        .otherwise(pl.col("value").round(1).cast(pl.String))
    )
    return (
        p9.ggplot(summary, p9.aes("stage_label", "value"))
        # geom_col starts at zero, which is undefined on a log axis.
        + p9.geom_segment(
            p9.aes(x="stage_label", xend="stage_label", y=BAR_BASELINE, yend="value"),
            # Neutral colour because the bars pool datasets.
            colour="#A6A6A6",
            size=8,
            inherit_aes=False,
        )
        + p9.geom_linerange(
            data=summary.filter(~pl.col("single_dataset")),
            mapping=p9.aes(x="stage_label", ymin="low", ymax="high"),
            colour="#3D3D3D",
            size=0.6,
            inherit_aes=False,
        )
        + p9.geom_text(
            p9.aes(y="label_x", label="text"),
            ha="left",
            size=7.5,
            colour="#4D4D4D",
            nudge_y=0.02,
        )
        + p9.scale_x_discrete(limits=stages)
        # Explicit limits censor data; the segment baseline sets the floor.
        + p9.scale_y_log10(labels=labels, expand=(0.02, 0, 0.25, 0))
        + p9.facet_wrap("panel")
        + p9.coord_flip()
        + p9.labs(x="", y="per sample")
        # The last composed plot sets the final canvas size.
        + theme_mag()
        + p9.theme(axis_text_y=p9.element_blank(), axis_ticks_major_y=p9.element_blank())
    )


def fig_stage_compute(budget: pl.DataFrame, canonical_trace: pl.DataFrame):
    """Plot per-job and total compute use by stage."""
    labelled = budget.with_columns(
        stage_label=pl.col("stage")
        + pl.when(pl.col("single_dataset")).then(pl.lit(" *")).otherwise(pl.lit(""))
    )
    label_of = dict(zip(labelled["stage"], labelled["stage_label"]))
    stages = [label_of[s] for s in reversed(STAGE_ORDER) if s in label_of]

    per_job = canonical_trace.filter(pl.col("cpu_hours") > 0).with_columns(
        panel=pl.lit(PER_JOB_PANEL),
        stage_label=pl.col("stage").replace_strict(label_of),
    )
    big = _big_groups(per_job, ["stage_label"])
    jitter = _sample_for_jitter(per_job, ["stage_label"])

    left = (
        p9.ggplot(per_job, p9.aes("stage_label", "cpu_hours"))
        + p9.geom_boxplot(
            data=per_job.join(big, on="stage_label", how="semi"),
            outlier_alpha=0,
            size=0.4,
            fill="#F2F2F2",
            colour="#4D4D4D",
            width=0.55,
        )
        + p9.geom_jitter(
            data=jitter,
            mapping=p9.aes(colour="dataset"),
            size=STAGE_JITTER_SIZE,
            alpha=STAGE_JITTER_ALPHA,
            width=0.2,
            height=0,
            random_state=JITTER_SEED,
        )
        + p9.scale_x_discrete(limits=stages)
        + _log_y()
        + _dataset_scale()
        + p9.guides(colour=p9.guide_legend(override_aes={"size": 3, "alpha": 1}))
        + p9.facet_wrap("panel")
        + p9.coord_flip()
        + p9.labs(x="Analysis stage", y="")
        + theme_mag()
    )

    jobs_rows = labelled.select(
        "stage_label",
        value="jobs_per_sample",
        low="jobs_low",
        high="jobs_high",
        single_dataset="single_dataset",
    )
    total_rows = labelled.select(
        "stage_label",
        value="value_per_sample",
        low="low",
        high="high",
        single_dataset="single_dataset",
    )
    return _WeightedBeside(
        [
            left,
            _summary_bar_panel(jobs_rows, stages, JOBS_PANEL, si_labels),
            _summary_bar_panel(total_rows, stages, TOTAL_PANEL, si_labels),
        ],
        width_ratios=(4, 1, 1),
    )


def fig_storage(budget: pl.DataFrame, storage: pl.DataFrame, storage_budget: pl.DataFrame) -> p9.ggplot:
    """Plot data written and retained per stage and sample."""
    long = to_long(
        storage.filter(pl.col("stage").cast(pl.String) != "Assembly QC (ALE)"),
        {"written_gb": "Bytes written by the jobs", "workdir_gb": "Left on disk (work directory)"},
        ["dataset", "stage", "sample"],
    )
    present = set(long["stage"].cast(pl.String))
    stages = [s for s in reversed(STAGE_ORDER) if s in present]
    big = _big_groups(long, ["stage", "metric"])

    def distribution(metric_label: str, *, show_axis: bool) -> p9.ggplot:
        data = long.filter(pl.col("metric") == metric_label)
        jitter = _sample_for_jitter(data, ["stage"])
        plot = (
            p9.ggplot(data, p9.aes("stage", "value"))
            + p9.geom_boxplot(
                data=data.join(big, on=["stage", "metric"], how="semi"),
                outlier_alpha=0,
                size=0.4,
                fill="#F2F2F2",
                colour="#4D4D4D",
                width=0.55,
            )
            + p9.geom_jitter(
                data=jitter,
                mapping=p9.aes(colour="dataset"),
                size=STAGE_JITTER_SIZE,
                alpha=STAGE_JITTER_ALPHA,
                width=0.18,
                height=0,
                random_state=JITTER_SEED,
            )
            + p9.scale_x_discrete(limits=stages)
            + p9.scale_y_log10(labels=gb_labels)
            + p9.facet_wrap("metric")
            + _dataset_scale()
            + p9.coord_flip()
            + p9.labs(x="Analysis stage" if show_axis else "", y="")
            + theme_mag()
        )
        if not show_axis:
            plot = plot + p9.theme(
                legend_position="none", axis_text_y=p9.element_blank(), axis_ticks_major_y=p9.element_blank()
            )
        return plot

    written = distribution("Bytes written by the jobs", show_axis=True)
    left_on_disk = distribution("Left on disk (work directory)", show_axis=False)

    in_stages = pl.col("stage").cast(pl.String).is_in(stages)
    jobs_rows = budget.filter(in_stages).select(
        stage_label="stage",
        value="jobs_per_sample",
        low="jobs_low",
        high="jobs_high",
        single_dataset="single_dataset",
    )
    total_disk_rows = storage_budget.filter(in_stages).select(
        stage_label="stage",
        value="value_per_sample",
        low="low",
        high="high",
        single_dataset="single_dataset",
    )
    bars = _summary_bar_panel(jobs_rows, stages, JOBS_PANEL, si_labels) | _summary_bar_panel(
        total_disk_rows, stages, TOTAL_DISK_PANEL, gb_labels
    )
    return written | left_on_disk | bars


def fig_tool_footprint(tasks: pl.DataFrame) -> p9.ggplot:
    """Plot assembler and binner resource distributions."""
    # Exclude known retries without dropping the run that lacks attempt data.
    labelled = tasks.filter(
        pl.col("tool").is_in(ASSEMBLERS + BINNERS),
        ~(pl.col("attempt") > 1).fill_null(False),
    ).with_columns(
        tool=pl.col("tool").cast(pl.String),
        group=pl.when(pl.col("tool").cast(pl.String).is_in(ASSEMBLERS))
        .then(pl.lit("Assemblers"))
        .otherwise(pl.lit("Binners")),
    )
    metrics = _labels("peak_rss_gb", "realtime_h", "cpu_hours")
    long = to_long(labelled, metrics, ["dataset", "tool", "group"])
    big = _big_groups(long, ["tool", "metric", "dataset"])

    def panel(
        group: str, tools: list[str], metric: str, *, row_label: str, show_strip: bool, legend: bool
    ) -> p9.ggplot:
        data = long.filter((pl.col("group") == group) & (pl.col("metric") == metric))
        jitter = _sample_for_jitter(data, ["tool"])
        plot = (
            p9.ggplot(data, p9.aes("tool", "value"))
            + p9.geom_boxplot(
                data=data.join(big, on=["tool", "metric", "dataset"], how="semi"),
                outlier_alpha=0,
                size=0.4,
                fill="#F2F2F2",
                colour="#4D4D4D",
                width=0.55,
            )
            + p9.geom_jitter(
                data=jitter,
                mapping=p9.aes(colour="dataset"),
                size=1.5,
                alpha=0.75,
                width=0.18,
                height=0,
                random_state=JITTER_SEED,
            )
            + p9.scale_x_discrete(limits=tools)
            + _log_y()
            + _dataset_scale()
            + p9.facet_wrap("metric")
            + p9.coord_flip()
            + p9.labs(x=row_label, y="")
            + theme_mag(height=8)
        )
        if not legend:
            plot = plot + p9.theme(legend_position="none")
        if not show_strip:
            plot = plot + p9.theme(strip_background=p9.element_blank(), strip_text=p9.element_blank())
        return plot

    metric_names = list(metrics.values())
    top = (
        panel(
            "Assemblers",
            ASSEMBLERS,
            metric_names[0],
            row_label="Assemblers",
            show_strip=True,
            legend=False,
        )
        | panel("Assemblers", ASSEMBLERS, metric_names[1], row_label="", show_strip=True, legend=False)
        | panel("Assemblers", ASSEMBLERS, metric_names[2], row_label="", show_strip=True, legend=False)
    )
    bottom = (
        panel("Binners", BINNERS, metric_names[0], row_label="Binners", show_strip=False, legend=True)
        | panel("Binners", BINNERS, metric_names[1], row_label="", show_strip=False, legend=False)
        | panel("Binners", BINNERS, metric_names[2], row_label="", show_strip=False, legend=False)
    )
    return top / bottom


def _scaling(
    tasks_with_stats: pl.DataFrame,
    tools: list[str],
    x_col: str,
    x_label: str,
    scales: str,
) -> p9.ggplot:
    """Plot resource use against tool input size."""
    long = to_long(
        tasks_with_stats.filter(pl.col("tool").is_in(tools)),
        _labels("peak_rss_gb", "cpu_hours"),
        ["dataset", "tool", x_col],
    )

    return (
        p9.ggplot(long, p9.aes(x_col, "value", colour="dataset"))
        + p9.geom_point(size=1.6, alpha=0.8)
        + p9.geom_smooth(
            mapping=p9.aes(x_col, "value"),
            method="lm",
            colour="#4D4D4D",
            fill="#CCCCCC",
            size=0.6,
            alpha=0.3,
            inherit_aes=False,
        )
        + p9.scale_x_log10(labels=si_labels)
        + _log_y()
        + p9.facet_grid("metric ~ tool", scales=scales)
        + _dataset_scale()
        + p9.labs(x=x_label, y="")
        + theme_mag(height=6)
    )


def fig_scaling_assemblers(tasks_with_stats: pl.DataFrame) -> p9.ggplot:
    return _scaling(
        tasks_with_stats,
        ASSEMBLERS,
        x_col="input_gbp",
        x_label="Sequencing input (Gbp)",
        scales="free",
    )


def fig_scaling_binners(tasks_with_stats: pl.DataFrame) -> p9.ggplot:
    return _scaling(
        tasks_with_stats,
        BINNERS,
        x_col="assembly_length",
        x_label="Assembly size (bp)",
        scales="free_y",
    )


def rank_by_best_resource(canonical_trace: pl.DataFrame, columns: list[str], top_n: int) -> pl.DataFrame:
    """Rank processes by their best rank across the selected resources."""
    ranks = None
    for column in columns:
        per_metric = (
            canonical_trace.filter(pl.col(column) > 0)
            .group_by("process")
            .agg(median=pl.col(column).median())
            .sort("median", descending=True)
            .with_row_index(f"rank_{column}", offset=1)
            .select("process", f"rank_{column}")
        )
        ranks = (
            per_metric if ranks is None else ranks.join(per_metric, on="process", how="full", coalesce=True)
        )

    rank_columns = [f"rank_{column}" for column in columns]
    # Missing metrics rank last rather than first.
    last = canonical_trace["process"].n_unique()
    return (
        ranks.with_columns([pl.col(c).fill_null(last) for c in rank_columns])
        .with_columns(best_rank=pl.min_horizontal(rank_columns), rank_sum=pl.sum_horizontal(rank_columns))
        .sort("best_rank", "rank_sum")
        .head(top_n)
    )


def fig_top_processes(canonical_trace: pl.DataFrame, top_n: int = 15):
    """Plot the top processes across runtime, memory, and disk use."""
    metrics = _labels("realtime_h", "peak_rss_gb", "workdir_gb")
    ranked = rank_by_best_resource(canonical_trace, list(metrics), top_n)
    labelled = canonical_trace.join(ranked.select("process"), on="process", how="semi").with_columns(
        label=pl.format("{}  ({} jobs)", "process", pl.len().over("process")),
        granularity=pl.col("granularity").cast(pl.Enum(GRANULARITY_ORDER)),
    )
    order = (
        labelled.select("process", "label")
        .unique()
        .join(ranked.select("process", "best_rank", "rank_sum"), on="process")
        .sort("best_rank", "rank_sum", descending=True)["label"]
        .to_list()
    )

    def panel(column: str, scale, *, first: bool = False) -> p9.ggplot:
        """Plot one resource panel."""
        plot = (
            p9.ggplot(
                to_long(labelled, _labels(column), ["label", "granularity"]),
                p9.aes("label", "value", fill="granularity"),
            )
            + p9.geom_boxplot(outlier_size=0.4, outlier_alpha=0.3, size=0.4)
            + p9.scale_x_discrete(limits=order)
            + scale
            + p9.scale_fill_manual(values=GRANULARITY_COLOURS, name="This step runs", drop=False)
            + p9.facet_wrap("metric", nrow=1)
            + p9.coord_flip()
            + p9.labs(x="", y="")
            + theme_mag(height=6.5)
        )
        if first:
            return plot
        return plot + p9.theme(
            axis_text_y=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
            legend_position="none",
        )

    return (
        panel("realtime_h", _log_y(), first=True)
        | panel("peak_rss_gb", _log_y())
        | panel("workdir_gb", p9.scale_y_continuous(labels=si_labels))
    )
