"""Supplementary resource-usage figures for the nf-core/mag v5 manuscript.

Every figure is built from the tidy trace produced by :mod:`mag_bench.trace`
and the stage/tool labels from :mod:`mag_bench.stages`. Three conventions
are shared by all of them:

* **Log scales throughout.** Resource use spans four orders of magnitude
  (MetaBAT2 finishes a bin set in ~1 min, COMEBin takes ~5 h), so a linear
  axis collapses everything except the largest tool into the baseline.
* **Boxes only where a distribution exists.** The zymo dataset is a single
  mock sample, so a box drawn over its handful of points would suggest a
  spread that was never measured. Boxes are drawn only above
  :data:`MIN_N_FOR_BOX` observations; the individual points are always shown.
* **No claim in a caption that is not computed.** Sample sizes and retry
  counts are interpolated from the data rather than typed in, because a
  reviewer can check them against the trace we ship.

The summary bars in Figure S1 sit on a log axis. A bar encodes a distance
from zero and a log axis has no zero, so their baseline is pinned explicitly
(:data:`BAR_BASELINE`) and every bar carries its printed value, rather than
leaving a length to speak for a number it cannot.
"""

from __future__ import annotations

import polars as pl
import plotnine as p9

from .stages import ASSEMBLERS, BINNERS, GRANULARITY_ORDER, STAGE_ORDER

#: Below this many observations a box would imply a spread we did not measure.
MIN_N_FOR_BOX = 8

DATASET_COLOURS = {"maghini": "#1F78B4", "zymo": "#E8A33D"}
DATASET_LABELS = {"maghini": "maghini (human gut, 10 samples)", "zymo": "zymo (mock, 1 sample)"}
GRANULARITY_COLOURS = {
    "per sample": "#8FBFE0",
    "per assembly": "#2E9E7C",
    "per assembly x binner": "#E8A33D",
    "per bin": "#C4553B",
    "once per run": "#B0AFAC",
}

#: One name per measured column, so the same quantity is called the same
#: thing in every figure.
METRIC_LABELS = {
    "peak_rss_gb": "Peak memory (GB)",
    "realtime_h": "Wall-clock time (h)",
    "cpu_hours": "CPU core-hours",
    "workdir_gb": "Work-directory size (GB)",
}


def _labels(*columns: str) -> dict[str, str]:
    """Ordered label subset; dict order fixes the panel order."""
    return {column: METRIC_LABELS[column] for column in columns}


#: Every figure is drawn at the same width so that, once scaled to the page,
#: type is the same size in all of them. Only the height varies with content.
FIGURE_WIDTH = 13
BASE_SIZE = 11


def theme_mag(height: float = 7) -> p9.theme:
    """Shared look: light grid, no panel clutter, readable strip labels."""
    base_size = BASE_SIZE
    return p9.theme_bw(base_size=base_size) + p9.theme(
        figure_size=(FIGURE_WIDTH, height),
        legend_position="bottom",
        panel_grid_minor=p9.element_blank(),
        panel_grid_major_y=p9.element_line(colour="#E6E6E6", size=0.4),
        panel_grid_major_x=p9.element_line(colour="#E6E6E6", size=0.4),
        panel_border=p9.element_rect(colour="#4D4D4D", size=0.5),
        strip_background=p9.element_rect(fill="#F0F0F0", colour="#4D4D4D", size=0.5),
        strip_text=p9.element_text(weight="bold", size=base_size),
        plot_title=p9.element_text(weight="bold", size=base_size + 3, ha="left"),
        plot_subtitle=p9.element_text(size=base_size, colour="#4D4D4D", ha="left"),
        plot_caption=p9.element_text(size=base_size - 2, colour="#4D4D4D", ha="right"),
        # anchor title and caption to the figure, not to the panel: with long
        # categorical labels the panel starts far to the right and a
        # panel-anchored title runs off the canvas
        plot_title_position="plot",
        plot_caption_position="plot",
        legend_key=p9.element_blank(),
        legend_title=p9.element_text(weight="bold"),
        axis_title=p9.element_text(size=base_size),
    )


def _dataset_scale():
    return p9.scale_colour_manual(values=DATASET_COLOURS, labels=DATASET_LABELS, name="Dataset")


def si_labels(breaks) -> list[str]:
    """Axis labels without scientific notation.

    A log axis produces breaks spanning milli to millions here, and mizani's
    default renders the large end as "1e6". Readers of a resource figure
    should not have to decode an exponent, so large values get an SI suffix
    ("200k", "1M") and small ones a plain decimal.
    """
    labels = []
    for value in breaks:
        if value is None:
            labels.append("")
            continue
        magnitude, suffix = next(
            (m, s) for m, s in ((1e9, "G"), (1e6, "M"), (1e3, "k"), (1, "")) if abs(value) >= m or m == 1
        )
        scaled = value / magnitude
        # keep one decimal only when the break is not a whole number of units
        text = f"{scaled:g}" if scaled == int(scaled) else f"{scaled:.10g}"
        labels.append(f"{text}{suffix}")
    return labels


def _log_y():
    # a fresh scale per call: plotnine mutates scales when they join a plot
    return p9.scale_y_log10(labels=si_labels)


def _big_groups(df: pl.DataFrame, group_cols: list[str]) -> pl.DataFrame:
    """Keys of the groups large enough to summarise with a box."""
    return (
        df.group_by(group_cols).len().filter(pl.col("len") >= MIN_N_FOR_BOX).select(group_cols)
    )


def to_long(df: pl.DataFrame, metrics: dict[str, str], id_vars: list[str]) -> pl.DataFrame:
    """Melt selected columns into ``metric`` / ``value``, keeping panel order.

    plotnine facets alphabetise their strip labels unless the column is an
    Enum, and the panels here have a deliberate order (memory, then time,
    then CPU), so the mapping is applied as a categorical. Non-positive and
    missing values are dropped, since every figure using this puts the
    value on a log axis.
    """
    return (
        df.select(*id_vars, *metrics)
        .unpivot(index=id_vars, variable_name="metric", value_name="value")
        .with_columns(
            metric=pl.col("metric").replace_strict(metrics).cast(pl.Enum(list(metrics.values())))
        )
        .drop_nulls("value")
        .filter(pl.col("value") > 0)
    )


#: Where the log-scale bars start. Below every value plotted, so no bar is
#: clipped, and fixed rather than inferred so the baseline does not move when
#: the data change.
BAR_BASELINE = 0.01

PER_JOB_PANEL = "CPU core-hours per job"
JOBS_PANEL = "Jobs launched"
TOTAL_PANEL = "Total core-hours"


def fig_stage_compute(budget: pl.DataFrame, canonical_trace: pl.DataFrame):
    """How much compute each stage uses, and which of its two factors drives that.

    A stage's total is (jobs launched) x (core-hours of a job), and the panels are
    those three quantities left to right, so the figure reads as the identity
    itself: assembly is heavy because each job is huge, annotation
    because there are thousands of them.

    The left panel keeps every job and both datasets, since that is a
    distribution. The two summary panels report a single pooled number per
    stage as a bar, with a thin line marking what each dataset gave on its
    own where both ran the stage.

    Built as two plots side by side rather than one three-facet plot: the
    per-job panel carries 15 000 points and needs half the width, and
    plotnine exposes no way to give facets unequal widths. Composing a
    half-width plot with a half-width two-facet plot gives 1/2, 1/4, 1/4.
    Both use the same discrete scale, so the rows line up.
    """
    labelled = budget.with_columns(
        stage_label=pl.col("stage")
        + pl.when(pl.col("single_dataset")).then(pl.lit(" *")).otherwise(pl.lit(""))
    )
    label_of = dict(zip(labelled["stage"], labelled["stage_label"]))
    stages = [label_of[s] for s in reversed(STAGE_ORDER) if s in label_of]
    pooled = budget.filter(~pl.col("single_dataset"))["n_samples"].max()

    per_job = (
        canonical_trace.filter(pl.col("cpu_hours") > 0)
        .with_columns(
            panel=pl.lit(PER_JOB_PANEL),
            stage_label=pl.col("stage").replace_strict(label_of),
        )
    )
    big = _big_groups(per_job, ["stage_label"])

    left = (
        p9.ggplot(per_job, p9.aes("stage_label", "cpu_hours"))
        # one box per stage over both datasets, with every job drawn on top.
        # The stages that run thousands of jobs become a dense band, which is
        # the point: their cost is repetition, not any one job being large.
        + p9.geom_boxplot(
            data=per_job.join(big, on="stage_label", how="semi"),
            outlier_alpha=0,
            size=0.4,
            fill="#F2F2F2",
            colour="#4D4D4D",
            width=0.55,
        )
        + p9.geom_jitter(
            p9.aes(colour="dataset"), size=0.9, alpha=0.35, width=0.2, height=0, random_state=0
        )
        + p9.scale_x_discrete(limits=stages)
        + _log_y()
        + _dataset_scale()
        + p9.guides(colour=p9.guide_legend(override_aes={"size": 3, "alpha": 1}))
        + p9.facet_wrap("panel")
        + p9.coord_flip()
        + p9.labs(
            title="Figure S1.  How much compute each analysis stage uses, and why",
            subtitle=(
                "Total = jobs launched x core-hours per job. Core-hours are (%CPU / 100) x run time;\n"
                "the two right-hand panels pool both datasets into one figure per sample."
            ),
            x="Analysis stage",
            y="",
            caption=(
                "Left: every job, none sampled away. Right: summed core-hours divided by the "
                f"samples that ran the\nstage ({pooled} where both datasets did), so the larger "
                "dataset is not outweighed by the single mock; the\ndark line is the figure each "
                "dataset gave alone. * measured from one dataset only -- annotation from zymo,\n"
                "the only run with Prokka; DAS Tool from maghini, which ran the batched version "
                "now in dev."
            ),
        )
        + theme_mag()
    )

    summary = pl.concat(
        [
            labelled.select(
                "stage_label",
                panel=pl.lit(JOBS_PANEL),
                value="jobs_per_sample",
                low="jobs_low",
                high="jobs_high",
                single_dataset="single_dataset",
            ),
            labelled.select(
                "stage_label",
                panel=pl.lit(TOTAL_PANEL),
                value="core_hours_per_sample",
                low="low",
                high="high",
                single_dataset="single_dataset",
            ),
        ]
    ).with_columns(
        panel=pl.col("panel").cast(pl.Enum([JOBS_PANEL, TOTAL_PANEL])),
        label_x=pl.max_horizontal("value", "high"),
    ).with_columns(
        text=pl.when(pl.col("value") >= 10)
        .then(pl.col("value").round(0).cast(pl.Int64).cast(pl.String))
        .otherwise(pl.col("value").round(1).cast(pl.String))
    )

    right = (
        p9.ggplot(summary, p9.aes("stage_label", "value"))
        # Drawn as thick segments, not geom_col: a column is built from y = 0,
        # which is minus infinity once the axis is logarithmic, so every bar
        # is discarded as non-finite. Starting them at BAR_BASELINE states the
        # baseline the log axis cannot supply.
        + p9.geom_segment(
            p9.aes(x="stage_label", xend="stage_label", y=BAR_BASELINE, yend="value"),
            # grey, not the maghini blue: these bars pool both datasets, and
            # borrowing a colour from the legend would read as one of them
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
        # These bars sit on a log axis, where there is no zero for them to
        # start from, so the baseline is whatever the lower limit happens to
        # be. It is pinned here rather than inferred, and the value is
        # printed on every bar, so a length is never the only thing carrying
        # the number.
        # no explicit limits: an open-ended `limits` censors the data here.
        # The axis floor comes from where the segments start instead.
        + p9.scale_y_log10(labels=si_labels, expand=(0.02, 0, 0.25, 0))
        + p9.facet_wrap("panel", scales="free_x", nrow=1)
        + p9.coord_flip()
        # blank title and subtitle of the same height as the left plot's, so
        # the two panel areas stay vertically aligned
        + p9.labs(x="", y="per sample", title=" ", subtitle=" \n ")
        # the last plot in a composition sets the canvas for the whole thing,
        # so this is the size of the finished figure, not of the right half
        + theme_mag()
        + p9.theme(axis_text_y=p9.element_blank(), axis_ticks_major_y=p9.element_blank())
    )
    return left | right


# ---------------------------------------------------------------------------
# Storage
# ---------------------------------------------------------------------------
def fig_storage(totals: pl.DataFrame) -> p9.ggplot:
    """What each stage leaves behind, per sample."""
    long = to_long(
        totals,
        {
            "written_gb_per_sample": "Bytes written by the jobs",
            "workdir_gb_per_sample": "Left on disk (work directory)",
        },
        ["dataset", "stage"],
    )
    return (
        p9.ggplot(long, p9.aes("stage", "value", colour="dataset"))
        + p9.geom_line(p9.aes(group="stage"), colour="#BDBDBD", size=0.6)
        # the two datasets often land on the same value; without the dodge
        # one marker sits exactly on top of the other and reads as missing
        + p9.geom_point(size=2.8, position=p9.position_dodge(width=0.55))
        + p9.scale_x_discrete(limits=list(reversed(STAGE_ORDER)))
        + _log_y()
        + p9.facet_wrap("metric", scales="free_x", nrow=1)
        + _dataset_scale()
        + p9.coord_flip()
        + p9.labs(
            title="Figure S4.  Storage per analysis stage",
            subtitle=(
                "GB per sample. 'Bytes written' counts every write a job makes, including "
                "temporary files it later deletes, so it exceeds what is left on disk."
            ),
            x="Analysis stage",
            y="GB per sample",
            caption=(
                "ALE is absent because it writes a single score file (< 1 MB); GUNC and GTDB-Tk "
                "were only run on zymo.\nSame run attribution as Figure S1."
            ),
        )
        + theme_mag()
    )
def _footprint(tasks: pl.DataFrame, tools: list[str], title: str, subtitle: str, caption: str):
    subset = tasks.filter(pl.col("tool").is_in(tools))
    # Retries hold more cores than first attempts, so mixing them would put
    # two different allocations under one axis label. Only genuine retries are
    # dropped: one run's report carries no attempt number, and excluding those
    # too would silently remove that whole dataset from the comparison.
    first_try = subset.filter(~(pl.col("attempt") > 1).fill_null(False))
    dropped = subset.height - first_try.height

    # The core counts are what make the wall-clock and CPU panels disagree,
    # so they have to be stated -- but one per tick label turns the axis into
    # a wall of text, so they are listed once in the caption instead. Both
    # this and the efficiencies are read off the data: typed in, they go stale
    # the moment the allocations do.
    summary = (
        first_try.group_by("tool")
        .agg(
            cpus=pl.col("cpus").mode().first(),
            efficiency=pl.col("cpu_efficiency").median(),
        )
        .sort("tool")
        .with_columns(name=pl.col("tool").cast(pl.String).str.replace_all("\n", " "))
    )
    allocation = ", ".join(f"{n} {c}" for n, c in zip(summary["name"], summary["cpus"]))
    worst, best = summary.sort("efficiency").row(0), summary.sort("efficiency").row(-1)
    efficiency = (
        f"Median CPU efficiency spans {worst[3]} at {worst[2]:.0%} of its reservation "
        f"to {best[3]} at {best[2]:.0%}"
    )

    labelled = first_try.with_columns(tool=pl.col("tool").cast(pl.String))
    long = to_long(labelled, _labels("peak_rss_gb", "realtime_h", "cpu_hours"), ["dataset", "tool"])

    return (
        p9.ggplot(long, p9.aes("tool", "value"))
        + p9.geom_boxplot(
            data=long.join(_big_groups(long, ["tool", "metric", "dataset"]),
                           on=["tool", "metric", "dataset"], how="semi"),
            outlier_alpha=0,
            size=0.4,
            fill="#F2F2F2",
            colour="#4D4D4D",
            width=0.55,
        )
        + p9.geom_jitter(
            p9.aes(colour="dataset"), size=1.5, alpha=0.75, width=0.18, height=0, random_state=0
        )
        + p9.scale_x_discrete(limits=tools)
        + _log_y()
        + _dataset_scale()
        + p9.facet_wrap("metric", scales="free_y", nrow=1)
        + p9.labs(
            title=title,
            subtitle=subtitle,
            x="",
            y="",
            caption=(
                f"Cores reserved per job: {allocation}. {efficiency}.\n{caption}\nFailed jobs "
                "carry no metrics in the Nextflow trace and are excluded; so are "
                f"{dropped} retried jobs, which are given more cores than a first attempt."
            ),
        )
        + theme_mag(height=6)
        + p9.theme(axis_text_x=p9.element_text(rotation=30, ha="right"))
    )


def fig_assembler_footprint(tasks: pl.DataFrame) -> p9.ggplot:
    return _footprint(
        tasks,
        ASSEMBLERS,
        title="Figure S2.  Resource footprint per assembler",
        subtitle=(
            "One point per assembly. Wall-clock time is not comparable between assemblers:\n"
            "they are allocated different numbers of cores (listed below), so CPU "
            "core-hours is the like-for-like measure."
        ),
        caption=(
            "Log scales. Boxes summarise the maghini assemblies; zymo is a single assembly per "
            "assembler and is shown as points only."
        ),
    )


def fig_binner_footprint(tasks: pl.DataFrame) -> p9.ggplot:
    return _footprint(
        tasks,
        BINNERS,
        title="Figure S3.  Resource footprint per binner",
        subtitle=(
            "One point per assembly. Binners do not receive the same allocation (listed below),\n"
            "so the CPU core-hours panel, not wall-clock time, is what compares them."
        ),
        caption=(
            "Log scales. Boxes summarise the maghini assemblies; the zymo assemblies are shown "
            "as points only."
        ),
    )


# ---------------------------------------------------------------------------
# Scaling with assembly size
# ---------------------------------------------------------------------------
def _scaling(
    tasks_with_stats: pl.DataFrame,
    tools: list[str],
    title: str,
    subtitle: str,
    caption: str,
    x_label: str,
    scales: str,
) -> p9.ggplot:
    """Resource use against the size of the assembly involved.

    Assemblers and binners are plotted separately on purpose. For a binner
    the contig count is the *input* it has to cluster, so the slope is a
    genuine scaling exponent and a shared x axis compares the six tools
    directly. For an assembler the contig count is its own *output*, and
    each assembler occupies its own narrow band of that axis, so those
    panels get free scales instead.
    """
    subset = tasks_with_stats.filter(pl.col("tool").is_in(tools))
    long = to_long(subset, _labels("peak_rss_gb", "cpu_hours"), ["dataset", "tool", "n_contigs"])
    n_per_panel = long.group_by("tool", "metric").len()["len"].min()

    return (
        p9.ggplot(long, p9.aes("n_contigs", "value", colour="dataset"))
        + p9.geom_point(size=1.6, alpha=0.8)
        + p9.geom_smooth(
            mapping=p9.aes("n_contigs", "value"),
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
        + p9.labs(
            title=title,
            subtitle=subtitle,
            x=x_label,
            y="",
            caption=(
                f"Grey band = linear fit on the log-log scale, both datasets pooled "
                f"(n = {n_per_panel} per panel).\n{caption}"
            ),
        )
        + theme_mag(height=6)
    )


def fig_scaling_assemblers(tasks_with_stats: pl.DataFrame) -> p9.ggplot:
    return _scaling(
        tasks_with_stats,
        ASSEMBLERS,
        title="Figure S5.  Assembler resource use against the size of the assembly produced",
        subtitle=(
            "Both axes log, each panel on its own scale: contig count here is the assembler's\n"
            "own output, and the long-read assemblers produce ~100x fewer contigs."
        ),
        caption="Panels are not on a common x axis: see Figure S6 for the comparable version.",
        x_label="Contigs in the assembly produced",
        scales="free",
    )


def fig_scaling_binners(tasks_with_stats: pl.DataFrame) -> p9.ggplot:
    return _scaling(
        tasks_with_stats,
        BINNERS,
        title="Figure S6.  Binner resource use scales with the complexity of the input assembly",
        subtitle=(
            "Both axes log, shared x axis. Every binner sees the same five assemblies per "
            "sample, so the slope is a scaling exponent: 1 means resource use grows in proportion to "
            "the number of contigs."
        ),
        caption="This is what explains the within-tool spread in Figure S3.",
        x_label="Contigs in the input assembly",
        scales="free_y",
    )


def rank_by_best_resource(
    canonical_trace: pl.DataFrame, columns: list[str], top_n: int
) -> pl.DataFrame:
    """Rank processes by their strongest showing on any single resource.

    Ranking on one resource only measures that resource: ordering by CPU
    drops GTDB-Tk, which is the largest memory consumer, and MetaBinner,
    which leaves the most on disk. A weighted score across the three is not
    an option either -- hours, GB of memory and GB of disk have no common
    unit, so any weighting would be arbitrary. Taking each process's *best*
    rank asks a question that needs no weights: "is this step extreme in
    any one respect?". Ties break on the sum of the ranks, which favours
    the process that is heavy across the board.
    """
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
        ranks = per_metric if ranks is None else ranks.join(per_metric, on="process", how="full", coalesce=True)

    rank_columns = [f"rank_{column}" for column in columns]
    # a process absent from one metric ranks last on it rather than winning by default
    last = canonical_trace["process"].n_unique()
    return (
        ranks.with_columns([pl.col(c).fill_null(last) for c in rank_columns])
        .with_columns(
            best_rank=pl.min_horizontal(rank_columns), rank_sum=pl.sum_horizontal(rank_columns)
        )
        .sort("best_rank", "rank_sum")
        .head(top_n)
    )
def fig_top_processes(canonical_trace: pl.DataFrame, top_n: int = 12):
    """The processes that dominate the bill, across all three resources.

    One ranking drives every panel -- each step's best rank on any single
    resource -- so the rows line up and a process can be read across.

    Time and memory are logarithmic because they span four orders of
    magnitude. Work-directory size is not: most of these steps leave 5-15 GB
    behind, and a log axis compresses that cluster against the right edge to
    make room for two outliers. plotnine applies one transform to all facets
    of a plot, so the linear panel is a separate plot composed alongside.
    """
    metrics = _labels("realtime_h", "peak_rss_gb", "workdir_gb")
    ranked = rank_by_best_resource(canonical_trace, list(metrics), top_n)
    labelled = (
        canonical_trace.join(ranked.select("process"), on="process", how="semi")
        .with_columns(
            label=pl.format("{}  ({} jobs)", "process", pl.len().over("process")),
            granularity=pl.col("granularity").cast(pl.Enum(GRANULARITY_ORDER)),
        )
    )
    order = (
        labelled.select("process", "label")
        .unique()
        .join(ranked.select("process", "best_rank", "rank_sum"), on="process")
        .sort("best_rank", "rank_sum", descending=True)["label"]
        .to_list()
    )

    title = f"Figure S7.  The {top_n} heaviest processes per individual job"
    subtitle = (
        "Ordered by each step's best rank across the three resources, so a step that is "
        "extreme in\nonly one of them still appears. Colour is how often the step runs."
    )
    caption = (
        "Steps that are small per job but launched thousands of times dominate the totals "
        "instead: see Figure S1.\nBoth datasets pooled, hybrid-run attribution as in Figure S1. "
        "Time and memory on log scales, disk linear."
    )

    def panel(column: str, scale, *, first: bool = False, caption: str = "") -> p9.ggplot:
        """One resource. Only the first panel keeps the y labels and legend.

        A composition splits its width evenly, so three single-panel plots
        give three equal thirds. Every plot needs a title and subtitle block
        of the same height -- blank on all but the first -- or the panel
        areas start at different heights and the rows stop lining up.
        """
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
            + p9.labs(
                x="", y="", caption=caption,
                title=title if first else " ",
                subtitle=subtitle if first else " \n ",
            )
            + theme_mag(height=6.5)
        )
        if first:
            return plot
        return plot + p9.theme(
            axis_text_y=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
            legend_position="none",
        )

    # legend under the first panel and caption under the last, so the two
    # do not land on the same strip of the composition and overlap
    return (
        panel("realtime_h", _log_y(), first=True)
        | panel("peak_rss_gb", _log_y())
        | panel("workdir_gb", p9.scale_y_continuous(labels=si_labels), caption=caption)
    )
