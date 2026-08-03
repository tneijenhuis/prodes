"""Plots the timings collected by scripts/benchmark.py.

Kept separate from the benchmark itself so that re-plotting never means
re-running a measurement that can take hours on a large structure. Reads the
published dataset written by scripts/publish_benchmark.py and writes the figures
beside it under docs/benchmark/::

    python scripts/plot_benchmark.py

Reading the published archive rather than the raw benchmark_results/ directory
is deliberate: every committed figure is then reproducible from committed data
alone, by anyone, without access to the machine the numbers were measured on.
"""

import zipfile
from pathlib import Path

import matplotlib
import pandas as pd

# Selected before pyplot is imported, so that the script also runs on a headless
# server where no display is attached.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
PUBLISHED_DIR = REPO_ROOT / "docs" / "benchmark"
DATA_ARCHIVE = PUBLISHED_DIR / "benchmark_data.zip"

TOTAL_TIME_PLOT = PUBLISHED_DIR / "benchmark_total_time.png"
PHASE_PLOT = PUBLISHED_DIR / "benchmark_phases.png"
SPEEDUP_PLOT = PUBLISHED_DIR / "benchmark_speedup.png"
SPEEDUP_BY_SIZE_PLOT = PUBLISHED_DIR / "benchmark_speedup_by_size.png"

# Plotted against the residue count rather than the file size, which is what the
# first version of this script used. File size is what a user has in front of
# them, but it is not what the calculation costs: hydrogens double the size of a
# PDB file and are discarded before any work starts, so 1GDW and 1GDW_h take the
# same time despite one file being twice as large. Residues rather than heavy
# atoms because the two are equivalent for this purpose, at a near constant 7.6
# to 8.1 heavy atoms per residue across the set, and sequence length is the
# number a user already knows about their protein.
X_COLUMN = "residues"
X_LABEL = "residues in structure (all chains)"

# 1GDW_h is 1GDW with hydrogens added. Once they are stripped the two are the
# same structure, produce byte-identical descriptors and land on the same x, so
# plotting both would just overprint one point on another. It stays in the
# published dataset as evidence that hydrogens make no difference; it is only
# the figures it is left out of.
DUPLICATE_STRUCTURES = {"1GDW_h"}

# The public name of each implementation, as written by publish_benchmark.py.
FORK = "datacatalysis/prodes"
UPSTREAM = "tneijenhuis/prodes"

# One series is one implementation, on one machine, at one worker count. The
# hostname has to be part of that key even though the published data currently
# holds a single machine: the same implementation at the same worker count runs
# several times faster on one machine than another, and averaging two machines
# together would invent a number that neither of them produced.
SERIES_COLUMNS = ["implementation", "hostname", "n_workers"]

# The house style is green, never the
# matplotlib default blue and orange, and the figures here need at most three
# series: upstream reads as the muted grey-teal baseline and the fork's two
# worker counts climb from dark teal to the primary brand green.
BRAND_GREEN = "#059456"
DARK_TEAL = "#285855"
MUTED_TEAL = "#628482"
BORDER = "#D8E8E1"

# The phase figure is the one deliberate exception to the house palette. Six
# phases in six shades of the brand green were not distinguishable, so it uses a
# qualitative palette instead, where the colours carry no meaning beyond telling
# the series apart. Every series also gets its own marker, so the figure still
# reads when printed in greyscale.
PHASE_COLOURS = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B"]
PHASE_MARKERS = ["o", "s", "^", "D", "v", "P"]

# Type sizes for the speedup figure, which is the one that gets put in a slide
# or skimmed on a phone rather than read on screen beside its explanation, and
# which carries three bars and their labels where the others carry a dozen
# points. The remaining figures keep the matplotlib defaults, where the smaller
# type is what lets them hold that much detail.
SPEEDUP_TITLE_SIZE = 15
SPEEDUP_LABEL_SIZE = 13
SPEEDUP_TICK_SIZE = 12


def use_house_style():
    """Applies the house style to every plot in this script.

    Roboto is the house font and is used when it happens to be installed, which
    on most machines it is not, so the first available fallback is taken rather
    than letting matplotlib warn once per text object and silently substitute.
    """

    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next(name for name in ("Roboto", "Arial", "DejaVu Sans") if name in available)

    plt.rcParams.update(
        {
            "font.family": family,
            "text.color": DARK_TEAL,
            "axes.labelcolor": DARK_TEAL,
            "axes.edgecolor": MUTED_TEAL,
            "xtick.color": DARK_TEAL,
            "ytick.color": DARK_TEAL,
            "axes.titlecolor": DARK_TEAL,
            "grid.color": BORDER,
            "grid.linestyle": ":",
            "grid.linewidth": 0.5,
            "grid.alpha": 0.6,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_results():
    """Reads and concatenates every CSV in the published benchmark archive.

    Raises:
        FileNotFoundError: if the archive has not been written yet.
    """

    if not DATA_ARCHIVE.exists():
        raise FileNotFoundError(f"{DATA_ARCHIVE} not found, run scripts/publish_benchmark.py first")

    with zipfile.ZipFile(DATA_ARCHIVE) as archive:
        names = sorted(name for name in archive.namelist() if name.endswith(".csv"))
        frames = [pd.read_csv(archive.open(name)) for name in names]

    print(f"read {len(names)} CSV file(s) from {DATA_ARCHIVE.name}: {', '.join(names)}")

    return pd.concat(frames, ignore_index=True)


def machine_caption(results):
    """Returns a one line description of the machine the timings were measured on."""

    machines = results[["hostname", "logical_cpus"]].drop_duplicates().sort_values("hostname")

    return ", ".join(f"{row.hostname}, {row.logical_cpus} logical CPUs" for row in machines.itertuples())


def series_label(implementation, hostname, workers, with_hostname):
    """Returns the legend label for one implementation, machine and worker count.

    Args:
        implementation: public implementation name.
        hostname: anonymised machine name.
        workers: number of worker processes the run resolved to.
        with_hostname: whether to name the machine, which is only worth the
            space in the legend once more than one machine has been measured.
    """

    if implementation == UPSTREAM:
        # Upstream has no parallel module at all, so quoting a worker count
        # would suggest a choice that does not exist there.
        label = f"{UPSTREAM}, serial"
    else:
        label = f"{implementation}, {workers} worker" + ("" if workers == 1 else "s")

    return f"{label} ({hostname})" if with_hostname else label


def bar_label(implementation, hostname, workers, with_hostname):
    """Returns the axis label for one bar of the speedup figure.

    Deliberately not series_label, which is written for a reader who already
    knows what prodes is and what a worker is. This figure is the one shown to
    someone who knows neither, so the repository names come off the bars
    entirely, leaving which codebase and how many CPUs, and the surrounding text
    says which repositories those are.
    """

    if implementation == UPSTREAM:
        # No core count: the original has no parallel mode, so the only thing a
        # number here could do is imply it was run with a handicap.
        label = "original codebase"
    else:
        label = f"improved codebase, {workers} CPU" + ("" if workers == 1 else "s")

    return f"{label}\n({hostname})" if with_hostname else label


def series_colour(implementation, workers):
    """Returns the house colour for one series.

    Upstream is the baseline and takes the muted teal; the fork darkens to the
    brand green as the worker count rises, so that the fastest series is the one
    carrying the brand colour.
    """

    if implementation == UPSTREAM:
        return MUTED_TEAL

    return DARK_TEAL if workers == 1 else BRAND_GREEN


def format_seconds(seconds):
    """Returns a short label for a time in seconds, at a sensible precision."""

    return f"{seconds:.0f}s" if seconds >= 10 else f"{seconds:.1f}s"


def format_computation_time(seconds):
    """Returns a time in minutes, spelled out over two lines.

    Used on the speedup figure only. Minutes rather than seconds because the
    point of that figure is the difference between a wait you sit through and a
    wait you plan a coffee around, and 5146 is a number a reader has to divide
    before it means anything. A run under a minute keeps one decimal, so the
    fastest bar reads 0.5 minutes rather than rounding to nothing.
    """

    minutes = seconds / 60
    value = f"{minutes:.0f}" if minutes >= 1 else f"{minutes:.1f}"
    unit = "minute" if value == "1" else "minutes"

    return f"{value} {unit}\ncomputation time"


def format_speedup(speedup):
    """Returns a speedup rounded to the nearest whole number, labelled.

    The decimal is false precision on a figure like this: it is one measurement
    on one machine, so 27.9x is not meaningfully different from 28x, and the
    digit after the point invites a reader to treat it as exact. Nothing
    coarser than the nearest whole number, though. An earlier version rounded to
    the nearest ten above 10x, which reads better but no longer matches what
    anyone dividing two columns of the published CSVs gets, and a headline
    number that cannot be recalculated is not worth the tidier figure.
    """

    return f"{round(speedup):.0f}x speedup"


def series_totals(results):
    """Returns the total time per structure and series, averaging repeated measurements.

    Averaging matters off Linux, where every requested worker count resolves to
    a single serial run, so the same series is measured more than once.
    """

    totals = results[results["phase"] == "TOTAL"]
    grouped = totals.groupby(SERIES_COLUMNS + ["structure", X_COLUMN, "heavy_atoms", "pdb_kb"], as_index=False)

    return grouped["seconds"].mean()


def plotted_totals(results):
    """Returns the per structure totals that go into the figures."""

    totals = series_totals(results)

    return totals[~totals["structure"].isin(DUPLICATE_STRUCTURES)]


def shared_structures(totals):
    """Returns the structures that every series measured.

    Comparisons are restricted to these: a series that was only ever run on the
    small structures would otherwise look fastest purely for having skipped the
    expensive ones, and a speedup can only be quoted where both implementations
    have a number.

    Raises:
        ValueError: if no single structure was measured by every series.
    """

    measured = totals.groupby(SERIES_COLUMNS)["structure"].apply(set)
    shared = set.intersection(*measured)
    if not shared:
        raise ValueError("no structure was measured by every series, so they cannot be compared")

    return shared


def fastest_series(totals):
    """Returns the series key with the lowest total time over the shared structures."""

    common = totals[totals["structure"].isin(shared_structures(totals))]
    ranked = common.groupby(SERIES_COLUMNS)["seconds"].sum().sort_values()

    return ranked.index[0]


def speedups(totals):
    """Returns the speedup of every series over upstream, per structure.

    Args:
        totals: per structure totals from series_totals.

    Returns:
        DataFrame of the shared structures only, with a speedup column giving
        how many times faster each series is than upstream on that structure.

    Raises:
        ValueError: if upstream was not measured on the shared structures, since
            there is then nothing to compare against.
    """

    common = totals[totals["structure"].isin(shared_structures(totals))].copy()

    baseline = common[common["implementation"] == UPSTREAM].set_index("structure")["seconds"]
    if baseline.empty:
        raise ValueError(f"{UPSTREAM} was not measured on any shared structure, so there is no baseline")

    common["speedup"] = common["structure"].map(baseline) / common["seconds"]

    return common


def apply_scale(axes):
    """Sets both axes to a logarithmic scale.

    Logarithmic on both: the times span from under two seconds to well over an
    hour, and the structures from 60 to 823 residues, so a linear axis either
    way flattens most of the data into one corner of the figure. An earlier
    version offered a linear x as well, which showed the four smaller
    structures crowded into the left hand edge and was dropped.
    """

    axes.set_xscale("log")
    axes.set_yscale("log")


def annotate_points(axes, group, y_column, formatter):
    """Writes the value of every point above its marker."""

    for _, row in group.iterrows():
        axes.annotate(
            formatter(row[y_column]),
            (row[X_COLUMN], row[y_column]),
            textcoords="offset points",
            xytext=(0, 9),
            ha="center",
            fontsize=7.5,
        )


def finish(figure, axes, title, results, output_path, which_grid="both", with_machine=True, title_size=11):
    """Applies the shared title, grid and legend treatment and writes the figure.

    Args:
        title_size: point size of the title. The default suits the figures that
            are read on screen next to their explanation; the speedup figure
            raises it, see plot_speedup.
        with_machine: whether to name the machine under the title. On by
            default, since a timing without the machine it was measured on is
            not a number anyone can check. The speedup figure turns it off: it
            plots a ratio between two runs on the same machine, which the
            machine caption does not qualify, and it is the one figure shown to
            people for whom a logical CPU count is noise.
    """

    axes.set_title(f"{title}\n{machine_caption(results)}" if with_machine else title, fontsize=title_size)
    axes.grid(True, which=which_grid, axis="both")
    axes.set_axisbelow(True)

    for side in ("top", "right"):
        axes.spines[side].set_visible(False)

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    print(f"wrote {output_path}")


def plot_total_time(results, output_path):
    """Plots total calculation time against structure size, one line per series."""

    totals = plotted_totals(results)
    with_hostname = totals["hostname"].nunique() > 1

    figure, axes = plt.subplots(figsize=(8, 5.5))

    for key, group in totals.groupby(SERIES_COLUMNS):
        implementation, _, workers = key
        group = group.sort_values(X_COLUMN)
        axes.plot(
            group[X_COLUMN],
            group["seconds"],
            marker="o",
            color=series_colour(implementation, workers),
            label=series_label(*key, with_hostname),
        )
        annotate_points(axes, group, "seconds", format_seconds)

    apply_scale(axes)
    axes.set_xlabel(X_LABEL)
    axes.set_ylabel("total calculation time (s)")
    axes.legend(fontsize="small", frameon=False)
    finish(figure, axes, "prodes calculation time by structure size", results, output_path)


def plot_phase_breakdown(results, output_path):
    """Plots each phase against structure size for the fastest series measured.

    Shows which phase dominates as structures grow, which is what decides where
    any further optimisation is worth spending effort.
    """

    implementation, hostname, workers = fastest_series(plotted_totals(results))

    selected = results[
        (results["implementation"] == implementation)
        & (results["hostname"] == hostname)
        & (results["n_workers"] == workers)
        & (results["phase"] != "TOTAL")
        & (~results["structure"].isin(DUPLICATE_STRUCTURES))
    ]
    per_phase = selected.groupby(["phase", X_COLUMN], as_index=False)["seconds"].mean()

    figure, axes = plt.subplots(figsize=(8, 5.5))

    for number, (phase, group) in enumerate(per_phase.groupby("phase")):
        group = group.sort_values(X_COLUMN)
        axes.plot(
            group[X_COLUMN],
            group["seconds"],
            marker=PHASE_MARKERS[number % len(PHASE_MARKERS)],
            markersize=5,
            color=PHASE_COLOURS[number % len(PHASE_COLOURS)],
            label=phase,
        )

    apply_scale(axes)
    axes.set_xlabel(X_LABEL)
    axes.set_ylabel("phase time (s)")
    axes.legend(fontsize="small", frameon=False)
    label = series_label(implementation, hostname, workers, with_hostname=False)
    finish(figure, axes, f"prodes time per phase: {label}", results, output_path)


def plot_speedup(results, output_path):
    """Plots the speedup over upstream on the largest structure measured, as bars.

    The largest structure is the one worth quoting: on a small structure the
    fixed costs dominate and the parallel phases have too little to do to pay
    for the workers.
    """

    measured = speedups(plotted_totals(results))
    largest = measured.loc[measured[X_COLUMN].idxmax()]
    measured = measured[measured["structure"] == largest["structure"]].sort_values("speedup")
    with_hostname = measured["hostname"].nunique() > 1

    labels = [bar_label(row.implementation, row.hostname, row.n_workers, with_hostname) for row in measured.itertuples()]
    colours = [series_colour(row.implementation, row.n_workers) for row in measured.itertuples()]

    figure, axes = plt.subplots(figsize=(8, 5.5))
    bars = axes.bar(labels, measured["speedup"], color=colours, width=0.55)

    for bar, (_, row) in zip(bars, measured.iterrows(), strict=True):
        axes.annotate(
            f"{format_speedup(row['speedup'])}\n{format_computation_time(row['seconds'])}",
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            textcoords="offset points",
            xytext=(0, 5),
            ha="center",
            fontsize=SPEEDUP_LABEL_SIZE,
        )

    # Deliberately no note about how much of the gain came from the extra cores.
    # It is the sort of thing that reads as clutter on a figure whose whole job
    # is one comparison, and the speedup by size figure already carries it.

    # Headroom for the value labels above the tallest bar, which tight_layout
    # does not reserve because annotations are not part of the data limits.
    # Three lines of label now rather than two, so more of it than before.
    axes.set_ylim(0, measured["speedup"].max() * 1.30)
    axes.set_ylabel("Speed-up in computational time (x)", fontsize=SPEEDUP_LABEL_SIZE)
    axes.tick_params(labelsize=SPEEDUP_TICK_SIZE)

    # The PDB code names the structure precisely and means nothing to anyone who
    # does not work with the PDB, where the size of the protein is the thing that
    # makes the number below meaningful. The code stays on the other figures,
    # which are read by people comparing structures.
    title = f"Example speed-up for {largest[X_COLUMN]}-residue protein"
    finish(figure, axes, title, results, output_path, which_grid="major", with_machine=False, title_size=SPEEDUP_TITLE_SIZE)


def parallel_gain_note(measured):
    """Returns a sentence on how the gain from the extra workers varies with size, or None.

    This is the half of the picture the total speedup hides. The speedup over
    upstream is flat because two effects cancel: the vectorisation gain shrinks
    with size while the parallel gain grows, since a bigger structure gives the
    workers more to do per unit of pool startup cost.
    """

    fork = measured[measured["implementation"] == FORK]
    serial = fork[fork["n_workers"] == 1].set_index("structure")
    parallel = fork[fork["n_workers"] > 1]

    if serial.empty or parallel.empty:
        return None

    gains = parallel.assign(gain=parallel["structure"].map(serial["seconds"]) / parallel["seconds"]).sort_values(X_COLUMN)
    smallest, largest = gains.iloc[0], gains.iloc[-1]
    workers = int(largest["n_workers"])

    return (
        f"the parallel gain does grow with size: {workers} workers over 1 worker is\n"
        f"{smallest['gain']:.1f}x at {smallest[X_COLUMN]} residues, {largest['gain']:.1f}x at {largest[X_COLUMN]}"
    )


def plot_speedup_by_size(results, output_path):
    """Plots the speedup over upstream against structure size, one line per fork series."""

    measured = speedups(plotted_totals(results))
    measured = measured[measured["implementation"] == FORK]
    with_hostname = measured["hostname"].nunique() > 1

    figure, axes = plt.subplots(figsize=(8, 5.5))

    for key, group in measured.groupby(SERIES_COLUMNS):
        implementation, _, workers = key
        group = group.sort_values(X_COLUMN)
        axes.plot(
            group[X_COLUMN],
            group["speedup"],
            marker="o",
            color=series_colour(implementation, workers),
            label=series_label(*key, with_hostname),
        )
        annotate_points(axes, group, "speedup", lambda speedup: f"{speedup:.0f}x")

    # The baseline is 1x by construction, so drawing it makes the figure read as
    # "how far above the original", which is the question being asked.
    axes.axhline(1, color=MUTED_TEAL, linewidth=1, linestyle="--")
    axes.annotate(
        f"{UPSTREAM} = 1x",
        (0.99, 1),
        # Pinned to y=1 in data space so the label stays on the line it belongs
        # to, and to the right hand end, which is the corner nothing else uses.
        xycoords=("axes fraction", "data"),
        textcoords="offset points",
        xytext=(0, 5),
        ha="right",
        fontsize=8,
        color=MUTED_TEAL,
    )

    note = parallel_gain_note(measured)
    if note is not None:
        axes.text(
            0.03,
            0.06,
            note,
            transform=axes.transAxes,
            va="bottom",
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "#EEF6F3", "edgecolor": BORDER},
        )

    apply_scale(axes)
    axes.set_xlabel(X_LABEL)
    axes.set_ylabel(f"speedup over {UPSTREAM} (x)")

    # Headroom above the highest point for its value label, and the legend clear
    # of both the series and the note box.
    axes.set_ylim(top=measured["speedup"].max() * 1.8)
    axes.legend(fontsize="small", frameon=False, loc="center right")

    finish(figure, axes, "prodes speedup over the original, by structure size", results, output_path)


def print_summary(results):
    """Prints the total time per structure and series as a plain text table."""

    totals = series_totals(results).sort_values([X_COLUMN, "structure"] + SERIES_COLUMNS)
    with_hostname = totals["hostname"].nunique() > 1

    print(f"\n{'structure':<12}{'res':>6}{'atoms':>7}{'KB':>7}  {'series':<38}{'seconds':>10}")
    for _, row in totals.iterrows():
        label = series_label(row["implementation"], row["hostname"], row["n_workers"], with_hostname)
        print(f"{row['structure']:<12}{row[X_COLUMN]:>6}{row['heavy_atoms']:>7}{row['pdb_kb']:>7.0f}  {label:<38}{row['seconds']:>10.2f}")


def main():
    """Loads the published results, prints a summary and writes every figure."""

    use_house_style()

    results = load_results()
    print_summary(results)

    plot_total_time(results, TOTAL_TIME_PLOT)
    plot_phase_breakdown(results, PHASE_PLOT)
    plot_speedup(results, SPEEDUP_PLOT)
    plot_speedup_by_size(results, SPEEDUP_BY_SIZE_PLOT)


if __name__ == "__main__":
    main()
