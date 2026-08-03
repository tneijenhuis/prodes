"""Times the prodes pipeline phase by phase across a range of structure sizes.

Deliberately not a pytest test. It is far too slow for CI and is meant to be
started by hand, locally or on a server, by someone who expects it to take a
while::

    python scripts/benchmark.py

Results are written to a CSV under benchmark_results/, one row per phase, and
flushed after every structure and worker count, so that a run killed part way
through still leaves usable data. The file is named after the host and is
rewritten from scratch on each run, so move or rename it first if a previous
run's numbers are worth keeping. Every phase is also printed as it finishes, so
even a run that dies mid-structure leaves its numbers in the console log.

The feature values every run produces are also written out, one CSV per
structure and worker count, under benchmark_results/features/. They cost nothing
to save, since the pipeline has already calculated them by the time a run is
timed, and they are what makes this script useful beyond timing: the same
structures run through both implementations give a direct check that the fork
returns what the original returns. Copy the ones worth keeping into
tests/data/reference_output/ and compare against them from a test;
benchmark_results/ itself is git ignored and is overwritten on every run.

The test structures are unpacked once into temp/benchmark_structures, which is
git ignored, and everything after that works with ordinary .pdb files. That
keeps archive handling out of the measured time and lets this script also
measure the original prodes at https://github.com/tneijenhuis/prodes, which has
no zip support but exposes the same phase functions under the same names. Run it
with the interpreter of whichever install is to be measured; the resolved
package path is written into every row, so the two cannot be confused. The
original has no parallelism, so it is timed once rather than at several worker
counts.

Deliberately no pKa file: every run uses the pKa values built into prodes, so
nothing here needs PROPKA or any other external tool installed. Supplying one
would make both the timings and the feature values depend on the version of that
tool, which nobody reproducing this would have. See PKA_FILE below.
"""

import csv
import importlib.util
import os
import platform
import socket
import time
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from prodes.run import (
    calculate_average_chargesurface_grid_features,
    calculate_shell_features,
    calculate_structure_features,
    calculate_surface_grid_features,
    construct_surface_grid,
    prepare_structure,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "tests" / "data"
RESULTS_DIR = REPO_ROOT / "benchmark_results"
FEATURE_DIR = RESULTS_DIR / "features"
STRUCTURE_DIR = REPO_ROOT / "temp" / "benchmark_structures"

# Structures to time, smallest first, so that a run killed early still covers the
# cheap end of the size range. Add the large structures to time on a server here,
# zipped or plain; any path that does not exist is reported and skipped.
STRUCTURES = [
    DATA_DIR / "ARH96693.pdb.zip",  # 479 atoms
    DATA_DIR / "1GDW.pdb.zip",  # 1022 atoms, same protein as 1GDW_h but without hydrogens
    DATA_DIR / "1GDW_h.pdb.zip",  # 1022 atoms, same protein as 1GDW but with hydrogens
    DATA_DIR / "ARH98503.pdb.zip",  # 3106 atoms
    DATA_DIR / "1GPB.pdb.zip",  # 6699 atoms
    # Path("/data/structures/large_complex.pdb"),  # 4000 KB, server only
]

PH = 7
R_PROBE = 1.4
HYDRO_SCALE = "mj_scaled"

# Passed to prepare_structure as the pKa file, and deliberately left as None so
# that prodes uses its own pKa table. Named rather than passed as a bare None at
# the call site because setting it is tempting and would quietly break the
# benchmark: a pKa file has to come from PROPKA, H++ or pypka, so anyone
# reproducing these numbers would need that tool installed, at the same version,
# and both the timings and the feature values would become predictions of that
# tool rather than properties of prodes. Leave it as None.
PKA_FILE = None

# Worker counts to time. 1 is the serial baseline and None uses the
# PRODES_N_WORKERS default of half the logical cores. Drop the serial baseline
# when timing a structure so large that running it twice is not worth the wait.
WORKER_COUNTS = [1, None]

# The original prodes has no parallel module, which is what distinguishes the two
# implementations without having to be told which one is installed.
IS_FORK = importlib.util.find_spec("prodes.parallel") is not None
IMPLEMENTATION = "datacatalysis-prodes" if IS_FORK else "tneijenhuis-prodes"
PRODES_PATH = str(Path(importlib.util.find_spec("prodes").origin).parent)

# Every row carries the settings it was measured under. PRODES_MEM_LIMIT_MB in
# particular changes how large the vectorised chunks are and so how fast the
# calculation runs, which makes a timing meaningless without it. Importing
# prodes.run has already called load_dotenv(), so a value set in .env rather
# than exported into the shell is picked up here too. Both are recorded as the
# empty string when unset, which is not the same as their defaults: the original
# prodes has no memory limit at all, so the column is empty for every row of it.
CSV_COLUMNS = [
    "timestamp",
    "implementation",
    "prodes_path",
    "hostname",
    "platform",
    "logical_cpus",
    "mem_limit_mb",
    "full_features",
    "structure",
    "pdb_kb",
    "atom_records",
    "residues",
    "heavy_atoms",
    "surface_points",
    "requested_workers",
    "n_workers",
    "phase",
    "seconds",
]


def extract_structures(archives, directory: Path):
    """Unpacks the benchmark structures into a working directory and returns their paths.

    Args:
        archives: paths to the structures, each either a .pdb.zip or a plain .pdb.
        directory: where to unpack to; created if it does not exist.

    Returns:
        Paths to plain .pdb files, in the order given, skipping any that are missing.
    """

    directory.mkdir(parents=True, exist_ok=True)
    structures = []

    for archive in archives:
        if not archive.exists():
            print(f"skipping missing structure: {archive}")
            continue

        if archive.suffix != ".zip":
            structures.append(archive)
            continue

        with zipfile.ZipFile(archive) as zipped:
            name = next(member for member in zipped.namelist() if member.endswith(".pdb"))
            zipped.extract(name, directory)
            structures.append(directory / name)

    return structures


def timed(function, *args, **kwargs):
    """Calls a function and returns its result together with the elapsed seconds."""

    start = time.perf_counter()
    result = function(*args, **kwargs)

    return result, time.perf_counter() - start


def structure_size(pdb_file: Path):
    """Returns the size in KB, the ATOM record count and the residue count of a PDB file.

    All three are recorded because they measure different things. Size is what a
    user sees but is the worst predictor of cost, since hydrogens inflate it
    without adding any work. The residue count, summed over every chain, is what
    the calculation really scales with and is what the figures are plotted
    against; the ATOM count is kept so that a protonated file can be told apart
    from an unprotonated one.

    Residues are counted as distinct (chain, residue number) pairs rather than
    with a running counter, so that insertion codes and non-contiguous numbering
    cannot inflate the total.
    """

    lines = [line for line in pdb_file.read_text().splitlines() if line.startswith("ATOM")]
    residues = {(line[21], line[22:27]) for line in lines}

    return pdb_file.stat().st_size / 1024, len(lines), len(residues)


def worker_counts():
    """Returns the worker counts to time, which is a single run for the original prodes."""

    return WORKER_COUNTS if IS_FORK else [None]


def set_workers(workers):
    """Pins PRODES_N_WORKERS for the next run, or clears it to fall back to the default.

    Args:
        workers: worker count to request, or None to use the default.
    """

    if workers is None:
        os.environ.pop("PRODES_N_WORKERS", None)
    else:
        os.environ["PRODES_N_WORKERS"] = str(workers)


def resolved_workers():
    """Returns the number of workers the fork will actually use, or 1 for the original."""

    if not IS_FORK:
        return 1

    from prodes import parallel

    return parallel.n_workers()


def time_pipeline(pdb_file: Path):
    """Runs the whole pipeline on one structure and times each phase separately.

    Mirrors what calculate() does, but times the phases individually and calls
    them with only the arguments the original prodes also accepts, so that the
    two implementations are measured doing the same work. That means the full
    feature set, including the average-charge phase, since the phases default to
    full_features=True and the original has no reduced mode to compare against.
    PRODES_FULL_FEATURES is not consulted here; it belongs to the command line
    entry point. The SASA calculation
    is charged to construct_surface_grid, which is where shrake_rupley runs; by
    the time structure_features is reached the surface is already done.

    The feature dict is returned as well as the timings. calculate() seeds it
    with the structure name and then calls these same phases in this same order,
    so what comes back here has the columns a full_features=True run of
    calculate() would write, and the two can be compared directly.

    Args:
        pdb_file: path to the structure to process.

    Returns:
        features: the calculated features, ID first, in output column order.
        counts: dict with the heavy atom and surface point counts.
        timings: list of (phase, seconds) in execution order, ending with TOTAL.
    """

    features: dict = {}
    timings = []
    counts = {}

    structure, seconds = timed(prepare_structure, str(pdb_file), PKA_FILE)
    features["ID"] = structure.name
    counts["heavy_atoms"] = len(structure.heavy_atoms)
    timings.append(("prepare_structure", seconds))
    print(f"    {'prepare_structure':<32}{seconds:9.2f}s", flush=True)

    surface_points, seconds = timed(construct_surface_grid, structure, R_PROBE)
    counts["surface_points"] = len(surface_points)
    timings.append(("construct_surface_grid (SASA)", seconds))
    print(f"    {'construct_surface_grid (SASA)':<32}{seconds:9.2f}s", flush=True)

    phases = [
        ("structure_features", calculate_structure_features, (structure, PH, R_PROBE, features)),
        ("surface_grid_features", calculate_surface_grid_features, (structure, surface_points, PH, HYDRO_SCALE, features)),
        ("average_charge_grid_features", calculate_average_chargesurface_grid_features, (structure, surface_points, PH, features)),
        ("shell_features", calculate_shell_features, (structure, surface_points, PH, features)),
    ]

    for phase, function, arguments in phases:
        _, seconds = timed(function, *arguments)
        timings.append((phase, seconds))
        print(f"    {phase:<32}{seconds:9.2f}s", flush=True)

    total = sum(seconds for _, seconds in timings)
    timings.append(("TOTAL", total))
    print(f"    {'TOTAL':<32}{total:9.2f}s", flush=True)

    return features, counts, timings


def write_features(pdb_file: Path, features: dict, workers: int):
    """Writes one run's feature values to their own CSV and returns the path.

    One row, one column per feature, written the way calculate() writes its
    output so that a reference file and a pipeline output can be compared
    without either having to be reformatted first.

    The worker count is in the filename rather than in a column, so that the
    fork's serial and parallel runs land in separate files and can be diffed
    against each other as well as against the original.

    Args:
        pdb_file: structure the features were calculated from.
        features: the feature dict returned by time_pipeline.
        workers: number of workers the run actually used.
    """

    path = FEATURE_DIR / f"features_{IMPLEMENTATION}_{pdb_file.stem}_{workers}w.csv"
    pd.Series(features).to_frame().transpose().to_csv(path, index=False)

    return path


def benchmark_structure(pdb_file: Path, writer, output):
    """Times one structure at every worker count and writes the rows to the CSV.

    Rows are written once a worker count has finished rather than after every
    phase, because the atom and surface point counts that go into each row are
    only known once the structure has been parsed. The feature values of each
    run go to their own file alongside, see write_features.

    Args:
        pdb_file: structure to benchmark.
        writer: csv.DictWriter already primed with the header.
        output: the underlying file handle, flushed after each worker count.
    """

    pdb_kb, atom_records, residues = structure_size(pdb_file)
    print(f"\n{pdb_file.name}  ({pdb_kb:.0f} KB, {residues} residues, {atom_records} ATOM records)")

    for workers in worker_counts():
        set_workers(workers)
        actual_workers = resolved_workers()
        print(f"  {IMPLEMENTATION}, {actual_workers} worker(s)", flush=True)

        features, counts, timings = time_pipeline(pdb_file)
        print(f"    features -> {write_features(pdb_file, features, actual_workers)}", flush=True)

        for phase, seconds in timings:
            writer.writerow(
                {
                    "timestamp": datetime.now(UTC).isoformat(timespec="seconds"),
                    "implementation": IMPLEMENTATION,
                    "prodes_path": PRODES_PATH,
                    "hostname": socket.gethostname(),
                    "platform": platform.platform(),
                    "logical_cpus": os.cpu_count(),
                    "mem_limit_mb": os.getenv("PRODES_MEM_LIMIT_MB", ""),
                    "full_features": os.getenv("PRODES_FULL_FEATURES", ""),
                    "structure": pdb_file.stem,
                    "pdb_kb": round(pdb_kb, 1),
                    "atom_records": atom_records,
                    "residues": residues,
                    "heavy_atoms": counts["heavy_atoms"],
                    "surface_points": counts["surface_points"],
                    "requested_workers": "default" if workers is None else workers,
                    "n_workers": actual_workers,
                    "phase": phase,
                    "seconds": round(seconds, 4),
                }
            )
        output.flush()


def main():
    """Benchmarks every structure in STRUCTURES and writes one CSV for this machine."""

    RESULTS_DIR.mkdir(exist_ok=True)
    FEATURE_DIR.mkdir(exist_ok=True)
    csv_path = RESULTS_DIR / f"benchmark_{IMPLEMENTATION}_{socket.gethostname()}.csv"

    print(f"implementation:  {IMPLEMENTATION}")
    print(f"prodes package:  {PRODES_PATH}")
    print(f"logical CPUs:    {os.cpu_count()}")
    print(f"worker counts:   {worker_counts()}")
    print(f"pka file:        {PKA_FILE}")
    print(f"structures:      {STRUCTURE_DIR}")
    print(f"writing:         {csv_path}")
    print(f"features:        {FEATURE_DIR}")

    structures = extract_structures(STRUCTURES, STRUCTURE_DIR)
    if not structures:
        raise FileNotFoundError(f"none of the benchmark structures exist, looked in {DATA_DIR}")

    with open(csv_path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        for pdb_file in structures:
            benchmark_structure(pdb_file, writer, output)

    print(f"\nwrote {csv_path}")


if __name__ == "__main__":
    main()
