"""Turns the raw benchmark CSVs into the anonymised dataset that ships with the repo.

Raw runs land in benchmark_results/, which is git ignored, and carry the file
system path of the package that was measured and the real hostname of the
machine it ran on. Neither belongs in a public repository, so this script strips
them, renames the implementations to their public names and writes the result as
a zip under docs/benchmark/::

    python scripts/publish_benchmark.py

That zip is what plot_benchmark.py reads and what anyone re-plotting the data
unpacks, so every committed figure can be reproduced from committed data alone.

Only the CSVs directly in benchmark_results/ are published. Anything moved into
benchmark_results/archive/ is ignored, which is how a run measured on a machine
that should not appear in the published figures is left out without deleting it.
"""

import zipfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "benchmark_results"
DATA_DIR = REPO_ROOT / "tests" / "data"
PUBLISHED_DIR = REPO_ROOT / "docs" / "benchmark"
DATA_ARCHIVE = PUBLISHED_DIR / "benchmark_data.zip"

# Written into every row by benchmark.py and useful only on the machine that
# produced it. It is an absolute path into someone's home directory or a mounted
# share, so it is dropped rather than anonymised.
INTERNAL_COLUMNS = ["prodes_path"]

# The implementations have each been recorded under two names: the earliest runs
# wrote "fork" and "original", and benchmark.py now takes the name from the
# installed package. Both spellings map onto the public name of the repository.
PUBLIC_IMPLEMENTATIONS = {
    "fork": "datacatalysis/prodes",
    "datacatalysis-prodes": "datacatalysis/prodes",
    "original": "tneijenhuis/prodes",
    "tneijenhuis-prodes": "tneijenhuis/prodes",
}


def load_raw():
    """Reads every benchmark CSV directly in the results directory.

    Not recursive, so benchmark_results/archive/ is skipped.

    Raises:
        FileNotFoundError: if no benchmark CSV has been written yet.
    """

    csv_paths = sorted(RESULTS_DIR.glob("benchmark_*.csv"))
    if not csv_paths:
        raise FileNotFoundError(f"no benchmark CSVs found in {RESULTS_DIR}, run scripts/benchmark.py first")

    print(f"read {len(csv_paths)} CSV file(s) from {RESULTS_DIR}")

    return pd.concat([pd.read_csv(path) for path in csv_paths], ignore_index=True)


def anonymise(raw):
    """Returns the raw results with the internal details removed.

    Hostnames are replaced by server1, server2 and so on, numbered in the order
    they appear. The mapping is deliberately not written down anywhere: this
    script is committed, so a lookup table in it would publish exactly the names
    it is meant to remove.

    Raises:
        ValueError: if a CSV names an implementation with no public name.
    """

    unknown = sorted(set(raw["implementation"]) - set(PUBLIC_IMPLEMENTATIONS))
    if unknown:
        raise ValueError(f"no public name for implementation(s) {unknown}, add them to PUBLIC_IMPLEMENTATIONS")

    published = raw.drop(columns=INTERNAL_COLUMNS)
    published["implementation"] = published["implementation"].map(PUBLIC_IMPLEMENTATIONS)

    hostnames = {hostname: f"server{number}" for number, hostname in enumerate(published["hostname"].unique(), start=1)}
    published["hostname"] = published["hostname"].map(hostnames)
    print(f"anonymised {len(hostnames)} machine(s) to {', '.join(sorted(hostnames.values()))}")

    # Earlier runs recorded the file name ("1GDW.pdb") where later ones record
    # the stem, and the same structure has to carry the same name in both.
    published["structure"] = published["structure"].str.removesuffix(".pdb")

    return published


def residue_counts():
    """Returns the residue count of every benchmark structure in tests/data.

    Counted as distinct (chain, residue number) pairs, matching how
    benchmark.py counts them, so that a backfilled row and a freshly measured
    one agree.

    Returns:
        dict mapping structure name to the number of residues over all chains.
    """

    counts = {}

    for archive in sorted(DATA_DIR.glob("*.pdb.zip")):
        with zipfile.ZipFile(archive) as zipped:
            name = next(member for member in zipped.namelist() if member.endswith(".pdb"))
            lines = zipped.read(name).decode().splitlines()

        residues = {(line[21], line[22:27]) for line in lines if line.startswith("ATOM")}
        counts[Path(name).stem] = len(residues)

    return counts


def backfill_residues(published):
    """Adds the residue count to results measured before benchmark.py recorded it.

    The figures are plotted against sequence length, so the column has to be
    there for every row. Runs made after benchmark.py started writing it are
    left alone; older ones are filled in from the structures in tests/data,
    which are the same files the benchmark unpacks.

    Raises:
        KeyError: if a structure in the results has no matching file, since
            guessing a sequence length would silently corrupt the x axis.
    """

    if "residues" in published.columns and published["residues"].notna().all():
        return published

    counts = residue_counts()
    missing = sorted(set(published["structure"]) - set(counts))
    if missing:
        raise KeyError(f"no structure file in {DATA_DIR} for {missing}, cannot determine their sequence length")

    published["residues"] = published["structure"].map(counts)
    print(f"backfilled residue counts for {len(counts)} structure(s) from {DATA_DIR.name}")

    # Put it where benchmark.py writes it, so a backfilled CSV and a freshly
    # measured one have the same column order.
    columns = [column for column in published.columns if column != "residues"]
    position = columns.index("atom_records") + 1

    return published[columns[:position] + ["residues"] + columns[position:]]


def archive_name(implementation, hostname):
    """Returns the file name to store one implementation and machine under.

    The public implementation names contain a slash, which cannot appear in a
    file name inside the archive.
    """

    return f"benchmark_{implementation.replace('/', '-')}_{hostname}.csv"


def write_archive(published, output_path):
    """Writes one CSV per implementation and machine into a single zip.

    Splitting them keeps the published layout the same as the raw one, where
    each run writes its own file, so an unpacked archive can be dropped straight
    back into benchmark_results/ and re-plotted.
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for (implementation, hostname), group in published.groupby(["implementation", "hostname"]):
            name = archive_name(implementation, hostname)
            archive.writestr(name, group.to_csv(index=False))
            print(f"  {name}  ({len(group)} rows)")

    print(f"wrote {output_path}")


def main():
    """Anonymises the raw benchmark CSVs and writes the published archive."""

    published = backfill_residues(anonymise(load_raw()))
    write_archive(published, DATA_ARCHIVE)


if __name__ == "__main__":
    main()
