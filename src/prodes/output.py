"""Assembles the output of a Prodes run into one self contained bundle.

Prodes used to append a row to a growing CSV. That made a run depend on whatever
happened to be at the output path already, which is unusable under
multiprocessing, and it threw away everything except the summary numbers. A run
now takes one structure and produces one bundle, which is what the rest of a
bioinformatics workflow expects.

The bundle holds the features, the surface points they were calculated from, and
viewer files that open those points directly. Because all of it comes out of the
same run, a figure and a feature value cannot disagree.
"""

import json
import zipfile
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from prodes.io.parser import read_pdb_text
from prodes.viz import (
    HYDROPHOBICITY_DESCRIPTION,
    charged_residue_underlay,
    chimerax_ep_colouring,
    chimerax_hydrophobicity_colouring,
    chimerax_script,
    ep_colouring,
    hydrophobic_residue_underlay,
    hydrophobicity_colouring,
    pymol_script,
    symmetric_limit,
    write_point_pdb,
)

README = """Prodes output bundle for {name}
{underline}

This bundle holds the calculated features, the surface points they were
calculated from, and two ready-made views of that surface.


LOOKING AT THE SURFACE
----------------------

The .pml files are PyMOL scripts: plain text files of PyMOL commands. Opening
one loads the structure and colours its surface in a single step, so there is
nothing to set up by hand.

Open a terminal in this directory and run:

    pymol {name}_ep.pml               (electrostatic potential)
    pymol {name}_hydrophobicity.pml   (hydrophobicity)

The paths inside the scripts are relative, so they only work when opened from
inside this directory. If PyMOL is already running, use File, Run Script instead
and pick the .pml file.

  Electrostatic potential   red is negative, blue positive, white near zero.
                            The scale is set from this protein's own range, and
                            is written in the script if you want to match it
                            across two proteins.

  Hydrophobicity            grey is hydrophilic, pale green hydrophobic, dark
                            green strongly so. These cutoffs are fixed, so the
                            same green means the same thing on any protein.

The .cxc files are the same two views for ChimeraX, if you prefer it to PyMOL.


WHAT IS IN HERE
---------------

  {name}_features.csv              the calculated features, one row
  {name}_surface_points.csv        every surface point: x, y, z, potential,
                                   hydrophobicity
  {name}_ep.pml                    PyMOL, coloured by electrostatic potential
  {name}_hydrophobicity.pml        PyMOL, coloured by hydrophobicity
  {name}_ep.cxc                    ChimeraX, electrostatic potential
  {name}_hydrophobicity.cxc        ChimeraX, hydrophobicity
  {name}_ep.pdb                    the surface points, potential in the
                                   B factor column
  {name}_hydrophobicity.pdb        the surface points, hydrophobicity in the
                                   B factor column
  {name}.pdb                       the structure this run was given
  prodes_run.json                  version, settings and time of the run

The point clouds and the features come out of the same calculation, so a picture
here can never disagree with a number in the feature table.


WHAT THE POTENTIAL IS, AND IS NOT
---------------------------------

The electrostatic potential is a Coulomb sum at a uniform relative permittivity
of 4, damped with distance to account for the mobile ions in the buffer. The
ionic strength used is recorded in prodes_run.json; 0 means no damping.

It has no dielectric boundary, and the damping length is the one for bulk water
rather than the value consistent with a permittivity of 4, so although it is
reported in volts it is not comparable to a Poisson-Boltzmann calculation such
as APBS. Treat it as a relative description of one protein's surface.

Note the ShellEp features in the feature table are computed by a different route
and are not damped. That route already weights the solvent part of the path far
more heavily, so it did not have the problem the damping fixes. See the
repository README.
"""


def prodes_version():
    """Returns the installed package version, or 'unknown' outside an install."""

    try:
        return version("prodes")
    except PackageNotFoundError:
        return "unknown"


def run_metadata(pdb_file, settings, n_points, ep):
    """Returns the record of what was run, written into every bundle."""

    return {
        "prodes_version": prodes_version(),
        "utc": datetime.now(UTC).isoformat(),
        "input_file": str(pdb_file),
        "settings": settings,
        "surface_points": int(n_points),
        "ep_min_volts": round(float(ep.min()), 3) if n_points else None,
        "ep_max_volts": round(float(ep.max()), 3) if n_points else None,
    }


def bundle_files(directory, name, features, coords, ep, lipo, pdb_file, metadata, hydro_scale):
    """Writes every member of the bundle into an existing directory."""

    structure_file, point_file = f"{name}.pdb", f"{name}_ep.pdb"

    features.to_csv(directory / f"{name}_features.csv", index=False)
    pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1], "z": coords[:, 2], "ep_volts": ep, "hydrophobicity": lipo}).to_csv(
        directory / f"{name}_surface_points.csv", index=False
    )

    # Two parallel views of the same points: the electrostatic potential and the
    # projected hydrophobicity. Both are already in the point table; without a
    # point cloud carrying each in its B factor column there is no way to look at
    # them in a viewer.
    ep_limit = symmetric_limit(ep)
    write_point_pdb(coords, ep, directory / point_file)
    (directory / f"{name}_ep.pml").write_text(
        pymol_script(structure_file, point_file, ep_colouring("surface_ep", ep_limit), charged_residue_underlay("surface_ep"))
    )
    (directory / f"{name}_ep.cxc").write_text(chimerax_script(structure_file, point_file, chimerax_ep_colouring(ep_limit)))

    phobic_file = f"{name}_hydrophobicity.pdb"
    write_point_pdb(coords, lipo, directory / phobic_file)
    (directory / f"{name}_hydrophobicity.pml").write_text(
        pymol_script(
            structure_file,
            phobic_file,
            hydrophobicity_colouring("hydrophobicity"),
            hydrophobic_residue_underlay("hydrophobicity", hydro_scale),
            name="hydrophobicity",
            description=HYDROPHOBICITY_DESCRIPTION,
        )
    )
    (directory / f"{name}_hydrophobicity.cxc").write_text(chimerax_script(structure_file, phobic_file, chimerax_hydrophobicity_colouring()))

    # The structure travels with the bundle, or the viewer scripts have nothing
    # to load once it has been copied to another machine. read_pdb_text reads
    # straight through a .pdb.zip, so a zipped input needs no special case.
    (directory / structure_file).write_text(read_pdb_text(pdb_file))

    (directory / "prodes_run.json").write_text(json.dumps(metadata, indent=2))
    (directory / "README.txt").write_text(README.format(name=name, underline="=" * (len(name) + 26)))


def read_member(bundle, suffix):
    """Reads the one bundle member whose name ends in suffix, as a dataframe."""

    with zipfile.ZipFile(bundle) as archive:
        try:
            name = next(member for member in archive.namelist() if member.endswith(suffix))
        except StopIteration:
            raise ValueError(f"{bundle} holds no member ending in {suffix}") from None

        with archive.open(name) as handle:
            return pd.read_csv(handle)


def read_features(bundle):
    """Reads the feature table out of a bundle.

    The counterpart to write_bundle, for a pipeline that wants the numbers
    rather than the viewer files.
    """

    return read_member(bundle, "_features.csv")


def read_surface_points(bundle):
    """Reads the surface point table out of a bundle: coordinates, potential, hydrophobicity."""

    return read_member(bundle, "_surface_points.csv")


def write_bundle(out_file, name, features, coords, ep, lipo, pdb_file, metadata, hydro_scale):
    """Writes the bundle to out_file, which must be a .zip path.

    One structure in, one archive out. Supporting a directory as well would mean
    every caller and every test branching on which shape it got, for something a
    caller can do in one line with unzip.

    Raises:
        ValueError: if out_file does not end in .zip.
    """

    out_file = Path(out_file)
    if out_file.suffix.lower() != ".zip":
        raise ValueError(f"the output must be a .zip path, got {out_file}")

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as scratch:
        staged = Path(scratch) / name
        staged.mkdir()
        bundle_files(staged, name, features, coords, ep, lipo, pdb_file, metadata, hydro_scale)

        with zipfile.ZipFile(out_file, "w", zipfile.ZIP_DEFLATED) as archive:
            for member in sorted(staged.iterdir()):
                archive.write(member, f"{name}/{member.name}")

    return out_file
