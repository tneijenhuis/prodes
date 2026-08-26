"""Full-pipeline regression tests for the Prodes calculate() function.

Also holds the unit tests for the SASA task building at the bottom of the file,
which is the other half of what prodes.calculations.sasa does.

How this file works:
  1. Run calculate() on a small test PDB (ARH96693) twice — once with the full
     105-feature set and once with the reduced default set.
  2. Load a saved reference CSV (ARH96693_prodes_orig_output.csv) that was
     produced by the original, pre-refactoring code.
  3. Compare the freshly calculated output against the reference to catch
     regressions in any calculation (SASA, electrostatic potential,
     lipophilicity, shell features, etc.).

The expensive calculate() calls happen once at module load time (not per test).
Each test function then just reads from the pre-computed DataFrames below.
"""

from math import pi

import numpy as np
import pandas as pd

from prodes.calculations import sasa
from prodes.calculations.geometry import Sunflower_sphere
from prodes.calculations.grid_wizard import Grid
from prodes.io import parser as ps
from tests.pipeline_output import pipeline_output

PDB_PATH = "tests/data/ARH96693.pdb.zip"
ORIG_OUTPUT = "tests/data/ARH96693_prodes_orig_output.csv"

# Columns that hold integer counts (not floats) — given a looser tolerance.
COUNT_COLUMNS = {
    "NSurfPoints",
    "NSurfPosEpFormal",
    "NSurfNegEpFormal",
    "NSurfPosMhp",
    "NSurfNegMhp",
    "NSurfPosEpAverage",
    "NSurfNegEpAverage",
    "NShellPosEpFormal",
    "NShellNegEpFormal",
}

# --- Run the pipeline once at import time ---
# These are plain variables, not pytest fixtures. calculate() is expensive
# (SASA + 120 shell planes + electrostatic/lipophilic projections), so it runs
# once per feature set for the whole session, shared with the other module that
# needs the same two runs.

# Unscreened, because this module's job is to prove the fork still reproduces the
# original implementation, and the reference file was generated before screening
# existed. Screening is exercised by tests/test_screening.py.
FULL_OUTPUT = pipeline_output(PDB_PATH, full_features=True, ionic_strength_molar=0)
REDUCED_OUTPUT = pipeline_output(PDB_PATH, full_features=False, ionic_strength_molar=0)

ORIGINAL_OUTPUT = pd.read_csv(ORIG_OUTPUT)


def test_all_columns_present():
    """Every column in the original output should appear in the calculated output."""
    calc_cols = set(FULL_OUTPUT.columns)
    orig_cols = set(ORIGINAL_OUTPUT.columns)
    assert calc_cols == orig_cols, f"Column mismatch: missing={orig_cols - calc_cols}, extra={calc_cols - orig_cols}"


def test_feature_values_match_original():
    """All feature values should match the original output within tolerance."""
    calc = FULL_OUTPUT.iloc[0].drop("ID")
    orig = ORIGINAL_OUTPUT.iloc[0].drop("ID")

    mismatches = []
    for col in orig.index:
        calc_val = float(calc[col])
        orig_val = float(orig[col])

        if col in COUNT_COLUMNS:
            tolerance = 5
        elif "Shell" in col:
            tolerance = 0.05
        else:
            tolerance = 0.01

        if abs(calc_val - orig_val) > tolerance:
            mismatches.append(f"{col}: orig={orig_val}, calc={calc_val}, " f"diff={abs(calc_val - orig_val):.6f}, tol={tolerance}")

    assert len(mismatches) == 0, f"Value mismatches for {len(mismatches)} columns:\n" + "\n".join(mismatches[:20])


def test_reduced_output_is_a_strict_column_subset_of_full():
    """The reduced run must drop columns and invent none.

    Which columns specifically is not asserted here. That split lives in the YAML
    dictionaries and is checked against a real reduced run by
    tests/test_feature_dictionary.py, so repeating the names here would be a third
    copy of the same list.
    """
    reduced, full = set(REDUCED_OUTPUT.columns), set(FULL_OUTPUT.columns)

    assert reduced < full, f"reduced run has columns the full run does not: {reduced - full}"


# --- SASA task building ---
#
# The neighbour arrays are built once per grid cell and shared by every atom in
# it. The oracle below is the earlier per-atom version, kept here so the sharing
# is checked against an implementation that cannot have the same mistake in it.


def atom_grid(pdb_path=PDB_PATH, probe_r=1.4):
    """Returns a grid filled with the heavy atoms of a structure, as the pipeline builds it."""

    structure = ps.PDBparser().parse(pdb_path)
    grid = Grid(10)
    grid.construct_cells(structure.heavy_atoms)
    grid.fill_cells(structure.heavy_atoms)

    return structure, grid


def exposed_per_atom_oracle(grid, probe_r=1.4):
    """Shrake-Rupley with one neighbour array pair per atom, keyed by atom identity."""

    results = {}
    for cell in grid.cells.flatten():
        if cell.empty:
            continue

        neighbourhood = grid.grid_content("Atom", cells=grid.find_surrounding_cells(cell))
        for atom in cell.filtered_content("Atom"):
            if atom.radius is None:
                continue

            valid = [n for n in neighbourhood if n.radius is not None and n is not atom]
            atom_coord = np.array([atom.x, atom.y, atom.z])
            total_radius = atom.radius + probe_r

            sphere_coords = Sunflower_sphere(atom.x, atom.y, atom.z, total_radius, int(total_radius**2 * 4 * pi) * 2).raw_coordinates
            mask = np.ones(sphere_coords.shape[0], dtype=bool)

            if valid:
                neighbour_coords = np.array([[n.x, n.y, n.z] for n in valid])
                neighbour_radii = np.array([n.radius + probe_r for n in valid])
                close = np.sqrt(np.sum((neighbour_coords - atom_coord) ** 2, axis=1)) < (neighbour_radii + total_radius)

                if np.any(close):
                    diff = sphere_coords[:, np.newaxis, :] - neighbour_coords[close][np.newaxis, :, :]
                    distances = np.sqrt(np.sum(diff**2, axis=2))
                    mask = np.all(distances >= neighbour_radii[close][np.newaxis, :], axis=1)

            results[id(atom)] = (mask.sum() / len(mask), sphere_coords[np.where(mask)[0]])

    return results


def test_shared_neighbourhoods_match_a_per_atom_build():
    """Sharing one neighbour array pair per cell gives every atom the same answer as its own copy."""

    structure, grid = atom_grid()
    expected = exposed_per_atom_oracle(grid)

    sasa.shrake_rupley(grid)

    for atom in structure.heavy_atoms:
        if atom.radius is None:
            continue

        want_exposed, want_cloud = expected[id(atom)]
        assert atom.exposed == want_exposed, f"{atom.name}{atom.residue_number}: exposed {atom.exposed} != {want_exposed}"

        cloud = np.array([[p.x, p.y, p.z] for p in atom.cloud]) if len(atom.cloud) else np.empty((0, 3))
        assert np.array_equal(cloud, want_cloud), f"{atom.name}{atom.residue_number}: surface points differ"


def test_an_atom_does_not_occlude_itself():
    """A lone atom is fully exposed, so its own sphere is not counted against it.

    Sharing the neighbour arrays means an atom now finds itself in its own
    neighbourhood and has to mask its own row out. Missing that would not fail
    outright: a sphere point sits at exactly total_radius from its own centre and
    the surface test is >=, so most points survive and only the few that land a
    few ULP short are lost.
    """

    _structure, grid = atom_grid()
    _atoms, neighbourhoods, tasks = sasa.sasa_task_data(grid)
    neighbourhood_index, self_index, atom_coord, total_radius = tasks[0]
    coords, radii = neighbourhoods[neighbourhood_index]

    # A neighbourhood holding nothing but the atom itself, so anything it can
    # occlude, it occludes with its own sphere.
    sasa.SASA_TASKS["neighbourhoods"] = [(coords[[self_index]], radii[[self_index]])]
    sasa.SASA_TASKS["tasks"] = [(0, 0, atom_coord, total_radius)]
    try:
        exposed, _points = sasa.process_atom(0)
    finally:
        sasa.SASA_TASKS.clear()

    assert exposed == 1.0


def test_every_atom_knows_its_own_row_in_the_shared_arrays():
    """Each task points at the atom's own entry in its cell's neighbour arrays."""

    _structure, grid = atom_grid()
    atoms, neighbourhoods, tasks = sasa.sasa_task_data(grid)

    for atom, (neighbourhood_index, self_index, atom_coord, _total_radius) in zip(atoms, tasks, strict=True):
        neighbour_coords, _radii = neighbourhoods[neighbourhood_index]

        assert self_index >= 0, f"{atom.name}{atom.residue_number} is missing from its own neighbourhood"
        assert np.array_equal(neighbour_coords[self_index], atom_coord)


def test_neighbour_arrays_are_built_once_per_cell():
    """There is one neighbour array pair per occupied cell, not one per atom."""

    _structure, grid = atom_grid()
    atoms, neighbourhoods, _tasks = sasa.sasa_task_data(grid)
    occupied = sum(1 for cell in grid.cells.flatten() if not cell.empty)

    assert len(neighbourhoods) == occupied
    assert len(neighbourhoods) < len(atoms)
