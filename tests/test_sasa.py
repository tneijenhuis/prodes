import time
import numpy as np
import pytest

from prodes.io.parser import PDBparser
from prodes.calculations import grid_wizard
from prodes.calculations.sasa import shrake_rupley, shrake_rupley_original


PDB_PATH = "tests/data/ARH96693.pdb"


def build_grid_and_run(pdb_path, sasa_func, probe_r=1.4):
    """Parse PDB, build grid, run SASA calculation, return structure and surface."""
    structure = PDBparser().parse(pdb_path)
    grid_size = 10 + (probe_r - 1.4) * 2
    grid = grid_wizard.Grid(grid_size)
    grid.construct_cells(structure.heavy_atoms)
    grid.fill_cells(structure.heavy_atoms)
    surface = sasa_func(grid, probe_r)
    return structure, surface


def test_vectorized_matches_original_point_count():
    """The vectorized and original shrake_rupley should produce the same number of surface points."""
    _, surface_orig = build_grid_and_run(PDB_PATH, shrake_rupley_original)
    _, surface_vec = build_grid_and_run(PDB_PATH, shrake_rupley)
    assert len(surface_orig) == len(surface_vec), (
        f"Surface point count mismatch: original={len(surface_orig)}, vectorized={len(surface_vec)}"
    )


def test_vectorized_matches_original_total_area():
    """Total surface area should match within a small tolerance."""
    structure_orig, _ = build_grid_and_run(PDB_PATH, shrake_rupley_original)
    structure_vec, _ = build_grid_and_run(PDB_PATH, shrake_rupley)
    area_orig = structure_orig.surface_area(1.4)
    area_vec = structure_vec.surface_area(1.4)
    assert abs(area_orig - area_vec) < 1.0, (
        f"Total area mismatch: original={area_orig:.3f}, vectorized={area_vec:.3f}, "
        f"diff={abs(area_orig - area_vec):.3f}"
    )


def test_vectorized_matches_original_per_atom_exposed():
    """Per-atom exposed fraction should match within 1% tolerance."""
    structure_orig, _ = build_grid_and_run(PDB_PATH, shrake_rupley_original)
    structure_vec, _ = build_grid_and_run(PDB_PATH, shrake_rupley)
    for i, (atom_o, atom_v) in enumerate(
        zip(structure_orig.heavy_atoms, structure_vec.heavy_atoms)
    ):
        if atom_o.radius is None:
            continue
        assert abs(atom_o.exposed - atom_v.exposed) < 0.01, (
            f"Exposed fraction mismatch at atom {i} ({atom_o.name} {atom_o.residue_name}): "
            f"original={atom_o.exposed:.6f}, vectorized={atom_v.exposed:.6f}"
        )


def test_vectorized_matches_original_per_atom_cloud_size():
    """Per-atom cloud size (number of surface points per atom) should match exactly or within 1."""
    structure_orig, _ = build_grid_and_run(PDB_PATH, shrake_rupley_original)
    structure_vec, _ = build_grid_and_run(PDB_PATH, shrake_rupley)
    mismatches = []
    for i, (atom_o, atom_v) in enumerate(
        zip(structure_orig.heavy_atoms, structure_vec.heavy_atoms)
    ):
        if atom_o.radius is None:
            continue
        diff = abs(len(atom_o.cloud) - len(atom_v.cloud))
        if diff > 1:
            mismatches.append(
                f"atom {i} ({atom_o.name} {atom_o.residue_name}): "
                f"original={len(atom_o.cloud)}, vectorized={len(atom_v.cloud)}"
            )
    assert len(mismatches) == 0, (
        f"Cloud size mismatches (>1 point difference) for {len(mismatches)} atoms:\n"
        + "\n".join(mismatches[:10])
    )


def test_vectorized_matches_original_surface_point_coords():
    """Surface point coordinates should match within 0.01 Angstrom."""
    _, surface_orig = build_grid_and_run(PDB_PATH, shrake_rupley_original)
    _, surface_vec = build_grid_and_run(PDB_PATH, shrake_rupley)
    coords_orig = np.array([[p.x, p.y, p.z] for p in surface_orig])
    coords_vec = np.array([[p.x, p.y, p.z] for p in surface_vec])
    if len(coords_orig) == 0:
        assert len(coords_vec) == 0
        return
    max_diff = np.max(np.abs(coords_orig - coords_vec))
    assert max_diff < 0.01, f"Max coordinate difference: {max_diff:.6f}"


def test_vectorized_is_faster_or_equal():
    """The vectorized version should not be slower than the original."""
    t0 = time.time()
    build_grid_and_run(PDB_PATH, shrake_rupley_original)
    t_orig = time.time() - t0

    t0 = time.time()
    build_grid_and_run(PDB_PATH, shrake_rupley)
    t_vec = time.time() - t0

    print(f"\nOriginal: {t_orig:.2f}s, Vectorized: {t_vec:.2f}s, speedup: {t_orig/t_vec:.1f}x")
    assert t_vec <= t_orig * 1.5, (
        f"Vectorized is slower: original={t_orig:.2f}s, vectorized={t_vec:.2f}s"
    )
