"""Unit tests for prodes.calculations.geometry.

The vectorised functions in geometry are the hot path of the whole pipeline, and
until now they were only covered indirectly through the ARH96693 regression
fixture. The tests here pin them down directly.

Where a function was vectorised, it is checked against the original loop based
implementation from https://github.com/Tneijenhuis/prodes, reproduced below.
Where a batch function was added alongside a scalar one that still exists, the
two are checked against each other instead, so no code has to be duplicated.
"""

import numpy as np
import pytest

from prodes.calculations import geometry
from prodes.calculations.grid_wizard import Grid
from prodes.calculations.sasa import shrake_rupley
from prodes.core.point import Point
from prodes.io.parser import PDBparser

PDB_PATH = "tests/data/ARH96693.pdb.zip"
PH = 7

# --- Run the expensive setup once at import time, as the other test modules do ---

structure = PDBparser().parse(PDB_PATH)
_grid = Grid(12)
_grid.construct_cells(structure.heavy_atoms)
_grid.fill_cells(structure.heavy_atoms)

surface_points = shrake_rupley(_grid)
surface_coords = np.array([[p.x, p.y, p.z] for p in surface_points])

charged_atoms = [atom for atom in structure.heavy_atoms if atom.charge(PH) != 0]
charged_coords = np.array([[a.x, a.y, a.z] for a in charged_atoms])
charged_charges = np.array([a.charge(PH) for a in charged_atoms])

shell_points = geometry.Sunflower_sphere(structure.x, structure.y, structure.z, 1, 8).points


# --- The original loop based implementations, used as the oracle ---


def reference_maximal_distance(normal_vector, vector_on_plane, points):
    """The original per-point loop that geometry.maximal_distance replaced.

    Copied from the upstream prodes geometry module so that the vectorised
    version has something independent to be checked against.
    """

    unit_vector = normal_vector / np.linalg.norm(normal_vector)
    maximum = 0
    for point in points:
        point_vector = geometry.make_vector(point)
        dot_prod = np.dot(unit_vector, point_vector - vector_on_plane)
        if dot_prod > maximum:
            maximum = dot_prod

    return maximum


def reference_required_distance(point_for_plane, structure, points):
    """The original required_distance, including the +1 margin."""

    vector_on_plane = geometry.make_vector(structure)
    normal_vector = geometry.make_vector(point_for_plane) - vector_on_plane

    return reference_maximal_distance(normal_vector, vector_on_plane, points) + 1


# --- make_vector, find_plane, move_point ---


def test_make_vector_returns_the_coordinates():
    """make_vector turns anything with x, y and z into a 3 element array."""
    assert list(geometry.make_vector(Point(1.5, -2.0, 3.25))) == [1.5, -2.0, 3.25]


def test_find_plane_passes_through_its_point():
    """The plane returned for a point satisfies ax + by + cz = d at that point."""
    point = Point(structure.x, structure.y, structure.z + 10)

    a, b, c, d = geometry.find_plane(point, structure)

    assert a * point.x + b * point.y + c * point.z == pytest.approx(d)


def test_find_plane_normal_points_from_the_plane_to_the_structure():
    """The plane normal is the vector from the plane point to the structure centroid."""
    point = Point(structure.x, structure.y, structure.z + 10)

    a, b, c, _ = geometry.find_plane(point, structure)

    assert [a, b, c] == pytest.approx(list(geometry.make_vector(structure) - geometry.make_vector(point)))


def test_move_point_sets_the_distance_from_the_origin():
    """move_point places the point at the requested distance from the origin."""
    point = Point(structure.x + 3, structure.y, structure.z)

    geometry.move_point(point, structure, 25)

    assert np.linalg.norm(geometry.make_vector(point) - geometry.make_vector(structure)) == pytest.approx(25)


def test_move_point_keeps_the_direction():
    """move_point changes only the magnitude, not the direction, of the point."""
    point = Point(structure.x + 3, structure.y + 4, structure.z + 12)
    before = geometry.make_vector(point) - geometry.make_vector(structure)

    geometry.move_point(point, structure, 25)
    after = geometry.make_vector(point) - geometry.make_vector(structure)

    assert after / np.linalg.norm(after) == pytest.approx(before / np.linalg.norm(before))


# --- Vectorised against the original loop ---


@pytest.mark.parametrize("index", range(len(shell_points)))
def test_maximal_distance_matches_the_original_loop(index):
    """The vectorised maximal_distance agrees with the original per-point loop."""
    point = shell_points[index]
    origin = geometry.make_vector(structure)
    normal = geometry.make_vector(point) - origin

    assert geometry.maximal_distance(normal, origin, surface_coords) == pytest.approx(reference_maximal_distance(normal, origin, surface_points), rel=1e-12)


@pytest.mark.parametrize("index", range(len(shell_points)))
def test_required_distance_matches_the_original_loop(index):
    """The vectorised required_distance agrees with the original, +1 margin included."""
    point = shell_points[index]

    assert geometry.required_distance(point, structure, surface_coords) == pytest.approx(
        reference_required_distance(point, structure, surface_points), rel=1e-12
    )


# --- Batch against the scalar implementation it was added beside ---


def test_project_point_batch_matches_the_scalar_version():
    """project_point_batch reproduces project_point for every atom."""
    a, b, c, d = geometry.find_plane(shell_points[0], structure)

    batch = geometry.project_point_batch(a, b, c, d, charged_coords)
    scalar = np.array([geometry.project_point(a, b, c, d, *coord) for coord in charged_coords])

    assert batch == pytest.approx(scalar)


def test_projected_points_lie_on_the_plane():
    """Every projected point satisfies the plane equation it was projected onto."""
    a, b, c, d = geometry.find_plane(shell_points[0], structure)

    projected = geometry.project_point_batch(a, b, c, d, charged_coords)

    assert projected @ np.array([a, b, c]) == pytest.approx(np.full(len(projected), d))


def test_map_ep_to_plane_batch_matches_the_scalar_version():
    """map_ep_to_plane_batch reproduces map_ep_to_plane for every charged atom."""
    a, b, c, d = geometry.find_plane(shell_points[0], structure)
    projected = geometry.project_point_batch(a, b, c, d, charged_coords)
    exits, has_exit = geometry.find_exit_batch(charged_coords, projected, surface_coords)

    batch = geometry.map_ep_to_plane_batch(charged_coords, projected, exits, has_exit, charged_charges)
    scalar = np.array([geometry.map_ep_to_plane(atom, projected[i], exits[i], PH) for i, atom in enumerate(charged_atoms)])

    assert batch[has_exit] == pytest.approx(scalar[has_exit])


def test_map_ep_to_plane_batch_is_zero_without_an_exit():
    """An atom whose ray never leaves the protein contributes no potential."""
    a, b, c, d = geometry.find_plane(shell_points[0], structure)
    projected = geometry.project_point_batch(a, b, c, d, charged_coords)
    exits, _ = geometry.find_exit_batch(charged_coords, projected, surface_coords)

    no_exit = np.zeros(len(charged_coords), dtype=bool)
    potentials = geometry.map_ep_to_plane_batch(charged_coords, projected, exits, no_exit, charged_charges)

    assert not np.any(potentials)


# --- find_exit_batch, including the memory limited chunking ---


@pytest.mark.parametrize("mem_limit_mb", [0.001, 0.05, 1, 2048])
def test_find_exit_batch_is_unaffected_by_the_memory_limit(mem_limit_mb):
    """Chunking to stay within a memory budget must not change the result.

    The smallest limits here force many chunks, including one atom per chunk,
    which is what the low RAM configuration does on a large protein.
    """
    a, b, c, d = geometry.find_plane(shell_points[0], structure)
    projected = geometry.project_point_batch(a, b, c, d, charged_coords)

    unchunked, unchunked_has_exit = geometry.find_exit_batch(charged_coords, projected, surface_coords, mem_limit_mb=1e6)
    chunked, chunked_has_exit = geometry.find_exit_batch(charged_coords, projected, surface_coords, mem_limit_mb=mem_limit_mb)

    assert np.array_equal(chunked_has_exit, unchunked_has_exit)
    assert chunked[chunked_has_exit] == pytest.approx(unchunked[unchunked_has_exit])


def test_find_exit_batch_handles_no_atoms():
    """A structure with no charged atoms gives empty results rather than an error."""
    exits, has_exit = geometry.find_exit_batch(np.empty((0, 3)), np.empty((0, 3)), surface_coords)

    assert exits.shape == (0, 3)
    assert has_exit.shape == (0,)


# --- Sunflower_sphere ---


def test_raw_coordinates_match_the_point_objects():
    """raw_coordinates is the array form of the same points that .points builds."""
    sphere = geometry.Sunflower_sphere(1.0, 2.0, 3.0, 4.0, 50)

    from_points = np.array([[p.x, p.y, p.z] for p in sphere.points])

    assert sphere.raw_coordinates == pytest.approx(from_points)


def test_sphere_points_lie_on_the_sphere():
    """Every generated point sits at the requested radius from the centre."""
    centre, radius = np.array([1.0, 2.0, 3.0]), 4.0
    sphere = geometry.Sunflower_sphere(*centre, radius, 200)

    distances = np.linalg.norm(sphere.raw_coordinates - centre, axis=1)

    assert distances == pytest.approx(np.full(len(distances), radius))


def test_sphere_produces_the_requested_number_of_points():
    """The sphere holds exactly as many points as were asked for."""
    assert len(geometry.Sunflower_sphere(0, 0, 0, 1, 37).raw_coordinates) == 37
