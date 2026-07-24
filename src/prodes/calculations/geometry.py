import os

import numpy as np

from prodes.calculations.distance_functions import atom_charge_coulomb, potential_multiple_media


# inspired by https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
class Sunflower_sphere:

    def __init__(self, x, y, z, r, n_points):

        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.n_points = n_points
        self._points = None
        self._raw_coords = None

    @property
    def points(self):

        if self._points is None:

            from prodes.core.point import Point
            indices = np.arange(0, self.n_points, dtype=float) + 0.5

            phi = np.arccos(1 - 2*indices/self.n_points)
            theta = np.pi * (1 + 5**0.5) * indices

            x, y, z = (
                np.cos(theta) * np.sin(phi) * self.r + self.x,
                np.sin(theta) * np.sin(phi) * self.r + self.y,
                np.cos(phi) * self.r + self.z
            )

            self.angles = np.array([theta, phi])

            points = np.array([Point(xi, yi, zi) for xi, yi, zi in zip(x, y, z)])

            self._points = points

        return self._points

    @property
    def raw_coordinates(self):
        """Returns sphere point coordinates as an (N, 3) numpy array without wrapping in Point objects."""

        if self._raw_coords is None:
            indices = np.arange(0, self.n_points, dtype=float) + 0.5
            phi = np.arccos(1 - 2*indices/self.n_points)
            theta = np.pi * (1 + 5**0.5) * indices
            x = np.cos(theta) * np.sin(phi) * self.r + self.x
            y = np.sin(theta) * np.sin(phi) * self.r + self.y
            z = np.cos(phi) * self.r + self.z
            self._raw_coords = np.column_stack([x, y, z])

        return self._raw_coords


def make_vector(point):
    return np.array([point.x, point.y, point.z])


def find_plane(point_on_plane, normal_vector_point):
    """finds the linear equation of a plane"""

    point_on_plane_vector = make_vector(point_on_plane)

    a, b, c = make_vector(normal_vector_point) - point_on_plane_vector
    x_0, y_0, z_0 = point_on_plane_vector

    product = a * x_0 + b * y_0 + c * z_0

    return a, b, c, product


def move_point(point, origin, magnitude):
    """changes the magnitude of a point from a specific origin"""

    point_vector, origin_vector = make_vector(point), make_vector(origin)

    vector = point_vector - origin_vector
    unit_vector = vector/np.linalg.norm(vector)
    new_vector = unit_vector * magnitude + origin_vector
    point.x, point.y, point.z = new_vector


def maximal_distance(normal_vector, vector_on_plane, surface_coords):
    """Determines the maximal distance by calculating the dot product
    of each surface point and the unit normal vector.

    Args:
        normal_vector: (3,) array, direction from origin to shell point.
        vector_on_plane: (3,) array, origin (structure centroid).
        surface_coords: (N, 3) array of surface point coordinates.
    """
    unit_vector = normal_vector / np.linalg.norm(normal_vector)
    dot_prods = (surface_coords - vector_on_plane) @ unit_vector
    return float(dot_prods.max())


def required_distance(point_for_plane, structure, surface_coords):
    """Determines the required distance for the plane to be formed on the protein surface.

    Args:
        point_for_plane: Point object, a shell point.
        structure: Structure object with x/y/z centroid.
        surface_coords: (N, 3) array of surface point coordinates.
    """
    vector_on_plane = make_vector(structure)
    normal_vector = make_vector(point_for_plane) - vector_on_plane

    return maximal_distance(normal_vector, vector_on_plane, surface_coords) + 1


def project_point(a, b, c, d, x1, y1, z1):
    """projects point 1 onto plane ax+by+cz=-d"""

    k = (d - a * x1-b * y1-c * z1)/(a * a + b * b + c * c)
    x2 = a * k + x1
    y2 = b * k + y1
    z2 = c * k + z1
    return x2, y2, z2


def find_exit(point_vector, projected_point_vector, grid):
    """finds the position where the vector between two points leaves the protein"""
    from math import ceil

    from prodes.core.point import Point

    normal_vector = projected_point_vector - point_vector
    total_distance = np.linalg.norm(normal_vector)
    direction = normal_vector/total_distance
    highest = 0
    surface_exit = None

    cells = []
    for i in range(ceil(total_distance)*2):
        sample_point = point_vector+direction*i/2
        cell = grid.in_which_cell(Point(*sample_point))
        if cell not in cells:
            cells.append(cell)

    environment = []
    for cell in cells:
        x, y, z = cell
        try:
            environment.append(grid.cells[z, y, x])
        except IndexError:
            pass

    for cell in environment.copy():
        surrounding = grid.find_surrounding_cells(cell)
        for sur_cell in surrounding:
            if sur_cell not in cells:
                environment.append(sur_cell)

    surface_points = np.array(grid.grid_content(cells=environment))

    for surface_point in surface_points:
        surface_point_vector = make_vector(surface_point) - point_vector
        dot_prod = np.dot(direction, surface_point_vector)
        if dot_prod > 0:
            potential_exit = direction * dot_prod
            distance = round(np.linalg.norm(potential_exit-surface_point_vector), 1)

            if distance <= 1:
                if dot_prod > highest:
                    highest = dot_prod
                    surface_exit = potential_exit+point_vector

    return surface_exit


def _find_exit_chunk(atom_coords, projected_coords, surface_coords):
    """Memory-efficient find_exit for a chunk of atoms.

    Uses algebraic identities to avoid (M, N, 3) intermediate arrays:
      dot_prods[m,n]  = dir[m] · (surf[n] - atom[m])
                      = dir[m]·surf[n] - dir[m]·atom[m]
      perp_dist²[m,n] = ||surf[n]-atom[m]||² - dot_prods[m,n]²

    Only (M, N) float64 arrays are created instead of (M, N, 3).
    """
    normal_vectors = projected_coords - atom_coords
    total_distances = np.linalg.norm(normal_vectors, axis=1)
    directions = normal_vectors / total_distances[:, None]

    # dot_prods: (M, N) — avoids (M, N, 3) surface_vectors array
    dot_prods = directions @ surface_coords.T  # (M, N)
    dot_prods -= np.sum(atom_coords * directions, axis=1)[:, None]

    # surface_dist_sq: (M, N) — ||surf[n] - atom[m]||² via expansion
    surface_sq = np.sum(surface_coords ** 2, axis=1)  # (N,)
    atom_sq = np.sum(atom_coords ** 2, axis=1)        # (M,)
    surface_dist_sq = atom_sq[:, None] + surface_sq[None, :] - 2.0 * (atom_coords @ surface_coords.T)

    # perp_dist via Pythagorean decomposition (no (M, N, 3) array needed)
    perp_dist_sq = surface_dist_sq - dot_prods ** 2
    np.maximum(perp_dist_sq, 0, out=perp_dist_sq)
    perp_dist = np.sqrt(perp_dist_sq)

    positive = dot_prods > 0
    close = np.round(perp_dist, 1) <= 1
    valid = positive & close

    has_exit = valid.any(axis=1)
    dot_prods_masked = np.where(valid, dot_prods, -np.inf)
    max_indices = np.argmax(dot_prods_masked, axis=1)

    # exits[m] = directions[m] * dot_prods[m, max_idx] + atom_coords[m]
    exits = directions * dot_prods[np.arange(len(atom_coords)), max_indices][:, None] + atom_coords

    return exits, has_exit


def find_exit_batch(atom_coords, projected_coords, surface_coords, mem_limit_mb=None):
    """Vectorized find_exit for all charged atoms at once, with memory-bounded chunking.

    For each atom, finds where the ray from atom to its projected point on the
    shell plane exits the protein surface.

    Processes atoms in chunks to keep peak memory within *mem_limit_mb*.
    The limit defaults to the ``PRODES_MEM_LIMIT_MB`` env var (2048 MB if unset).

    Args:
        atom_coords: (M, 3) array of atom positions.
        projected_coords: (M, 3) array of projected positions on the plane.
        surface_coords: (N, 3) array of all surface point coordinates.
        mem_limit_mb: Maximum memory in MB for intermediate arrays per chunk.
            If None, reads from PRODES_MEM_LIMIT_MB env var (default 2048).

    Returns:
        exits: (M, 3) array of exit points (undefined where has_exit is False).
        has_exit: (M,) boolean array indicating which atoms found an exit.
    """
    n_atoms = len(atom_coords)
    if n_atoms == 0:
        return np.empty((0, 3)), np.empty(0, dtype=bool)

    if mem_limit_mb is None:
        mem_limit_mb = float(os.getenv("PRODES_MEM_LIMIT_MB", "2048"))

    n_surface = len(surface_coords)
    # Per atom we need ~5 (N,) float64 arrays + 1 bool array ≈ N * 48 bytes
    bytes_per_atom = n_surface * 48
    max_atoms_per_chunk = max(1, int(mem_limit_mb * 1024 * 1024 / bytes_per_atom))

    all_exits = np.empty((n_atoms, 3))
    all_has_exit = np.empty(n_atoms, dtype=bool)

    for start in range(0, n_atoms, max_atoms_per_chunk):
        end = min(start + max_atoms_per_chunk, n_atoms)
        chunk_exits, chunk_has_exit = _find_exit_chunk(
            atom_coords[start:end], projected_coords[start:end], surface_coords
        )
        all_exits[start:end] = chunk_exits
        all_has_exit[start:end] = chunk_has_exit

    return all_exits, all_has_exit


def map_ep_to_plane(atom, projected_point_vector, surface_exit, ph=7):
    """Calculates the electrostatic potential of an atom projected onto a plane"""
    atom_vector = make_vector(atom)

    total_distance = np.linalg.norm(projected_point_vector - atom_vector)
    protein_distance = np.linalg.norm(surface_exit - atom_vector)

    distances = [(total_distance - protein_distance) * 10**-10, protein_distance * 10**-10]
    atom_charge = atom_charge_coulomb(atom.charge(ph))

    point_potential = potential_multiple_media(atom_charge, distances, [80, 4])

    return point_potential


def map_ep_to_plane_batch(atom_coords, projected_coords, surface_exits, has_exit, charges):
    """Vectorized map_ep_to_plane for all atoms at once.

    Args:
        atom_coords: (M, 3) array of atom positions.
        projected_coords: (M, 3) array of projected positions.
        surface_exits: (M, 3) array of exit points.
        has_exit: (M,) boolean mask for valid exits.
        charges: (M,) array of atom charges (in elementary charge units).

    Returns:
        (M,) array of potentials (0 where no exit found).
    """
    total_distances = np.linalg.norm(projected_coords - atom_coords, axis=1)
    protein_distances = np.linalg.norm(surface_exits - atom_coords, axis=1)

    dist_water = (total_distances - protein_distances) * 1e-10
    dist_protein = protein_distances * 1e-10

    absolute_permittivity = 8.854e-12
    perm_water = 80 * absolute_permittivity
    perm_protein = 4 * absolute_permittivity

    denominator = perm_water * dist_water + perm_protein * dist_protein
    coulomb_charges = charges * 1.6e-19

    potentials = coulomb_charges / (denominator * 4 * np.pi)
    potentials = np.where(has_exit, potentials, 0.0)

    return potentials


def project_point_batch(a, b, c, d, coords):
    """Projects an array of points onto plane ax+by+cz=-d.

    Args:
        a, b, c, d: plane equation coefficients.
        coords: (M, 3) array of point coordinates.

    Returns:
        (M, 3) array of projected coordinates.
    """
    x1, y1, z1 = coords[:, 0], coords[:, 1], coords[:, 2]
    k = (d - a * x1 - b * y1 - c * z1) / (a * a + b * b + c * c)
    return coords + k[:, None] * np.array([a, b, c])
