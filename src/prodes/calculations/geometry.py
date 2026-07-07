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

    k =(d -a * x1-b * y1-c * z1)/(a * a + b * b + c * c)
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
            distance = round(np.linalg.norm(potential_exit-surface_point_vector),1)

            if distance <= 1:
                if dot_prod > highest:
                    highest = dot_prod
                    surface_exit = potential_exit+point_vector

    return surface_exit


def find_exit_batch(atom_coords, projected_coords, surface_coords):
    """Vectorized find_exit for all charged atoms at once.

    For each atom, finds where the ray from atom to its projected point on the
    shell plane exits the protein surface.

    Args:
        atom_coords: (M, 3) array of atom positions.
        projected_coords: (M, 3) array of projected positions on the plane.
        surface_coords: (N, 3) array of all surface point coordinates.

    Returns:
        exits: (M, 3) array of exit points (undefined where has_exit is False).
        has_exit: (M,) boolean array indicating which atoms found an exit.
    """
    normal_vectors = projected_coords - atom_coords
    total_distances = np.linalg.norm(normal_vectors, axis=1)
    directions = normal_vectors / total_distances[:, None]

    surface_vectors = surface_coords[None, :, :] - atom_coords[:, None, :]
    dot_prods = np.einsum('mnd,md->mn', surface_vectors, directions)

    positive = dot_prods > 0
    potential_exits = directions[:, None, :] * dot_prods[:, :, None]
    perp = surface_vectors - potential_exits
    perp_dist = np.linalg.norm(perp, axis=2)
    close = np.round(perp_dist, 1) <= 1
    valid = positive & close

    has_exit = valid.any(axis=1)
    dot_prods_masked = np.where(valid, dot_prods, -np.inf)
    max_indices = np.argmax(dot_prods_masked, axis=1)

    exits = potential_exits[np.arange(len(atom_coords)), max_indices] + atom_coords

    return exits, has_exit

def map_ep_to_plane(atom, projected_point_vector, surface_exit, ph=7):
    """Calculates the electrostatic potential of an atom projected onto a plane"""
    atom_vector = make_vector(atom)

    total_distance = np.linalg.norm(projected_point_vector- atom_vector)
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

