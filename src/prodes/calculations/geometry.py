import numpy as np

from prodes.calculations.distance_functions import atom_charge_coulomb, potential_multiple_media
from math import sqrt


# inspired by https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
class Sunflower_sphere:

    def __init__(self, x, y, z, r, n_points):
        
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.n_points = n_points
        self._points = None

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

            points = np.empty([0])
            for i in range(len(x)):
                point = Point(x[i], y[i], z[i])
                points = np.concatenate([points, np.array([point])])

            self._points = points

        return self._points


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


def maximal_distance(normal_vector, vector_on_plane, points):
    """Ditermines the maximal distance by calculating the dot product
    of each point in the system and the unit normal vector"""

    unit_vector = normal_vector/np.linalg.norm(normal_vector)
    maximum = 0
    for point in points:
        point_vector = make_vector(point)
        dot_prod = np.dot(unit_vector, point_vector - vector_on_plane)
        if dot_prod > maximum:
            maximum = dot_prod

    return maximum


def required_distance(point_for_plane, structure, surface_points):
    """ditermines the required distance for the plane to be formed on the protein surface"""

    vector_on_plane = make_vector(structure)
    normal_vector = make_vector(point_for_plane) - vector_on_plane

    return maximal_distance(normal_vector, vector_on_plane, surface_points) + 1
    

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

def map_ep_to_plane(atom, projected_point_vector, surface_exit, ph=7):
    """Calculates the electrostatic potential of an atom projected onto a plane"""
    atom_vector = make_vector(atom)

    total_distance = np.linalg.norm(projected_point_vector- atom_vector)
    protein_distance = np.linalg.norm(surface_exit - atom_vector)
    
    distances = [(total_distance - protein_distance) * 10**-10, protein_distance * 10**-10]
    atom_charge = atom_charge_coulomb(atom.charge(ph))
    
    point_potential = potential_multiple_media(atom_charge, distances, [80, 4])

    return point_potential

