
import numpy as np
from prodes.io.parser import PDBparser, Builder, write_pdb
from prodes.calculations.geometry import Sunflower_sphere
from prodes.calculations.grid_wizard import Grid, property_points_on_surface
from prodes.calculations.sasa import shrake_rupley


from prodes.core.structure import Structure

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

    return maximal_distance(normal_vector, vector_on_plane, surface_points)
    

def project_point(a, b, c, d, x1, y1, z1):
    """projects point 1 onto plane ax+by+cz=-d"""
      
    k =(d -a * x1-b * y1-c * z1)/(a * a + b * b + c * c)
    x2 = a * k + x1
    y2 = b * k + y1
    z2 = c * k + z1
    return x2, y2, z2

def main():
    structure = PDBparser().parse("package/tests/data/1GDW.pdb")
    grid = Grid(12)
    grid.construct_cells(structure.atoms)
    grid.fill_cells(structure.atoms)

    surface = shrake_rupley(grid)

    small_grid = Grid(1)
    small_grid.construct_cells(surface)
    small_grid.fill_cells(surface)

    property_points = property_points_on_surface(small_grid)

    surf_struct = Structure()
    property_dummies = []
    for point in property_points:
        dummy = Builder().build_dummy_atom(*make_vector(point))
        property_dummies.append(dummy)
    surf_struct.atoms = np.array(property_dummies)
        
    write_pdb(surf_struct, f"1GDW_surf.pdb")
    surrounding_points = Sunflower_sphere(*make_vector(structure), 1, 120).points
        
    for i, point in enumerate(surrounding_points):    
        distance = required_distance(point, structure, surface)

        move_point(point, structure, distance)  

        plane = find_plane(point, structure)

        projected_points = []
        plane_structure = Structure("plane")
        for atom in structure.atoms:
            projected_point = project_point(*plane, *make_vector(atom))
            dummy = Builder().build_dummy_atom(*projected_point)
            projected_points.append(dummy)

        plane_structure.atoms = np.array(projected_points)
        
        # write_pdb(plane_structure, f"plane_{i}.pdb")

if __name__ == "__main__":
    main()