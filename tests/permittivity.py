from prodes.core.point import Point
from prodes.io.parser import PDBparser, Builder, write_pdb
from prodes.calculations.grid_wizard import Grid, property_points_on_surface
from prodes.calculations.sasa import shrake_rupley
from prodes.core.structure import Structure
from prodes.core.point import Point, Property_point

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

import local_descriptors as ld

structure = PDBparser().parse("package/tests/data/1GDW.pdb")
grid = Grid(12)
grid.construct_cells(structure.atoms)
grid.fill_cells(structure.atoms)

surface = shrake_rupley(grid)

small_grid = Grid(1)
small_grid.construct_cells(surface)
small_grid.fill_cells(surface)

property_points = property_points_on_surface(small_grid)

to_plot = [point for point in property_points]
# surf_struct = Structure()
# property_dummies = []
for point in property_points:
    point.set_ep(structure.heavy_atoms, 8)
#     print(point.ep)
#     dummy = Builder().build_dummy_atom(point.x, point.y, point.z)
#     property_dummies.append(dummy)
# surf_struct.atoms = np.array(property_dummies)
    
# write_pdb(surf_struct, f"1GDW_surf.pdb")



point = Point(structure.x, structure.y, structure.z-1)
distance = ld.required_distance(point, structure, surface)
ld.move_point(point, structure, distance)  

plane = ld.find_plane(point, structure)

# plane_structure = Structure("plane")
for atom in np.array([atom for atom in structure.atoms if atom.charge() != 0]):
    projected_point = ld.project_point(*plane, *ld.make_vector(atom))
    atom_vector = ld.make_vector(atom)

    normal_vector = projected_point - atom_vector
    highest = 0
    for surface_point in property_points:
        
        surface_point_vector = ld.make_vector(surface_point)- atom_vector
        dot_prod = np.dot(normal_vector/np.linalg.norm(normal_vector), surface_point_vector)
        if dot_prod > 0:
       
            

            potential_pass = normal_vector/np.linalg.norm(normal_vector)
            distance = np.linalg.norm(potential_pass-surface_point_vector)

            
            if distance < 2:
                if dot_prod > highest:
                    highest = dot_prod
    print(highest)
    to_plot.append(Property_point(*projected_point,-4))
#     dummy = Builder().build_dummy_atom(*projected_point)
#     projected_points.append(dummy)

# plane_structure.atoms = np.array(projected_points)






fig = plt.figure()
ax = fig.add_subplot(projection="3d")

f = ax.scatter3D([point.x for point in to_plot], [point.y for point in to_plot], [point.z for point in to_plot], c=[point.ep for point in to_plot], cmap='coolwarm', s=1)
fig.colorbar(f, ax=ax)
plt.axis('off')
plt.savefig("splane.png", dpi=300)