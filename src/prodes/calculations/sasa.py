import numpy as np

from prodes.calculations.geometry import Sunflower_sphere
from prodes.calculations.grid_wizard import Grid
from prodes.calculations.standard_equations import distance


def shrake_rupley(grid:Grid, probe_r=1.4, consider=["All"]):
    """takes a grid of atoms and will return a Structure object containing dummy atoms"""

    from prodes.core.point import Surface_point
    from math import pi
    surface = np.empty([0])
    grid_cells = grid.cells.flatten()
    for cell in grid_cells:

        if not cell.empty:
            surounding_cells = grid.find_surrounding_cells(cell)
            neighbourhood = grid.grid_content("Atom", cells=surounding_cells)
            for atom in cell.filtered_content("Atom"):
                if atom.radius is None or ("All" not in consider and atom.element not in consider):
                    pass

                else:
                    radius = atom.radius + probe_r
                    sphere_points = Sunflower_sphere(atom.x, atom.y, atom.z, radius, int(radius**2*4*pi)*2)

                    # finding a sub neighbourhood for neighbouringing atoms
                    sub_neighbourhood = np.empty([0])
                    for neighbour in neighbourhood:

                        if neighbour.radius is not None:
                            dist = distance(atom, neighbour)
                            neighbour_radius = neighbour.radius + probe_r
                            if dist < neighbour_radius + radius:
                                sub_neighbourhood = np.concatenate([sub_neighbourhood, np.array([neighbour])])

                    point_nr = 0
                    surf_points = 0
                    cloud = []
                    for point in sphere_points.points:
                        point_nr += 1
                        on_surf = True
                        for neighbour in sub_neighbourhood:
                            if atom != neighbour:
                                neighbour_radius = neighbour.radius + probe_r
                                dist = distance(point, neighbour)

                                if dist < neighbour_radius:
                                    on_surf = False
                                    break
                        if on_surf:
                            surf_points += 1
                            surf_point = Surface_point(point.x, point.y, point.z, atom)
                            cloud.append(surf_point)
                            surface = np.concatenate([surface, np.array([surf_point])])

                    atom.exposed = surf_points/point_nr
                    atom.cloud = np.array(cloud)

    atom.structure.surface_done = True
    return surface


def shape(surface_points, structure):
    """Describes the shape of a protein surface"""

    max_dist = 0
    min_dist = 10e100
    total_dist = 0

    for point in surface_points:
        dist = distance(point, structure)
        total_dist += dist
        if dist > max_dist:
            max_dist = dist

        if dist < min_dist:
            min_dist = dist 

    mean_dist = total_dist/len(surface_points)

    return mean_dist/max_dist, min_dist/mean_dist
