import numpy as np

from prodes.calculations.geometry import Sunflower_sphere
from prodes.calculations.grid_wizard import Grid
from prodes.calculations.standard_equations import distance


def shrake_rupley_original(grid:Grid, probe_r=1.4, consider=["All"]):
    """takes a grid of atoms and will return a Structure object containing dummy atoms
    NOTE: This function is not optimized for large grids And has been left in the repo for reference only.
    In order to verify that the vectorization gave the same result as the original function.
    For production use, use the optimized version shrake_rupley() below.
    """

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


def shrake_rupley(grid:Grid, probe_r=1.4, consider=["All"]):
    """Vectorized Shrake-Rupley SASA calculation.

    Takes a grid of atoms and returns a numpy array of Surface_point objects.
    Uses numpy broadcasting for the point-neighbour distance computation instead
    of per-point Python loops, giving identical results with significant speedup.
    """

    from prodes.core.point import Surface_point
    from math import pi
    surface_list = []
    grid_cells = grid.cells.flatten()
    structure_ref = None

    for cell in grid_cells:

        if not cell.empty:
            surounding_cells = grid.find_surrounding_cells(cell)
            neighbourhood = grid.grid_content("Atom", cells=surounding_cells)
            for atom in cell.filtered_content("Atom"):
                if atom.radius is None or ("All" not in consider and atom.element not in consider):
                    continue

                if structure_ref is None:
                    structure_ref = atom.structure

                radius = atom.radius + probe_r
                sphere = Sunflower_sphere(atom.x, atom.y, atom.z, radius, int(radius**2*4*pi)*2)
                sphere_coords = sphere.raw_coordinates
                n_points = sphere_coords.shape[0]

                # Build sub_neighbourhood: neighbors with valid radii, excluding self,
                # whose sphere could overlap with this atom's sphere
                valid_neighbors = [n for n in neighbourhood
                                   if n.radius is not None and n is not atom]
                if len(valid_neighbors) == 0:
                    on_surf_mask = np.ones(n_points, dtype=bool)
                else:
                    neighbor_coords = np.array([[n.x, n.y, n.z] for n in valid_neighbors])
                    neighbor_radii = np.array([n.radius + probe_r for n in valid_neighbors])

                    # Vectorized sub_neighbourhood filter: dist < neighbour_radius + radius
                    atom_coord = np.array([atom.x, atom.y, atom.z])
                    atom_neighbor_dists = np.sqrt(np.sum((neighbor_coords - atom_coord)**2, axis=1))
                    sub_mask = atom_neighbor_dists < (neighbor_radii + radius)

                    if not np.any(sub_mask):
                        on_surf_mask = np.ones(n_points, dtype=bool)
                    else:
                        sub_coords = neighbor_coords[sub_mask]
                        sub_radii = neighbor_radii[sub_mask]

                        # Compute all point-neighbour distances via broadcasting
                        # sphere_coords: (N, 3), sub_coords: (M, 3) -> diff: (N, M, 3)
                        diff = sphere_coords[:, np.newaxis, :] - sub_coords[np.newaxis, :, :]
                        dists = np.sqrt(np.sum(diff**2, axis=2))  # (N, M)

                        # Point is on surface if NOT inside any neighbour's sphere
                        on_surf_mask = np.all(dists >= sub_radii[np.newaxis, :], axis=1)

                surf_indices = np.where(on_surf_mask)[0]
                atom.exposed = len(surf_indices) / n_points

                cloud = []
                for i in surf_indices:
                    surf_point = Surface_point(
                        sphere_coords[i, 0], sphere_coords[i, 1], sphere_coords[i, 2], atom
                    )
                    cloud.append(surf_point)
                    surface_list.append(surf_point)

                atom.cloud = np.array(cloud)

    if structure_ref is not None:
        structure_ref.surface_done = True

    surface = np.array(surface_list) if surface_list else np.empty([0])
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
