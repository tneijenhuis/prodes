import numpy as np

from prodes.calculations.geometry import Sunflower_sphere
from prodes.calculations.grid_wizard import Grid


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

    coords = np.array([[p.x, p.y, p.z] for p in surface_points])
    center = np.array([structure.x, structure.y, structure.z])
    dists = np.sqrt(np.sum((coords - center)**2, axis=1))
    mean_dist = dists.mean()
    max_dist = dists.max()
    min_dist = dists.min()

    return mean_dist/max_dist, min_dist/mean_dist
