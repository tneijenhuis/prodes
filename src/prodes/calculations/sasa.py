from math import pi

import numpy as np

from prodes.calculations.geometry import Sunflower_sphere
from prodes.calculations.grid_wizard import Grid
from prodes.parallel import run_tasks

# Task descriptors for process_atom, populated by the parent before the worker
# pool is created so that the forked children inherit them without pickling.
# Cleared as soon as the tasks have run, and only ever written by one thread:
# see the note on thread safety in shrake_rupley.
SASA_TASKS: dict = {}


def sasa_task_data(grid: Grid, probe_r=1.4, consider=("All",)):
    """Collects the per-atom task descriptors for a Shrake-Rupley calculation.

    Walks the grid once in the parent and reduces every atom to plain numpy
    arrays, so that the workers never touch an Atom object. Atom carries a
    back-reference to its Structure, and pickling one atom would drag the whole
    structure along with it.

    Every atom in a cell sees the same neighbourhood, so the neighbour arrays are
    built once per cell and shared by every task in that cell. Building them per
    atom instead costs both time and memory in proportion to the atoms per cell,
    which is ~27 on a typical structure at the default grid size: on 1GPB that
    measured 1.6 s and 156 MB against 0.09 s and 4.5 MB. This step runs in the
    parent before the fork and cannot be parallelised, so it sets the floor on
    how fast the SASA phase can get.

    Args:
        grid: a Grid already filled with the atoms to consider.
        probe_r: radius of the solvent probe.
        consider: element symbols to include, or ("All",) for every element.

    Returns:
        atoms: the Atom objects in task order, used by the parent for reassembly.
        neighbourhoods: one (neighbour_coords, neighbour_radii) pair per cell that
            holds at least one atom to process.
        tasks: one (neighbourhood_index, self_index, atom_coord, total_radius)
            tuple per atom, in the same order as atoms. self_index is the atom's
            own row in its neighbourhood arrays, or -1 when it is absent.
    """

    atoms: list = []
    tasks: list = []
    neighbourhoods: list = []

    for cell in grid.cells.flatten():

        if cell.empty:
            continue

        selected = [atom for atom in cell.filtered_content("Atom") if atom.radius is not None and ("All" in consider or atom.element in consider)]
        if not selected:
            continue

        surrounding_cells = grid.find_surrounding_cells(cell)
        neighbourhood = grid.grid_content("Atom", cells=surrounding_cells)

        valid_neighbors = [n for n in neighbourhood if n.radius is not None]
        if len(valid_neighbors) == 0:
            neighbour_coords = np.empty((0, 3))
            neighbour_radii = np.empty(0)
        else:
            neighbour_coords = np.array([[n.x, n.y, n.z] for n in valid_neighbors])
            neighbour_radii = np.array([n.radius + probe_r for n in valid_neighbors])

        # Row of each atom in the shared arrays, so that a task can mask out its
        # own entry. Keyed by identity because two atoms can share coordinates.
        rows = {id(n): row for row, n in enumerate(valid_neighbors)}

        neighbourhood_index = len(neighbourhoods)
        neighbourhoods.append((neighbour_coords, neighbour_radii))

        for atom in selected:
            atoms.append(atom)
            tasks.append((neighbourhood_index, rows.get(id(atom), -1), np.array([atom.x, atom.y, atom.z]), atom.radius + probe_r))

    return atoms, neighbourhoods, tasks


def process_atom(index: int):
    """Calculates the exposed fraction and surface coordinates of a single atom.

    Reads its task descriptor from SASA_TASKS, which the parent fills in before
    the workers are forked. Returns raw coordinates rather than Surface_point
    objects because those carry a back-reference to their atom, which would make
    them expensive to send back to the parent.

    Returns:
        exposed: the fraction of the atom sphere that is solvent accessible.
        surface_coords: (N, 3) array of the accessible sphere points.
    """

    neighbourhood_index, self_index, atom_coord, total_radius = SASA_TASKS["tasks"][index]
    neighbour_coords, neighbour_radii = SASA_TASKS["neighbourhoods"][neighbourhood_index]

    sphere = Sunflower_sphere(atom_coord[0], atom_coord[1], atom_coord[2], total_radius, int(total_radius**2 * 4 * pi) * 2)
    sphere_coords = sphere.raw_coordinates
    n_points = sphere_coords.shape[0]

    if len(neighbour_coords) == 0:
        on_surf_mask = np.ones(n_points, dtype=bool)
    else:
        # Vectorized sub_neighbourhood filter: dist < neighbour_radius + radius
        atom_neighbor_dists = np.sqrt(np.sum((neighbour_coords - atom_coord) ** 2, axis=1))
        sub_mask = atom_neighbor_dists < (neighbour_radii + total_radius)

        # The atom shares its neighbourhood arrays with the rest of its cell, so
        # it has to drop its own row here. Keeping it would be a no-op in exact
        # arithmetic, since a sphere point sits at exactly total_radius from its
        # own centre and the test below is >=, but the recomputed distance can
        # land a few ULP short and discard a genuine surface point.
        if self_index >= 0:
            sub_mask[self_index] = False

        if not np.any(sub_mask):
            on_surf_mask = np.ones(n_points, dtype=bool)
        else:
            sub_coords = neighbour_coords[sub_mask]
            sub_radii = neighbour_radii[sub_mask]

            # Compute all point-neighbour distances via broadcasting
            # sphere_coords: (N, 3), sub_coords: (M, 3) -> diff: (N, M, 3)
            diff = sphere_coords[:, np.newaxis, :] - sub_coords[np.newaxis, :, :]
            dists = np.sqrt(np.sum(diff**2, axis=2))  # (N, M)

            # Point is on surface if NOT inside any neighbour's sphere
            on_surf_mask = np.all(dists >= sub_radii[np.newaxis, :], axis=1)

    surf_indices = np.where(on_surf_mask)[0]

    return len(surf_indices) / n_points, sphere_coords[surf_indices]


def shrake_rupley(grid: Grid, probe_r=1.4, consider=("All",)):
    """Vectorized Shrake-Rupley SASA calculation.

    Takes a grid of atoms and returns a numpy array of Surface_point objects.
    Uses numpy broadcasting for the point-neighbour distance computation instead
    of per-point Python loops, and spreads the atoms over worker processes where
    the platform allows it.

    Not thread-safe: the tasks are handed to the workers through the module-level
    SASA_TASKS, so two threads calling this at once would overwrite each other's
    tasks. Run one calculation per process, or serialise the calls.
    """

    from prodes.core.point import Surface_point

    atoms, neighbourhoods, tasks = sasa_task_data(grid, probe_r, consider)
    if not atoms:
        return np.empty([0])

    SASA_TASKS["neighbourhoods"] = neighbourhoods
    SASA_TASKS["tasks"] = tasks
    try:
        results = run_tasks(process_atom, len(tasks), "Shrake-Rupley SASA")
    finally:
        # Nothing reads the tasks after this point, and they are the largest
        # thing this module holds. Leaving them in place would keep one
        # structure's arrays alive until the next call overwrites them.
        SASA_TASKS.clear()

    surface_list = []
    for atom, (exposed, surface_coords) in zip(atoms, results, strict=True):
        atom.exposed = exposed

        cloud = [Surface_point(x, y, z, atom) for x, y, z in surface_coords]
        surface_list.extend(cloud)
        atom.cloud = np.array(cloud)

    structure_ref = atoms[0].structure
    if structure_ref is not None:
        structure_ref.surface_done = True

    surface = np.array(surface_list) if surface_list else np.empty([0])
    return surface


def shape(surface_points, structure):
    """Describes the shape of a protein surface"""

    coords = np.array([[p.x, p.y, p.z] for p in surface_points])
    center = np.array([structure.x, structure.y, structure.z])
    dists = np.sqrt(np.sum((coords - center) ** 2, axis=1))
    mean_dist = dists.mean()
    max_dist = dists.max()
    min_dist = dists.min()

    return mean_dist / max_dist, min_dist / mean_dist
