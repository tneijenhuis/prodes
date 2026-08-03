import argparse
import os

import numpy as np
import pandas as pd
from dotenv import load_dotenv

from prodes.calculations import geometry, grid_wizard
from prodes.calculations.sasa import shape, shrake_rupley
from prodes.calculations.standard_equations import trimean
from prodes.io import parser as ps
from prodes.parallel import run_tasks, validate_settings, worker_mem_limit_mb

load_dotenv()

# Shared state for the worker functions below. The parent fills these in before
# creating a worker pool, so that the forked children inherit the grid, atom and
# structure objects through copy-on-write instead of having them pickled. The
# workers are keyed by an integer index and return plain values only.
#
# Each dict is cleared as soon as its tasks have run, so that a structure's grid
# does not stay alive for the rest of the process. Because they are module-level,
# calculate() is safe to run in several processes but not in several threads:
# see the note in its docstring.
SURFACE_GRID_STATE: dict = {}
AVERAGE_CHARGE_GRID_STATE: dict = {}
SHELL_STATE: dict = {}


def parse_arguments():
    """Parses all arguments given in the commandline"""

    env_full = os.getenv("PRODES_FULL_FEATURES", "false").lower() in ("true", "1", "yes")

    parser = argparse.ArgumentParser(description="Calculate descriptors from atomic data")
    parser.add_argument("pdb_file", help="file location of a pdb or pqr file", type=str)
    parser.add_argument("out_file", help="file path of the output csv file", type=str)
    parser.add_argument("-p", "--pka", help="file location of the pka propka output", type=str, default=None)
    parser.add_argument("--probe", help="Radius of the surface probe", type=float, default=1.4)
    parser.add_argument("--ph", help="pH of the system", type=float, default=7)
    parser.add_argument("--hydro", help="Abbreviation of the hydrophobicity scale to be used", type=str, default="mj_scaled")
    parser.add_argument(
        "--full-features",
        action=argparse.BooleanOptionalAction,
        default=env_full,
        help="Calculate the full legacy feature set including redundant features (default: from PRODES_FULL_FEATURES env var)",
    )
    parser.add_argument(
        "--mem-limit",
        type=float,
        default=None,
        help="Memory budget in MB for the intermediate NumPy arrays of the whole run, divided among the workers (default: from PRODES_MEM_LIMIT_MB env var, or 2048)",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Number of worker processes, Linux only (default: from PRODES_N_WORKERS env var, or half the logical CPUs)",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=None,
        help="Tasks dispatched per chunk to each worker, advanced tuning (default: from PRODES_CHUNKSIZE env var, or 1)",
    )

    arg = parser.parse_args()

    # The worker settings are read from the environment at the point of use, deep
    # inside the first phase, so a CLI flag is applied by setting the variable it
    # overrides. Validating all of them here means a typo fails in the first
    # second rather than after the grid has been built.
    if arg.n_workers is not None:
        os.environ["PRODES_N_WORKERS"] = str(arg.n_workers)
    if arg.chunksize is not None:
        os.environ["PRODES_CHUNKSIZE"] = str(arg.chunksize)
    validate_settings(arg.mem_limit)

    pdb_file = arg.pdb_file
    out_file = arg.out_file
    pkas_file = arg.pka
    ph = arg.ph
    r_probe = arg.probe
    hydro_scale = arg.hydro
    full_features = arg.full_features
    mem_limit_mb = arg.mem_limit

    return pdb_file, out_file, pkas_file, ph, r_probe, hydro_scale, full_features, mem_limit_mb


def open_output_file(out_file):
    """opens the output file"""
    return pd.read_csv(out_file)


def write_output_file(dataframe, out_file):
    """writes a CSV file"""
    dataframe.to_csv(out_file, index=False)


def standard_features(values, name="", reduced=False):
    """calculates central-tendency and spread statistics of a np array of values.

    When reduced=True, only Mean and Std are returned (Trimean, Median, and Sum
    are dropped as they are near-perfectly correlated with Mean or recoverable
    from Mean x N).
    """
    features = {}

    if len(values) != 0:
        features[f"{name}Mean"] = round(values.mean(), 3)
        if not reduced:
            features[f"{name}Trimean"] = round(trimean(values), 3)
            features[f"{name}Median"] = round(np.median(values), 3)
            features[f"{name}Sum"] = round(np.sum(values), 3)
        features[f"{name}Std"] = round(values.std(), 3)

    else:
        features[f"{name}Mean"] = 0
        if not reduced:
            features[f"{name}Trimean"] = 0
            features[f"{name}Median"] = 0
            features[f"{name}Sum"] = 0
        features[f"{name}Std"] = 0

    return features


def cells_with_property_points(grid):
    """Returns the grid cells that hold at least one Property_point.

    Cells without property points contribute nothing, so collecting them up
    front turns the grid traversal into a flat list of equally shaped tasks.
    """

    return [cell for cell in grid.cells.flatten() if len(cell.filtered_content("Property_point")) > 0]


def process_surface_grid_cell(index: int) -> list:
    """Calculates the electrostatic potential and lipophilicity of one cell's property points.

    Reads the grid and charged atoms from SURFACE_GRID_STATE. Returns (ep, lipo)
    tuples rather than the points themselves, because a forked worker mutates
    its own copy of the point and the parent has to write the values back onto
    the originals.
    """

    grid = SURFACE_GRID_STATE["grid"]
    cell = SURFACE_GRID_STATE["cells"][index]

    environment = np.array([atoms for atoms in grid.grid_content("Atom", cells=grid.find_surrounding_cells(cell))])

    values = []
    for point in cell.filtered_content("Property_point"):
        point.set_ep(SURFACE_GRID_STATE["charged_atoms"], SURFACE_GRID_STATE["ph"])
        point.set_lipo(environment, 10, SURFACE_GRID_STATE["hydro_scale"])
        values.append((point.ep, point.lipo))

    return values


def calculate_surface_grid_features(structure, surface_points, ph, hydro_scale, features: dict, full_features=True):
    """calculates the features from the surface grid"""

    reduced = not full_features

    if full_features:
        features["NSurfPoints"] = len(surface_points)
    surf_shape_max, surf_shape_min = shape(surface_points, structure)
    features["Shape max"] = round(surf_shape_max, 3)
    features["Shape min"] = round(surf_shape_min, 3)

    concat = np.concatenate([structure.heavy_atoms, surface_points])
    grid = grid_wizard.Grid(12)
    grid.construct_cells(concat)
    grid.fill_cells(concat)

    cells = cells_with_property_points(grid)
    SURFACE_GRID_STATE.update(
        grid=grid,
        cells=cells,
        charged_atoms=np.array([atom for atom in structure.atoms if atom.charge(ph=ph) != 0]),
        ph=ph,
        hydro_scale=hydro_scale,
    )
    try:
        cell_values = run_tasks(process_surface_grid_cell, len(cells), "surface grid features")
    finally:
        SURFACE_GRID_STATE.clear()

    for cell, values in zip(cells, cell_values, strict=True):
        for point, (ep, lipo) in zip(cell.filtered_content("Property_point"), values, strict=True):
            point.set_values(ep=ep, lipo=lipo)

    # Electrostatic potential features
    eps = np.array([point.ep for point in surface_points])

    features["SurfEpMaxFormal"] = round(eps.max(), 3)
    features["SurfEpMinFormal"] = round(eps.min(), 3)
    features.update({f"{v}Formal": k for v, k in standard_features(eps, "SurfEp", reduced=reduced).items()})

    positive_eps = np.array([ep for ep in eps if ep > 0])
    features["NSurfPosEpFormal"] = len(positive_eps)
    features.update({f"{v}Formal": k for v, k in standard_features(positive_eps, "SurfEpPos", reduced=reduced).items()})

    negative_eps = np.array([ep for ep in eps if ep < 0])
    # NSurfNegEpFormal = NSurfPoints - NSurfPosEpFormal, near-exact: the two counts
    # partition the same point set, apart from points whose potential is exactly zero.
    # NSurfPoints itself is R2 0.9998 with Area, so NSurfNegEpFormal is 0.9998
    # predicted by Area + NSurfPosEpFormal, both of which are kept. NSurfNegEpFormal
    # is also the size proxy of the pair (R2 0.87 with Area, vs 0.03 for
    # NSurfPosEpFormal), so NSurfPosEpFormal is the one that carries signal.
    if full_features:
        features["NSurfNegEpFormal"] = len(negative_eps)
    features.update({f"{v}Formal": k for v, k in standard_features(negative_eps, "SurfEpNeg", reduced=reduced).items()})

    # Hydrophobic potential features
    lipos = np.array([point.lipo for point in surface_points])

    features["SurfMhpMax"] = round(lipos.max(), 3)
    features["SurfMhpMin"] = round(lipos.min(), 3)
    features.update(standard_features(lipos, "SurfMhp", reduced=reduced))

    positive_lipos = np.array([lipo for lipo in lipos if lipo > 0])
    features["NSurfPosMhp"] = len(positive_lipos)
    features.update(standard_features(positive_lipos, "SurfPosMhp", reduced=reduced))

    negative_lipos = np.array([lipo for lipo in lipos if lipo < 0])
    # NSurfNegMhp = NSurfPoints - NSurfPosMhp, exactly: the two counts partition the
    # same point set. NSurfPoints itself is R2 0.9998 with Area, so NSurfNegMhp is a
    # linear combination of Area + NSurfPosMhp, both of which are kept, and carries no
    # information of its own. SurfNegMhpMean and SurfNegMhpStd below are not size
    # proxies and stay.
    if full_features:
        features["NSurfNegMhp"] = len(negative_lipos)
    features.update(standard_features(negative_lipos, "SurfNegMhp", reduced=reduced))

    return features


def process_average_charge_grid_cell(index: int) -> list:
    """Calculates the average-charge electrostatic potential of one cell's property points.

    The average-charge counterpart of process_surface_grid_cell: it uses partial
    rather than formal charges and does not project lipophilicity, so it needs
    neither the grid neighbourhood nor a hydrophobicity scale.
    """

    cell = AVERAGE_CHARGE_GRID_STATE["cells"][index]

    values = []
    for point in cell.filtered_content("Property_point"):
        point.set_ep(AVERAGE_CHARGE_GRID_STATE["charged_atoms"], AVERAGE_CHARGE_GRID_STATE["ph"], formal=False)
        values.append(point.ep)

    return values


def calculate_average_chargesurface_grid_features(structure, surface_points, ph, features: dict):

    concat = np.concatenate([structure.heavy_atoms, surface_points])
    grid = grid_wizard.Grid(12)
    grid.construct_cells(concat)
    grid.fill_cells(concat)

    cells = cells_with_property_points(grid)
    AVERAGE_CHARGE_GRID_STATE.update(
        cells=cells,
        charged_atoms=np.array([atom for atom in structure.atoms if atom.charge(ph=ph, formal=False) != 0]),
        ph=ph,
    )
    try:
        cell_values = run_tasks(process_average_charge_grid_cell, len(cells), "average charge surface grid features")
    finally:
        AVERAGE_CHARGE_GRID_STATE.clear()

    for cell, values in zip(cells, cell_values, strict=True):
        for point, ep in zip(cell.filtered_content("Property_point"), values, strict=True):
            point.set_values(ep=ep)

    # Electrostatic potential features
    eps = np.array([point.ep for point in surface_points])

    features["SurfEpMaxAverage"] = round(eps.max(), 3)
    features["SurfEpMinAverage"] = round(eps.min(), 3)
    features.update({f"{v}Average": k for v, k in standard_features(eps, "SurfEp").items()})

    positive_eps = np.array([ep for ep in eps if ep > 0])
    features["NSurfPosEpAverage"] = len(positive_eps)
    features.update({f"{v}Average": k for v, k in standard_features(positive_eps, "SurfEpPos").items()})

    negative_eps = np.array([ep for ep in eps if ep < 0])
    features["NSurfNegEpAverage"] = len(negative_eps)
    features.update({f"{v}Average": k for v, k in standard_features(negative_eps, "SurfEpNeg").items()})
    return features


def charged_atom_arrays(atoms, ph, formal=True):
    """Returns the coordinates and charges of the atoms that actually carry a charge.

    Atoms with zero charge contribute nothing to an electrostatic potential, so
    dropping them here keeps the batch geometry calculations to the smallest
    array that gives the same answer.

    Args:
        atoms: iterable of Atom objects.
        ph: pH at which the charge is evaluated.
        formal: use formal charges when True and partial charges when False.

    Returns:
        coords: (M, 3) array of atom positions.
        charges: (M,) array of charges in elementary charge units.
    """

    charged = [atom for atom in atoms if atom.charge(ph, formal=formal) != 0]
    if not charged:
        return np.empty((0, 3)), np.empty(0)

    return np.array([[a.x, a.y, a.z] for a in charged]), np.array([a.charge(ph, formal=formal) for a in charged])


def process_shell_plane(index: int) -> float:
    """Calculates the summed electrostatic potential mapped onto a single shell plane.

    Reads the structure, surface and charged atom arrays from SHELL_STATE.
    move_point rewrites the plane point in place, but in a forked worker that
    happens to the child's own copy and never reaches the parent, which needs
    only the returned sum.
    """

    structure = SHELL_STATE["structure"]
    surface_coords = SHELL_STATE["surface_coords"]
    charged_coords = SHELL_STATE["charged_coords"]

    point = SHELL_STATE["plane_points"][index]
    distance = geometry.required_distance(point, structure, surface_coords)
    geometry.move_point(point, structure, distance)

    a, b, c, d = geometry.find_plane(point, structure)
    projected_coords = geometry.project_point_batch(a, b, c, d, charged_coords)
    exits, has_exit = geometry.find_exit_batch(charged_coords, projected_coords, surface_coords, mem_limit_mb=SHELL_STATE["mem_limit_mb"])
    potentials = geometry.map_ep_to_plane_batch(charged_coords, projected_coords, exits, has_exit, SHELL_STATE["charged_charges"])

    return float(potentials.sum())


def calculate_shell_features(structure, surface_points, ph: float, features: dict, numb_of_planes=120, full_features=True, mem_limit_mb=None):
    """Constructs a number of planes onto which charges are mapped

    This is the phase that allocates the large intermediate arrays, in
    find_exit_batch. Every worker runs it at the same time, so the memory budget
    is divided among them here rather than being handed to each worker whole.
    """

    reduced = not full_features

    surface_coords = np.array([[p.x, p.y, p.z] for p in surface_points])
    plane_points = geometry.Sunflower_sphere(structure.x, structure.y, structure.z, 1, numb_of_planes).points
    charged_coords, charged_charges = charged_atom_arrays(structure.heavy_atoms, ph)

    SHELL_STATE.update(
        structure=structure,
        surface_coords=surface_coords,
        plane_points=plane_points,
        charged_coords=charged_coords,
        charged_charges=charged_charges,
        mem_limit_mb=worker_mem_limit_mb(mem_limit_mb),
    )
    try:
        shell_potentials = np.array(run_tasks(process_shell_plane, len(plane_points), "shell features"))
    finally:
        SHELL_STATE.clear()

    features["ShellEpMaxFormal"] = round(shell_potentials.max(), 3)
    features["ShellEpminFormal"] = round(shell_potentials.min(), 3)
    features.update({f"{v}Formal": k for v, k in standard_features(shell_potentials, "ShellEp", reduced=reduced).items()})

    pos_potentials = np.array([potential for potential in shell_potentials if potential > 0])
    features["NShellPosEpFormal"] = len(pos_potentials)
    features.update({f"{v}Formal": k for v, k in standard_features(pos_potentials, "ShellEpPos", reduced=reduced).items()})

    neg_potentials = np.array([potential for potential in shell_potentials if potential < 0])
    if full_features:
        features["NShellNegEpFormal"] = len(neg_potentials)
    features.update({f"{v}Formal": k for v, k in standard_features(neg_potentials, "ShellEpNeg", reduced=reduced).items()})

    return features


def calculate_structure_features(structure, ph, r_probe, features: dict, full_features=True):
    """calculates general structure features, including"""

    features["Molecular weight"] = structure.mw
    features["Isoelectric point"] = structure.isoelectric_point()
    features["Dipole"] = structure.dipole(ph)
    features["Formal charge"] = structure.charge(ph)
    if full_features:
        features["Average charge"] = structure.charge(ph, formal=False)

    if structure.surface_done is False:
        grid_size = 10 + (r_probe - 1.4) * 2
        grid = grid_wizard.Grid(grid_size)
        grid.construct_cells(structure.heavy_atoms)
        grid.fill_cells(structure.heavy_atoms)

        shrake_rupley(grid, r_probe)

    features["Area"] = round(structure.surface_area(r_probe), 3)

    fractions = structure.residue_surf_fractions(r_probe)

    for residue, fraction in fractions.items():
        features[f"{residue}SurfFrac"] = fraction

    return features


def prepare_structure(pdb_file, pkas_file):
    """loads and prepares the structure for the calculations"""

    structure = ps.PDBparser().parse(pdb_file)
    if pkas_file:
        pkas = ps.read_pka(pkas_file)
        structure.redo_pkas(pkas)

    return structure


def construct_surface_grid(structure, r_probe):
    """Constructs the protein surface grid"""

    grid_size = 10 + (r_probe - 1.4) * 2
    grid = grid_wizard.Grid(grid_size)
    grid.construct_cells(structure.heavy_atoms)
    grid.fill_cells(structure.heavy_atoms)

    surface = shrake_rupley(grid, r_probe)

    surface_grid = grid_wizard.Grid(1)
    surface_grid.construct_cells(surface)
    surface_grid.fill_cells(surface)

    property_points = grid_wizard.property_points_on_surface(surface_grid)

    return property_points


def calculate(
    pdb_file,
    out_file,
    pkas_file=None,
    ph=7,
    r_probe=1.4,
    hydro_scale="mj_scaled",
    full_features=False,
    mem_limit_mb=None,
):
    """Calculates a list of supported features and returns a csv file

    Arguments
    required:
        pdb_file: path to the pdb file to be analysed
        out_file: path to the output csv file which should be written
    Optional:
        pka_file: an output file of PROPKA which is required for custom pKa assignment
        ph: the ph at which protonation states should be calculated
        r_probe: the radius of the probe used to calculate the solvent accessible surface area
        hydro_scale: the abbreviation of the hydrophobicity scale used (scales can be found in data/hydrophobicity)
        full_features: when True, calculates the full legacy 105-feature set including
            redundant features. Defaults to False for the reduced, non-redundant set.
        mem_limit_mb: memory budget in MB for the intermediate NumPy arrays of this
            run, divided among the workers. If None, falls back to the
            PRODES_MEM_LIMIT_MB env var (default 2048). Pass explicitly when calling
            as a library to avoid relying on environment variables.

    Thread safety
        Safe to run in several processes, one structure each, which is how the
        pipeline uses it. Not safe to run in several threads of one process: the
        phases hand their work to the workers through module-level state, so
        concurrent calls would overwrite each other and return wrong values
        rather than raising.
    """

    validate_settings(mem_limit_mb)

    structure = prepare_structure(pdb_file, pkas_file)
    surface_points = construct_surface_grid(structure, r_probe)
    features = {"ID": structure.name}
    print(f"calculating {structure.name}")
    features = calculate_structure_features(structure, ph, r_probe, features, full_features=full_features)
    features = calculate_surface_grid_features(structure, surface_points, ph, hydro_scale, features, full_features=full_features)
    if full_features:
        features = calculate_average_chargesurface_grid_features(structure, surface_points, ph, features)
    features = calculate_shell_features(structure, surface_points, ph, features, full_features=full_features, mem_limit_mb=mem_limit_mb)

    calculated_features = pd.Series(features).to_frame().transpose()

    try:
        calculated_features = pd.concat([open_output_file(out_file), calculated_features])

    except FileNotFoundError:
        pass

    write_output_file(calculated_features, out_file)


def main():
    parsed_arguments = parse_arguments()
    calculate(*parsed_arguments)


if __name__ == "__main__":
    main()
