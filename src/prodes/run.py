import argparse
import os

import numpy as np
import pandas as pd
from dotenv import load_dotenv

from prodes.calculations.sasa import shrake_rupley, shape
from prodes.io import parser as ps
from prodes.calculations import grid_wizard
from prodes.calculations.standard_equations import trimean
from prodes.calculations import geometry

load_dotenv()


def parse_arguments():
    """Parses all arguments given in the commandline"""

    env_full = os.getenv("PRODES_FULL_FEATURES", "false").lower() in ("true", "1", "yes")

    parser = argparse.ArgumentParser(description='Calculate descriptors from atomic data')
    parser.add_argument("pdb_file", help="file location of a pdb or pqr file", type=str)
    parser.add_argument("out_file", help="file path of the output csv file", type=str)
    parser.add_argument("-p", "--pka", help="file location of the pka propka output", type=str, default=None)
    parser.add_argument("--probe", help="Radius of the surface probe", type=float, default=1.4)
    parser.add_argument("--ph", help="pH of the system", type=float, default=7)
    parser.add_argument("--hydro", help="Abriviation of the hydrophobicity scale to be used", type=str, default="mj_scaled")
    parser.add_argument("--full-features", action=argparse.BooleanOptionalAction, default=env_full,
                        help="Calculate the full legacy feature set including redundant features (default: from PRODES_FULL_FEATURES env var)")

    arg = parser.parse_args()

    pdb_file = arg.pdb_file
    out_file = arg.out_file
    pkas_file = arg.pka
    ph = arg.ph
    r_probe = arg.probe
    hydro_scale = arg.hydro
    full_features = arg.full_features

    return  pdb_file, out_file, pkas_file, ph, r_probe, hydro_scale, full_features

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

def calculate_surface_grid_features(structure, surface_points, ph, hydro_scale, features:dict, full_features=True):
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

    charged_atoms = np.array([atom for atom in structure.atoms if atom.charge(ph=ph) != 0])
    for cell in grid.cells.flatten():
        cell_surface_points = cell.filtered_content("Property_point")
        enviroment = np.array([atoms for atoms in grid.grid_content("Atom", cells = grid.find_surrounding_cells(cell))])
        for point in cell_surface_points:
            point.set_ep(charged_atoms, ph)
            point.set_lipo(enviroment, 10, hydro_scale)

    # Electrostatic potential fratures
    eps = np.array([point.ep for point in surface_points])

    features["SurfEpMaxFormal"] = round(eps.max(), 3)
    features["SurfEpMinFormal"] = round(eps.min(), 3)
    features.update({f"{v}Formal": k for v, k in standard_features(eps, "SurfEp", reduced=reduced).items()})

    positive_eps = np.array([ep for ep in eps if ep > 0])
    features["NSurfPosEpFormal"] = len(positive_eps)
    features.update({f"{v}Formal": k for v, k in standard_features(positive_eps, "SurfEpPos", reduced=reduced).items()})

    negative_eps = np.array([ep for ep in eps if ep < 0])
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
    features["NSurfNegMhp"] = len(negative_lipos)
    features.update(standard_features(negative_lipos, "SurfNegMhp", reduced=reduced))

    return features

def calculate_average_chargesurface_grid_features(structure, surface_points, ph, features:dict):

    concat = np.concatenate([structure.heavy_atoms, surface_points])
    grid = grid_wizard.Grid(12)
    grid.construct_cells(concat)
    grid.fill_cells(concat)
    charged_atoms = np.array([atom for atom in structure.atoms if atom.charge(ph=ph, formal=False) != 0])
    for cell in grid.cells.flatten():
        cell_surface_points = cell.filtered_content("Property_point")
        for point in cell_surface_points:
            point.set_ep(charged_atoms, ph, formal=False)

# Electrostatic potential fratures
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


def calculate_shell_features(structure, surface_points, ph:float, features:dict, numb_of_planes=120, full_features=True ):
    """Constructs a number of planes onto which charges are mapped"""

    reduced = not full_features

    surface_coords = np.array([[p.x, p.y, p.z] for p in surface_points])
    distributed_points = geometry.Sunflower_sphere(*geometry.make_vector(structure), 1, numb_of_planes).points

    charged_atoms = [atom for atom in structure.heavy_atoms if atom.charge != 0]
    charged_coords = np.array([[a.x, a.y, a.z] for a in charged_atoms])
    charged_charges = np.array([a.charge(ph) for a in charged_atoms])

    shell_potentials = []
    for point in distributed_points:
        distance = geometry.required_distance(point, structure, surface_coords)
        geometry.move_point(point, structure, distance)

        plane = geometry.find_plane(point, structure)
        projected_coords = geometry.project_point_batch(*plane, charged_coords)
        exits, has_exit = geometry.find_exit_batch(charged_coords, projected_coords, surface_coords)
        potentials = geometry.map_ep_to_plane_batch(
            charged_coords, projected_coords, exits, has_exit, charged_charges
        )
        shell_potentials.append(float(potentials.sum()))

    shell_potentials = np.array(shell_potentials)

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

def calculate_structure_features(structure, ph, r_probe, features:dict, full_features=True):
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


def calculate(pdb_file, out_file, pkas_file=None, ph=7, r_probe=1.4, hydro_scale="mj_scaled", full_features=False):
    """Calculates a list of supported features and returns a csv file

    Arguments
    required:
        pdb_file: path to the pdb file to be analysed
        out_file: path to the output csv file which should be written
    Optional:
        pka_file: a output file of PROPKA which is required for costum pKa assignment
        ph: the ph at which protonation states should be calculated
        r_probe: the radius of the probe used to calculate the solvent accessible surface area
        hydro_scale: the abriviation of the hydrophibicity scale used (scales can be found in data/hydrophobicity)
        full_features: when True, calculates the full legacy 105-feature set including
            redundant features. Defaults to False for the reduced, non-redundant set.
    """

    structure = prepare_structure(pdb_file, pkas_file)
    surface_points = construct_surface_grid(structure, r_probe)
    features = {"ID": structure.name}
    print(f"calculating {structure.name}")
    features = calculate_structure_features(structure, ph, r_probe, features, full_features=full_features)
    features = calculate_surface_grid_features(structure, surface_points, ph, hydro_scale, features, full_features=full_features)
    if full_features:
        features = calculate_average_chargesurface_grid_features(structure, surface_points, ph, features)
    features = calculate_shell_features(structure, surface_points, ph, features, full_features=full_features)

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
