import argparse
import numpy as np
import pandas as pd

from prodes.calculations.sasa import shrake_rupley, shape
from prodes.io import parser as ps
from prodes.calculations import grid_wizard
from prodes.calculations.standard_equations import trimean
from prodes.calculations import geometry


def parse_arguments():
    """Parses all arguments given in the commandline"""

    parser = argparse.ArgumentParser(description='Calculate descriptors from atomic data')
    parser.add_argument("pdb_file", help="file location of a pdb or pqr file", type=str)
    parser.add_argument("out_file", help="file path of the output csv file", type=str)
    parser.add_argument("-p", "--pka", help="file location of the pka propka output", type=str, default=None)
    parser.add_argument("--probe", help="Radius of the surface probe", type=float, default=1.4)
    parser.add_argument("--ph", help="pH of the system", type=float, default=7)
    parser.add_argument("--hydro", help="Abriviation of the hydrophobicity scale to be used", type=str, default="mj_scaled")

    arg = parser.parse_args()

    pdb_file = arg.pdb_file
    out_file = arg.out_file
    pkas_file = arg.pka
    ph = arg.ph
    r_probe = arg.probe
    hydro_scale = arg.hydro


    return  pdb_file, out_file, pkas_file, ph, r_probe, hydro_scale 

def open_output_file(out_file):
    """opens the output file"""
    return pd.read_csv(out_file)

def write_output_file(dataframe, out_file):
    """writes a CSV file"""
    dataframe.to_csv(out_file, index=False)


def standard_features(values, name=""):
    """calculates the mean, trimean, median, sum and standard diviation of a np array of values"""
    features = {}
    
    if len(values) != 0:
        features[f"{name}Mean"] = round(values.mean(), 3)
        features[f"{name}Trimean"] = round(trimean(values), 3)
        features[f"{name}Median"] = round(np.median(values), 3)
        features[f"{name}Sum"] = round(np.sum(values), 3)
        features[f"{name}Std"] = round(values.std(), 3)

    else:
        features[f"{name}Mean"] = 0
        features[f"{name}Trimean"] = 0
        features[f"{name}Median"] = 0
        features[f"{name}Sum"] = 0
        features[f"{name}Std"] = 0

    return features

def calculate_surface_grid_features(structure, surface_points, ph, hydro_scale, features:dict):
    """calculates the features from the surface grid"""

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
    features.update({f"{v}Formal": k for v, k in standard_features(eps, "SurfEp").items()})

    positive_eps = np.array([ep for ep in eps if ep > 0])
    features["NSurfPosEpFormal"] = len(positive_eps)
    features.update({f"{v}Formal": k for v, k in standard_features(positive_eps, "SurfEpPos").items()})

    negative_eps = np.array([ep for ep in eps if ep < 0])
    features["NSurfNegEpFormal"] = len(negative_eps)
    features.update({f"{v}Formal": k for v, k in standard_features(negative_eps, "SurfEpNeg").items()})
    

    # Hydrophobic potential features
    lipos = np.array([point.lipo for point in surface_points])

    features["SurfMhpMax"] = round(lipos.max(), 3)
    features["SurfMhpMin"] = round(lipos.min(), 3)
    features.update(standard_features(lipos, "SurfMhp"))

    positive_lipos = np.array([lipo for lipo in lipos if lipo > 0])
    features["NSurfPosMhp"] = len(positive_lipos)
    features.update(standard_features(positive_lipos, "SurfPosMhp"))

    negative_lipos = np.array([lipo for lipo in lipos if lipo < 0])
    features["NSurfNegMhp"] = len(negative_lipos)
    features.update(standard_features(negative_lipos, "SurfNegMhp"))

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


def calculate_shell_features(structure, surface_points, ph:float, features:dict, numb_of_planes=120 ):
    """Constructs a number of planes onto which charges are mapped"""

    distributed_points = geometry.Sunflower_sphere(*geometry.make_vector(structure), 1, numb_of_planes).points
    shell_potentials = []
    surface_grid = grid_wizard.Grid(2)
    surface_grid.construct_cells(surface_points)
    surface_grid.fill_cells(surface_points)
    charged_atoms = np.array([atom for atom in structure.heavy_atoms if atom.charge != 0])
    for i,point in enumerate(distributed_points):
        distance = geometry.required_distance(point, structure, surface_points)
        geometry.move_point(point, structure, distance)  

        plane = geometry.find_plane(point, structure)
        shell_potential = 0
        for atom in charged_atoms:
            
            projected_atom = geometry.project_point(*plane, *geometry.make_vector(atom))
            surface_exit = geometry.find_exit(geometry.make_vector(atom), projected_atom, surface_grid)
            projected_potential = geometry.map_ep_to_plane(atom, projected_atom, surface_exit, ph)
            shell_potential += projected_potential

        shell_potentials.append(shell_potential)
        # print(shell_potential, shell_potentials)
    shell_potentials = np.array(shell_potentials)
    # print(shell_potentials)

    features["ShellEpMaxFormal"] = round(shell_potentials.max(), 3)
    features["ShellEpminFormal"] = round(shell_potentials.min(), 3)
    features.update({f"{v}Formal": k for v, k in standard_features(shell_potentials, "ShellEp").items()})

    pos_potentials = np.array([potential for potential in shell_potentials if potential > 0])
    features["NShellPosEpFormal"] = len(pos_potentials)
    features.update({f"{v}Formal": k for v, k in standard_features(pos_potentials, "ShellEpPos").items()})

    neg_potentials = np.array([potential for potential in shell_potentials if potential < 0])
    features["NShellNegEpFormal"] = len(neg_potentials)
    features.update({f"{v}Formal": k for v, k in standard_features(neg_potentials, "ShellEpNeg").items()})

    return features

def calculate_structure_features(structure, ph, r_probe, features:dict):
    """calculates general structure features, including"""
    
    features["Molecular weight"] = structure.mw
    features["Isoelectric point"] = structure.isoelectric_point()
    features["Dipole"] = structure.dipole(ph)
    features["Formal charge"] = structure.charge(ph)
    features["Average charge"] = structure.charge(ph, formal=False)

    if structure.surface_done is False:
        grid_size = 10 + (r_probe - 1.4) * 2
        grid = grid_wizard.Grid(grid_size)
        grid.construct_cells(structure.heavy_atoms)
        grid.fill_cels(structure.heavy_atoms)

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


def calculate(pdb_file, out_file, pkas_file=None, ph=7, r_probe=1.4, hydro_scale="mj_scaled"):
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
    """    
    
    structure = prepare_structure(pdb_file, pkas_file)
    surface_points = construct_surface_grid(structure, r_probe)
    features = {"ID": structure.name}
    print(f"calculating {structure.name}")
    features = calculate_structure_features(structure, ph, r_probe, features)
    features = calculate_surface_grid_features(structure, surface_points, ph, hydro_scale, features)
    features = calculate_average_chargesurface_grid_features(structure, surface_points, ph, features)
    features = calculate_shell_features(structure, surface_points, ph, features)

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
