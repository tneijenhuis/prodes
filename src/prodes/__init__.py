from prodes import run
from prodes.io import parser
from prodes.run import calculate


def run_prodes(
    pdb_file,
    out_file,
    pkas_file=None,
    ph=7,
    r_probe=1.4,
    hydro_scale="mj_scaled",
    full_features=False,
    mem_limit_mb=None,
    ionic_strength_molar=run.DEFAULT_IONIC_STRENGTH_MOLAR,
):
    """Runs Prodes similar as the commandline tool.

    A thin wrapper around prodes.run.calculate; every option the commandline
    accepts is forwarded, so the two entry points cannot diverge.

    Arguments
    required:
        pdb_file: path to the pdb file to be analysed
        out_file: path to the output csv file which should be written
    Optional:
        pkas_file: an output file of PROPKA which is required for custom pKa assignment
        ph: the ph at which protonation states should be calculated
        r_probe: the radius of the probe used to calculate the solvent accessible surface area
        hydro_scale: the abbreviation of the hydrophobicity scale used (scales can be found in data/hydrophobicity)
        full_features: when True, calculates the full legacy 105-feature set including
            redundant features. Defaults to False for the reduced, non-redundant set.
        ionic_strength_molar: ionic strength of the buffer in mol/L, which screens
            the projected electrostatic potential. 0 disables screening and
            gives the unscreened potential of versions before 5.0.
        mem_limit_mb: memory budget in MB for the intermediate NumPy arrays of this
            run, divided among the workers. Pass explicitly in worker processes
            rather than relying on the PRODES_MEM_LIMIT_MB environment variable.
    """

    calculate(
        pdb_file,
        out_file,
        pkas_file=pkas_file,
        ph=ph,
        r_probe=r_probe,
        hydro_scale=hydro_scale,
        full_features=full_features,
        mem_limit_mb=mem_limit_mb,
        ionic_strength_molar=ionic_strength_molar,
    )


def read(file):
    """Loads files
    List of supported extentions:
    .pdb
    .pdb.zip
    .pka: a pKa file in the json layout written by prodes.io.pka_converter,
        rather than the raw output of PROPKA or pypka"""

    if file.endswith((".pdb", ".pdb.zip")):
        structure = parser.PDBparser().parse(file)
        return structure

    elif file[-4:] == ".pka":
        pkas = parser.read_pka(file)
        return pkas

    else:
        raise ValueError("File extention not recognized by prodes")
