"""Writes the viewer files that let a user look at the Prodes surface properties.

The surface grid and its projected electrostatic potential are the most
informative things Prodes computes, and until now they were thrown away as soon
as the summary features had been calculated. These writers turn them into files
that open directly in PyMOL or ChimeraX.

Nothing here recomputes anything. The values written are the same objects the
features were calculated from, so a picture and a feature value can never
disagree about what the potential was.
"""

import numpy as np

# The ramp is always symmetric, so that white means zero. An asymmetric ramp
# would make an uncharged patch read as a charged one. The limit itself is taken
# from the data rather than fixed, so no point is ever clipped and the full
# colour range is used: a fixed 4 V clips glucose oxidase, which reaches -5.02 V.
EP_LIMIT_FLOOR_VOLTS = 1.0

# Projected hydrophobicity spans roughly -0.6 to +0.5, so it needs its own floor
# and an extra decimal. Blue to white to orange rather than red to white to blue,
# because a hydrophobicity map that used the electrostatic colours would be read
# as one at a glance.
MHP_LIMIT_FLOOR = 0.2
MHP_DECIMALS = 2
SPHERE_SCALE = 0.6

EP_PALETTE = "red_white_blue"
CHIMERAX_EP_PALETTE = "red:white:blue"

# Hydrophobicity is shown as two tiers of green over a neutral base rather than a
# ramp. Almost the whole surface of a soluble protein is hydrophilic, so a ramp
# spends most of its range on the part nobody is looking for. Fixed cutoffs also
# make two proteins directly comparable, which an auto-scaled ramp does not: the
# same green means the same hydrophobicity on every structure.
HYDROPHOBIC_CUTOFF = 0.03
STRONGLY_HYDROPHOBIC_CUTOFF = 0.10

EP_DESCRIPTION = "Electrostatic potential in volts. Red is negative, blue positive, white zero."
HYDROPHOBICITY_DESCRIPTION = (
    "Projected hydrophobicity. Green marks hydrophobic surface, pale above "
    f"{HYDROPHOBIC_CUTOFF} and dark above {STRONGLY_HYDROPHOBIC_CUTOFF}; grey is hydrophilic."
)

# The B factor column is six characters wide. A value of -100.00 or 1000.00
# needs seven and silently pushes the element symbol out of its own column,
# producing a line that still looks plausible. Real potentials are single digit
# volts, so anything outside this range is a bug worth surfacing.
B_FACTOR_MIN = -99.99
B_FACTOR_MAX = 999.99


def surface_point_arrays(surface_points):
    """Returns (coords, ep, lipo) as arrays, read off the property points.

    Call this before the average charge phase runs. That phase overwrites the
    electrostatic potential on these same objects with one computed from partial
    rather than formal charges, so arrays taken afterwards would hold a
    different physical quantity under the same name.
    """

    coords = np.array([[point.x, point.y, point.z] for point in surface_points])
    ep = np.array([point.ep for point in surface_points])
    lipo = np.array([point.lipo for point in surface_points])

    return coords, ep, lipo


def ep_colouring(name, limit):
    """Returns the PyMOL commands that colour points by electrostatic potential."""

    return f"spectrum b, {EP_PALETTE}, {name}, {-limit}, {limit}"


def hydrophobicity_colouring(name):
    """Returns the PyMOL commands that colour points by hydrophobicity, in two tiers.

    The cutoffs are absolute rather than derived from the structure, so the same
    green means the same hydrophobicity on any protein.
    """

    return f"color grey80, {name}\n" f"color palegreen, {name} and b > {HYDROPHOBIC_CUTOFF}\n" f"color forest, {name} and b > {STRONGLY_HYDROPHOBIC_CUTOFF}"


def chimerax_ep_colouring(limit):
    """Returns the ChimeraX command that colours points by electrostatic potential."""

    return f"color bfactor #2 palette {CHIMERAX_EP_PALETTE} range {-limit},{limit}"


def chimerax_hydrophobicity_colouring():
    """Returns the ChimeraX commands for the same two green tiers as PyMOL."""

    return "color #2 gray(200)\n" f"color #2 & @@bfactor>{HYDROPHOBIC_CUTOFF} palegreen\n" f"color #2 & @@bfactor>{STRONGLY_HYDROPHOBIC_CUTOFF} forestgreen"


def charged_residue_underlay(name):
    """Cartoon plus the charged side chains, coloured to match the surface ramp.

    Set up but invisible: the cloud is opaque by default, so this sits behind it
    until the user makes the spheres see-through. Everything is then already in
    place and revealing it is a single command.

    cartoon_color is set explicitly because the ribbon otherwise takes its colour
    from each residue's own atoms, so colouring a side chain would tint the
    backbone running through it.
    """

    return (
        f"set sphere_transparency, 0, {name}\n"
        "show cartoon, protein\n"
        "color grey80, protein\n"
        "set cartoon_color, grey80, protein\n"
        "show sticks, protein and resn ASP+GLU and sidechain\n"
        "show sticks, protein and resn LYS+ARG+HIS and sidechain\n"
        "color red, protein and resn ASP+GLU and sidechain\n"
        "color blue, protein and resn LYS+ARG+HIS and sidechain"
    )


def hydrophobic_residues(hydro_scale):
    """Returns the residue codes the given scale calls hydrophobic, most first.

    Taken from the scale the run actually used rather than a fixed textbook
    list, so the side chains shown always match the map above them even when a
    different scale is selected.
    """

    from prodes.data import hydrophobic_scale

    scale = hydrophobic_scale(hydro_scale)

    return [residue for residue, value in sorted(scale.items(), key=lambda item: -item[1]) if value > 0]


def hydrophobic_residue_underlay(name, hydro_scale):
    """Cartoon plus the hydrophobic side chains, in the same green as the surface.

    Green rather than the electrostatic red and blue: those mean charge, and
    putting them under a hydrophobicity map invites reading one property as the
    other.
    """

    residues = "+".join(hydrophobic_residues(hydro_scale))

    return (
        f"set sphere_transparency, 0, {name}\n"
        "show cartoon, protein\n"
        "color grey80, protein\n"
        "set cartoon_color, grey80, protein\n"
        f"show sticks, protein and resn {residues} and sidechain\n"
        f"color forest, protein and resn {residues} and sidechain"
    )


def symmetric_limit(values, floor=EP_LIMIT_FLOOR_VOLTS, decimals=1):
    """Returns the ramp limit to use for these values: symmetric, rounded, never clipping.

    Symmetric keeps white at zero. Taking the limit from the data means the full
    colour range is used and nothing is clipped, at the cost of two proteins
    getting different ramps. To compare proteins directly, pass the same limit
    for both; the value used is written into the generated scripts.
    """

    extreme = float(np.abs(values).max()) if len(values) else 0.0
    step = 10.0**-decimals

    return max(round(extreme + step / 2 - step / 100, decimals), floor)


def point_pdb_line(index, x, y, z, value):
    """Formats one surface point as a PDB pseudo atom carrying value in the B factor column.

    Serial and residue numbers wrap rather than overflow their columns. A large
    structure has more surface points than PDB numbering allows, and for a point
    cloud drawn as dots the numbering carries no meaning anyway.

    The residue is named EPT rather than HOH on purpose: both PyMOL and ChimeraX
    treat waters specially and several of their default presets hide them, which
    gives an empty screen that looks like a failure.
    """

    if not B_FACTOR_MIN <= value <= B_FACTOR_MAX:
        raise ValueError(f"value {value} does not fit the PDB B factor column; expected {B_FACTOR_MIN} to {B_FACTOR_MAX}")

    serial = index % 99999 + 1
    resnum = index % 9999 + 1

    return f"HETATM{serial:5d}  C   EPT X{resnum:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{value:6.2f}           C\n"


def write_point_pdb(coords, values, out_file):
    """Writes the surface points as pseudo atoms, one per point, value in the B factors."""

    with open(out_file, "w") as handle:
        for index, ((x, y, z), value) in enumerate(zip(coords, values, strict=True)):
            handle.write(point_pdb_line(index, x, y, z, value))
        handle.write("END\n")


def pymol_script(structure_file, point_file, colouring, underlay, name="surface_ep", description=EP_DESCRIPTION, sphere_scale=SPHERE_SCALE):
    """Returns a PyMOL script that loads the structure and its coloured surface points.

    Paths are written relative to the script, so the whole directory can be
    zipped and opened on another machine. PyMOL resolves a relative path against
    its own working directory rather than the script's, so the script has to be
    opened from inside the directory that holds it. The header says so.
    """

    return f"""# Prodes surface property map.
#
# Open from inside this directory, because the paths below are relative to it:
#
#     cd <this directory>
#     pymol {point_file.rsplit("/", 1)[-1].replace(".pdb", ".pml")}
#
# {description}

load ./{structure_file}, protein
load ./{point_file}, {name}

hide everything

# Opaque spheres. Transparency looks appealing on a handful of points and washes
# the whole map out on a real one: a protein surface is tens of thousands of
# overlapping spheres, and each transparent layer blends further toward the
# background until the cloud is nearly white.
show spheres, {name}
set sphere_scale, {sphere_scale}, {name}
{colouring}

# The structure underneath, already set up but hidden behind the opaque cloud.
{underlay}

bg_color white
set ray_opaque_background, 1
orient {name}

# To see the structure underneath, make the cloud see-through:
#
#     set sphere_transparency, 0.4, {name}
#
# Ray tracing this many spheres takes minutes. Frame the view first, then ray.
"""


def chimerax_script(structure_file, point_file, colouring, sphere_scale=SPHERE_SCALE):
    """Returns a ChimeraX script equivalent to the PyMOL one.

    ChimeraX sizes spheres by an absolute radius rather than a scale factor, so
    the scale is converted using the 1.7 Angstrom van der Waals radius of the
    carbon pseudo atoms, keeping both viewers visually consistent.
    """

    return f"""# Prodes surface property map.
# Open from inside this directory: the paths below are relative to it.

open {structure_file} name protein
open {point_file} name surface_ep

# The structure is hidden: the point cloud already describes the surface, and a
# translucent surface underneath only washes the colours out.
hide #1 target ams

style #2 sphere
size #2 atomRadius {sphere_scale * 1.7:.2f}
{colouring}

set bgColor white
view

# To read a patch back to the residues underneath:
#
# transparency #2 60
# show #1 target a
# style #1 stick
"""
