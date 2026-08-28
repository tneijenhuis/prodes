"""Tests the output bundle and the viewer files inside it.

A run now produces one bundle rather than a row appended to a growing CSV. What
these tests protect is the property that makes the bundle worth having: the
viewer files and the feature values come out of the same calculation, so a
picture and a number can never disagree about what the potential was.
"""

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from prodes.output import read_features, read_surface_points, run_metadata, write_bundle
from prodes.run import calculate
from prodes.viz import HYDROPHOBIC_CUTOFF, STRONGLY_HYDROPHOBIC_CUTOFF, ep_colouring, point_pdb_line, pymol_script, symmetric_limit, write_point_pdb

REPO_ROOT = Path(__file__).resolve().parent.parent
SMALL_STRUCTURE = REPO_ROOT / "tests" / "data" / "ARH96693.pdb.zip"
HYDROGEN_PAIR = (REPO_ROOT / "tests" / "data" / "1GDW.pdb.zip", REPO_ROOT / "tests" / "data" / "1GDW_h.pdb.zip")

EXPECTED_MEMBERS = {
    "ARH96693.pdb",
    "ARH96693_ep.cxc",
    "ARH96693_ep.pdb",
    "ARH96693_ep.pml",
    "ARH96693_features.csv",
    "ARH96693_hydrophobicity.cxc",
    "ARH96693_hydrophobicity.pdb",
    "ARH96693_hydrophobicity.pml",
    "ARH96693_surface_points.csv",
    "README.txt",
    "prodes_run.json",
}


@pytest.fixture(scope="module")
def bundle_zip(tmp_path_factory):
    """One zip bundle of the small structure, built once for the module."""

    out_file = tmp_path_factory.mktemp("zip") / "ARH96693.zip"
    calculate(str(SMALL_STRUCTURE), str(out_file))
    return out_file


@pytest.fixture(scope="module")
def bundle_dir(bundle_zip, tmp_path_factory):
    """The bundle unpacked, which is how a user reaches the viewer files."""

    directory = tmp_path_factory.mktemp("unpacked")
    with zipfile.ZipFile(bundle_zip) as archive:
        archive.extractall(directory)
    return directory / "ARH96693"


def test_the_zip_holds_every_member_under_one_top_level_directory(bundle_zip):
    """Unpacking the zip gives one directory, so it cannot scatter files."""

    with zipfile.ZipFile(bundle_zip) as archive:
        names = archive.namelist()

    assert {name.split("/", 1)[1] for name in names} == EXPECTED_MEMBERS
    assert {name.split("/", 1)[0] for name in names} == {"ARH96693"}


def test_a_path_that_is_not_a_zip_is_refused(tmp_path):
    """One output shape only, so nothing downstream has to branch on which it got."""

    with pytest.raises(ValueError, match=r"must be a \.zip path"):
        calculate(str(SMALL_STRUCTURE), str(tmp_path / "bundle"))


def test_the_tables_read_back_out_of_the_archive(bundle_zip):
    """read_features and read_surface_points are the counterparts to write_bundle."""

    features, points = read_features(bundle_zip), read_surface_points(bundle_zip)

    assert len(features) == 1
    assert features["ID"].iloc[0] == "ARH96693"
    assert list(points.columns) == ["x", "y", "z", "ep_volts", "hydrophobicity"]
    assert len(points) > 0


def test_reading_a_member_that_is_not_there_says_so(tmp_path):
    """A clear error beats a StopIteration from somewhere inside the reader."""

    empty = tmp_path / "empty.zip"
    with zipfile.ZipFile(empty, "w") as archive:
        archive.writestr("note.txt", "nothing here")

    with pytest.raises(ValueError, match="_features.csv"):
        read_features(empty)


def test_the_run_record_says_what_was_run(bundle_dir):
    """A bundle carries its own provenance, so a stray output can be identified."""

    record = json.loads((bundle_dir / "prodes_run.json").read_text())

    assert record["settings"] == {
        "ph": 7,
        "r_probe": 1.4,
        "hydro_scale": "mj_scaled",
        "full_features": False,
        "ionic_strength_molar": 0.15,
    }
    # ARH96693 is an AlphaFold model with three cysteines and no disulfide, so
    # the key has to be present and zero rather than absent.
    assert record["disulfides"] == 0
    assert record["surface_points"] > 0
    assert record["ep_min_volts"] <= record["ep_max_volts"]
    assert record["input_file"].endswith("ARH96693.pdb.zip")


def test_the_point_cloud_and_the_point_table_describe_the_same_points(bundle_dir):
    """The PDB a viewer opens and the CSV a script reads must not diverge."""

    points = pd.read_csv(bundle_dir / "ARH96693_surface_points.csv")
    lines = [line for line in (bundle_dir / "ARH96693_ep.pdb").read_text().splitlines() if line.startswith("HETATM")]

    assert len(lines) == len(points)
    assert float(lines[0][30:38]) == pytest.approx(points["x"].iloc[0], abs=5e-4)
    assert float(lines[0][60:66]) == pytest.approx(points["ep_volts"].iloc[0], abs=5e-3)


def test_the_structure_travels_with_the_bundle(bundle_dir):
    """A zipped input still yields a plain PDB the viewer scripts can load."""

    structure = bundle_dir / "ARH96693.pdb"

    assert structure.exists()
    assert any(line.startswith("ATOM") for line in structure.read_text().splitlines())


def test_the_hydrophobicity_view_uses_the_same_points_and_its_own_colours(bundle_dir):
    """Both maps describe the same surface, so only the B factor and the ramp differ.

    Sharing the electrostatic red and blue would make a hydrophobicity figure read
    as an electrostatic one at a glance.
    """

    points = pd.read_csv(bundle_dir / "ARH96693_surface_points.csv")
    ep_lines = [line for line in (bundle_dir / "ARH96693_ep.pdb").read_text().splitlines() if line.startswith("HETATM")]
    mhp_lines = [line for line in (bundle_dir / "ARH96693_hydrophobicity.pdb").read_text().splitlines() if line.startswith("HETATM")]

    assert len(mhp_lines) == len(ep_lines) == len(points)
    # Same coordinates, different value in the B factor column.
    assert [line[30:54] for line in mhp_lines] == [line[30:54] for line in ep_lines]
    assert float(mhp_lines[0][60:66]) == pytest.approx(points["hydrophobicity"].iloc[0], abs=5e-3)

    script = (bundle_dir / "ARH96693_hydrophobicity.pml").read_text()
    assert "forest" in script and "palegreen" in script
    assert "red_white_blue" not in script, "the electrostatic colours would make this read as a charge map"


def test_the_hydrophobicity_tiers_are_fixed_and_ordered(bundle_dir):
    """Absolute cutoffs are what make the same green mean the same thing on any protein.

    They are applied in increasing order, so the darker tier paints over the
    paler one rather than being overwritten by it.
    """

    lines = [line for line in (bundle_dir / "ARH96693_hydrophobicity.pml").read_text().splitlines() if line.startswith("color ")]

    assert lines[0].startswith("color grey80"), "the neutral base must be painted first"
    assert f"b > {HYDROPHOBIC_CUTOFF}" in lines[1]
    assert f"b > {STRONGLY_HYDROPHOBIC_CUTOFF}" in lines[2]
    assert HYDROPHOBIC_CUTOFF < STRONGLY_HYDROPHOBIC_CUTOFF


def test_the_pymol_script_uses_relative_paths_that_resolve_beside_it(bundle_dir):
    """The bundle must open after being copied to another machine."""

    script = (bundle_dir / "ARH96693_ep.pml").read_text()
    targets = [line.split(",")[0].split()[1] for line in script.splitlines() if line.startswith("load ")]

    assert targets, "the script loads nothing"
    for target in targets:
        assert not Path(target).is_absolute()
        assert "\\" not in target, "a Windows separator would break the bundle on Linux and vice versa"
        assert (bundle_dir / target).exists()
        # Pins the documented requirement to open the script from inside the
        # bundle: PyMOL resolves these against its own working directory.
        assert not (bundle_dir.parent / target).exists()


def test_the_ramp_is_symmetric_and_covers_the_data(bundle_dir):
    """White must mean zero, and no point may be clipped."""

    import re

    script = (bundle_dir / "ARH96693_ep.pml").read_text()
    low, high = (float(value) for value in re.search(r"spectrum b, red_white_blue, surface_ep, (\S+), (\S+)", script).groups())
    ep = pd.read_csv(bundle_dir / "ARH96693_surface_points.csv")["ep_volts"]

    assert low == -high, "an asymmetric ramp would stop white meaning zero"
    assert low <= ep.min() and ep.max() <= high, "the ramp clips a real value"
    assert "set sphere_scale, 0.6, surface_ep" in script


def test_the_structure_is_set_up_underneath_but_hidden(bundle_dir):
    """The cartoon and charged side chains are preloaded behind the opaque cloud.

    Revealing them is then one command rather than eight, which is the whole
    point of putting them in the script.
    """

    script = (bundle_dir / "ARH96693_ep.pml").read_text()

    assert "show cartoon, protein" in script
    assert "set cartoon_color, grey80, protein" in script, "the ribbon would otherwise be tinted by the side chain colours"
    assert "color red, protein and resn ASP+GLU and sidechain" in script
    assert "color blue, protein and resn LYS+ARG+HIS and sidechain" in script
    assert "set sphere_transparency, 0.4" in script, "the reveal command must be documented in the file"


def test_the_hydrophobicity_view_shows_hydrophobic_residues_in_green(bundle_dir):
    """Red and blue mean charge; under a hydrophobicity map they invite misreading."""

    from prodes.viz import hydrophobic_residues

    script = (bundle_dir / "ARH96693_hydrophobicity.pml").read_text()
    expected = "+".join(hydrophobic_residues("mj_scaled"))

    assert "show cartoon, protein" in script
    assert f"color forest, protein and resn {expected} and sidechain" in script
    assert "resn ASP+GLU" not in script
    assert "resn LYS+ARG+HIS" not in script


def test_the_residues_shown_come_from_the_scale_in_use():
    """A fixed textbook list would disagree with the map whenever --hydro changed."""

    from prodes.data import hydrophobic_scale
    from prodes.viz import hydrophobic_residues

    residues = hydrophobic_residues("mj_scaled")
    scale = hydrophobic_scale("mj_scaled")

    assert residues, "no residue scores as hydrophobic, which cannot be right"
    assert all(scale[residue] > 0 for residue in residues)
    assert all(scale[residue] <= 0 for residue in scale if residue not in residues)
    # Most hydrophobic first, so the selection reads in a sensible order.
    assert residues == sorted(residues, key=lambda residue: -scale[residue])


def test_each_script_is_titled_for_the_property_it_shows(bundle_dir):
    """The header must not claim to be the potential map on the hydrophobicity one."""

    phobic = (bundle_dir / "ARH96693_hydrophobicity.pml").read_text()

    assert "electrostatic potential" not in phobic.split("load ")[0].lower()
    assert "hydrophobic" in phobic.split("load ")[0].lower()


def test_the_spheres_are_opaque_by_default(bundle_dir):
    """Transparency washes a real surface out to near white.

    A protein surface is tens of thousands of overlapping spheres, so every
    transparent layer blends further toward the background. The setting is
    offered as a commented line for reading patches back to residues, but it
    must not be active by default.
    """

    active = [line for line in (bundle_dir / "ARH96693_ep.pml").read_text().splitlines() if not line.lstrip().startswith("#")]

    assert not any("sphere_transparency" in line and "0.4" in line for line in active)
    assert any("set sphere_transparency, 0, " in line for line in active), "opaque must be set explicitly, so it is easy to find and change"
    assert any("show spheres, surface_ep" in line for line in active)


def test_the_potential_does_not_change_meaning_under_full_features(tmp_path):
    """The average charge phase overwrites point.ep with a partial charge potential.

    The bundle is written before that phase runs. If it ever moves after it, the
    ep column would silently hold a different physical quantity whenever
    --full-features was passed, with no error and no column rename.
    """

    reduced = tmp_path / "reduced.zip"
    full = tmp_path / "full.zip"
    calculate(str(SMALL_STRUCTURE), str(reduced))
    calculate(str(SMALL_STRUCTURE), str(full), full_features=True)

    pd.testing.assert_series_equal(read_surface_points(reduced)["ep_volts"], read_surface_points(full)["ep_volts"])


def test_adding_hydrogens_does_not_change_the_map(tmp_path):
    """Prodes ignores hydrogens, so the same structure with and without them must agree."""

    plain, protonated = HYDROGEN_PAIR
    calculate(str(plain), str(tmp_path / "plain.zip"))
    calculate(str(protonated), str(tmp_path / "protonated.zip"))

    without = read_surface_points(tmp_path / "plain.zip")
    with_hydrogens = read_surface_points(tmp_path / "protonated.zip")

    assert len(without) == len(with_hydrogens)
    np.testing.assert_allclose(without["ep_volts"], with_hydrogens["ep_volts"])


def test_point_lines_keep_every_pdb_column_in_its_own_place():
    """Asserts field positions rather than line length.

    A B factor that overflows its column still produces a line of plausible
    length while pushing the element symbol out of columns 77 and 78, so a
    length check passes on a corrupt line.
    """

    line = point_pdb_line(0, -12.345, 6.789, 0.0, -3.21)

    assert line[0:6] == "HETATM"
    assert line[17:20] == "EPT", "PyMOL and ChimeraX hide waters by default, so this must not be HOH"
    assert float(line[30:38]) == pytest.approx(-12.345)
    assert float(line[38:46]) == pytest.approx(6.789)
    assert float(line[54:60]) == pytest.approx(1.00)
    assert float(line[60:66]) == pytest.approx(-3.21)
    assert line[76:78].strip() == "C"
    assert len(line.rstrip("\n")) == 78


@pytest.mark.parametrize(
    ("index", "serial", "resnum"),
    [(0, 1, 1), (9998, 9999, 9999), (9999, 10000, 1), (99998, 99999, 9), (99999, 1, 10)],
)
def test_serial_and_residue_numbers_wrap_rather_than_overflow(index, serial, resnum):
    """Catalase has more surface points than PDB numbering allows for either field."""

    line = point_pdb_line(index, 0.0, 0.0, 0.0, 0.0)

    assert int(line[6:11]) == serial
    assert int(line[22:26]) == resnum


def test_a_potential_too_large_for_the_b_factor_column_is_refused():
    """Silently corrupting the line would be worse than failing."""

    with pytest.raises(ValueError, match="B factor"):
        point_pdb_line(0, 0.0, 0.0, 0.0, 1000.0)


def test_write_point_pdb_ends_the_file_properly(tmp_path):
    """Some viewers need the END record to stop reading."""

    out_file = tmp_path / "points.pdb"
    write_point_pdb(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]), np.array([1.5, -1.5]), out_file)
    lines = out_file.read_text().splitlines()

    assert len(lines) == 3
    assert lines[-1] == "END"


@pytest.mark.parametrize(
    ("low", "high", "expected"),
    [(-3.36, 1.17, 3.4), (-5.02, -0.76, 5.1), (-0.2, 0.3, 1.0), (0.0, 0.0, 1.0)],
)
def test_the_ramp_limit_rounds_outward_and_has_a_floor(low, high, expected):
    """Rounding must never round down onto a real value, and a flat map needs a usable range."""

    assert symmetric_limit(np.array([low, high])) == pytest.approx(expected)


def test_the_pymol_script_tells_the_user_how_to_open_it():
    """The cd is a real requirement, not a nicety, so it belongs in the file."""

    script = pymol_script("thing.pdb", "thing_ep.pdb", ep_colouring("surface_ep", 4.0), "")

    assert "cd " in script
    assert "pymol thing_ep.pml" in script


def test_write_bundle_returns_the_path_it_wrote(tmp_path):
    """So a caller can chain on the result rather than rebuilding the path."""

    coords = np.array([[0.0, 0.0, 0.0]])
    features = pd.DataFrame([{"ID": "x", "Area": 1.0}])
    written = write_bundle(
        tmp_path / "out.zip", "x", features, coords, np.array([0.5]), np.array([0.1]), SMALL_STRUCTURE, {"prodes_version": "test"}, "mj_scaled"
    )

    assert written == tmp_path / "out.zip"
    assert zipfile.is_zipfile(written)


def test_the_run_record_carries_the_disulfide_count():
    """The count is how a user checks that detection saw what they expected.

    It is metadata rather than a feature column: the feature set is pinned by the
    two YAML files and asserted column by column elsewhere, and a structural
    property of the input does not belong among descriptors of its surface.
    """

    record = run_metadata("structure.pdb", {}, 10, np.array([0.0, 1.0]), disulfides=17)

    assert record["disulfides"] == 17
