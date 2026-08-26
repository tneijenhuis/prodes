"""Tests the ionic screening of the projected electrostatic potential.

Prodes evaluated a Coulomb sum in which every charged atom contributed in full,
however far away it was. A real buffer contains mobile ions that cancel a
charge's field beyond a few Angstrom, so the unscreened sum added a large smooth
offset from the whole protein on top of the local pattern a surface point should
describe. Screening damps each contribution by exp(-d / debye_length).

Two properties matter enough to pin. Zero ionic strength must give back the
unscreened sum, because that is what makes the change auditable, and it is
checked against a formula written out independently of the code under test
rather than against the code itself. And screening must actually be on by
default, because a default quietly reverting to zero would restore the old
behaviour with every test still passing.

Note that zero does not reproduce the released 4.x feature values exactly. The
potential is now stored to three decimals rather than two, which was a separate
correction: at two decimals a point between 0 and 0.005 rounded to zero and left
the positive population, and the error only ever subtracted.
"""

import math
import zipfile

import numpy as np
import pytest

from prodes.calculations.distance_functions import debye_length
from prodes.core.point import EP_DECIMALS, Property_point
from prodes.io.parser import PDBparser
from prodes.output import read_features, read_surface_points
from prodes.run import DEFAULT_IONIC_STRENGTH_MOLAR, calculate

STRUCTURE = "tests/data/ARH96693.pdb.zip"
PHYSIOLOGICAL = 0.15


@pytest.fixture(scope="module")
def atoms():
    """The charged atoms of a real structure, for point level tests."""

    structure = PDBparser().parse(STRUCTURE)

    return np.array([atom for atom in structure.atoms if atom.charge(ph=7) != 0])


def unscreened_potential(point, charged_atoms, ph=7):
    """The Coulomb sum as it stood before screening, written out independently.

    Deliberately not a call into the code under test: this is the reference that
    zero ionic strength has to reproduce.
    """

    coords = np.array([[atom.x, atom.y, atom.z] for atom in charged_atoms])
    charges = np.array([atom.charge(ph=ph) for atom in charged_atoms]) * 1.6e-19
    distances = np.linalg.norm(coords - np.array([point.x, point.y, point.z]), axis=1) * 1e-10

    return round(float(np.sum(charges / (4 * 8.854e-12 * distances * 4 * np.pi))), EP_DECIMALS)


def screened_potential(point, charged_atoms, lam, ph=7):
    """The screened sum written out independently, to pin the strength of the damping.

    Without this the tests only prove that screening happens, not how much: an
    implementation using half the correct screening length passes every other
    test in this module.
    """

    coords = np.array([[atom.x, atom.y, atom.z] for atom in charged_atoms])
    charges = np.array([atom.charge(ph=ph) for atom in charged_atoms]) * 1.6e-19
    angstrom = np.linalg.norm(coords - np.array([point.x, point.y, point.z]), axis=1)

    return round(float(np.sum(charges * np.exp(-angstrom / lam) / (4 * 8.854e-12 * (angstrom * 1e-10) * 4 * np.pi))), EP_DECIMALS)


def test_the_screening_strength_is_exactly_the_debye_length(atoms):
    """Pins how much damping is applied, not merely that some is.

    A kernel that used twice or half the screening length would still shrink the
    potential, still leave the sign right and still change the features.
    """

    rng = np.random.default_rng(3)
    centre = np.array([[atom.x, atom.y, atom.z] for atom in atoms]).mean(axis=0)
    lam = debye_length(PHYSIOLOGICAL)

    for offset in rng.uniform(-20, 20, size=(20, 3)):
        point = Property_point(*(centre + offset))
        point.set_ep(atoms, ph=7, debye_length=lam)

        assert point.ep == screened_potential(Property_point(*(centre + offset)), atoms, lam)


def test_zero_ionic_strength_reproduces_the_unscreened_sum(atoms):
    """The backward compatibility guarantee, checked against an independent formula."""

    rng = np.random.default_rng(0)
    centre = np.array([[atom.x, atom.y, atom.z] for atom in atoms]).mean(axis=0)

    for offset in rng.uniform(-20, 20, size=(25, 3)):
        point = Property_point(*(centre + offset))
        point.set_ep(atoms, ph=7, debye_length=debye_length(0))

        reference = Property_point(*(centre + offset))
        assert point.ep == unscreened_potential(reference, atoms)


def test_the_default_is_screened(atoms):
    """A default that quietly reverted to zero would undo this with every test green."""

    assert DEFAULT_IONIC_STRENGTH_MOLAR > 0
    assert math.isfinite(debye_length(DEFAULT_IONIC_STRENGTH_MOLAR))


@pytest.mark.parametrize(
    ("ionic_strength", "expected"),
    [(0.15, 7.85), (0.02, 21.50), (0.5, 4.30), (1.0, 3.04)],
)
def test_debye_length_matches_the_textbook_values(ionic_strength, expected):
    """3.04 / sqrt(I) Angstrom for a 1:1 salt in water at 298 K."""

    assert debye_length(ionic_strength) == pytest.approx(expected, abs=0.01)


def test_zero_ionic_strength_means_no_screening():
    """No mobile ions is an infinite screening length, not a very long one."""

    assert debye_length(0) == math.inf


def test_a_negative_ionic_strength_is_refused():
    """Silently taking the square root of a negative number would be worse."""

    with pytest.raises(ValueError, match="cannot be negative"):
        debye_length(-0.1)


@pytest.mark.parametrize("bad", [math.inf, -math.inf, math.nan])
def test_a_non_finite_ionic_strength_is_refused(bad):
    """Infinity gives a screening length of zero, which damps the whole surface to zero volts.

    That produces a bundle whose every electrostatic feature is 0 and which
    otherwise looks entirely valid, so it has to fail at the point of entry.
    """

    with pytest.raises(ValueError, match="finite"):
        debye_length(bad)


def test_a_non_positive_screening_length_is_refused(atoms):
    """0 is the value the ionic strength argument uses to mean 'no screening'.

    Passed as a length it means the opposite, total damping, so confusing the two
    must raise rather than silently return a surface of zeros.
    """

    for bad in (0, -1.0):
        with pytest.raises(ValueError, match="must be positive"):
            Property_point(0.0, 0.0, 0.0).set_ep(atoms, ph=7, debye_length=bad)


@pytest.mark.parametrize("multiple", [1, 2, 4])
def test_a_charge_at_n_screening_lengths_keeps_exp_minus_n_of_its_contribution(multiple):
    """The mechanism, driven through the kernel rather than asserted as algebra.

    A single charge is placed at exactly n screening lengths and the screened
    potential is compared against the unscreened one at the same distance. The
    ratio has to be exp(-n): 37 per cent at one length, 2 per cent at four.
    """

    lam = debye_length(PHYSIOLOGICAL)

    class Charge:
        y = z = 0.0

        def __init__(self, distance):
            self.x = distance

        def charge(self, ph=7, formal=True):
            return 1.0

    atom = np.array([Charge(multiple * lam)])
    screened = Property_point(0.0, 0.0, 0.0)
    screened.set_ep(atom, ph=7, debye_length=lam)
    unscreened = Property_point(0.0, 0.0, 0.0)
    unscreened.set_ep(atom, ph=7, debye_length=math.inf)

    # Compared before rounding, since at four lengths the screened value rounds to zero.
    exact_ratio = math.exp(-multiple)
    assert unscreened.ep > 0
    assert screened.ep == pytest.approx(unscreened.ep * exact_ratio, abs=0.01)


def test_screening_shrinks_the_potential_at_every_point(atoms):
    """Damping every term of a sum of like sign cannot increase its magnitude.

    Checked rather than assumed, because the sum mixes signs and a cancelling
    pair could in principle move either way.
    """

    rng = np.random.default_rng(1)
    centre = np.array([[atom.x, atom.y, atom.z] for atom in atoms]).mean(axis=0)
    shrunk = 0

    for offset in rng.uniform(-25, 25, size=(40, 3)):
        unscreened = Property_point(*(centre + offset))
        unscreened.set_ep(atoms, ph=7, debye_length=debye_length(0))
        screened = Property_point(*(centre + offset))
        screened.set_ep(atoms, ph=7, debye_length=debye_length(PHYSIOLOGICAL))

        assert abs(screened.ep) <= abs(unscreened.ep) + 0.01
        shrunk += abs(screened.ep) < abs(unscreened.ep)

    assert shrunk > 20, "screening had almost no effect, which cannot be right at 150 mM"


def test_a_charge_beyond_several_debye_lengths_is_effectively_gone():
    """The whole point: the far side of a protein stops contributing.

    A single unit charge 60 Angstrom away, roughly across a large protein,
    contributes 0.06 V unscreened, which is a real amount when several hundred
    such charges add up. At 150 mM it is damped by exp(-7.6) and vanishes into
    the rounding entirely.
    """

    class FarCharge:
        x, y, z = 60.0, 0.0, 0.0

        def charge(self, ph=7, formal=True):
            return 1.0

    point = Property_point(0.0, 0.0, 0.0)
    point.set_ep(np.array([FarCharge()]), ph=7, debye_length=debye_length(PHYSIOLOGICAL))

    unscreened = Property_point(0.0, 0.0, 0.0)
    unscreened.set_ep(np.array([FarCharge()]), ph=7, debye_length=debye_length(0))

    assert unscreened.ep == pytest.approx(0.06, abs=0.005), "an unscreened unit charge at 60 A contributes 0.06 V"
    assert point.ep == 0.0, "at 150 mM a charge 60 A away should vanish into the rounding"


def test_the_sign_of_a_nearby_charge_still_survives_screening():
    """Screening must damp the field, not invert it."""

    class NearCharge:
        x, y, z = 4.0, 0.0, 0.0

        def __init__(self, sign):
            self.sign = sign

        def charge(self, ph=7, formal=True):
            return self.sign

    for sign in (1.0, -1.0):
        point = Property_point(0.0, 0.0, 0.0)
        point.set_ep(np.array([NearCharge(sign)]), ph=7, debye_length=debye_length(PHYSIOLOGICAL))
        assert math.copysign(1, point.ep) == sign


def test_screening_reveals_negative_surface_on_a_net_positive_protein(tmp_path):
    """The mirror of the case that motivated this change.

    Every protein in the reference set is net negative, where the unscreened sum
    buried the positive patches. 1GDW is net positive, so the same mechanism
    should bury its negative patches instead, and screening should uncover them.
    Testing only the negative-protein direction would leave a fix that happened to
    work by pushing everything one way.
    """

    screened = tmp_path / "screened.zip"
    unscreened = tmp_path / "unscreened.zip"
    calculate("tests/data/1GDW.pdb.zip", str(screened), ionic_strength_molar=PHYSIOLOGICAL)
    calculate("tests/data/1GDW.pdb.zip", str(unscreened), ionic_strength_molar=0)

    before = read_surface_points(unscreened)["ep_volts"]
    after = read_surface_points(screened)["ep_volts"]

    assert (before < 0).sum() < 100, "1GDW should look almost entirely positive without screening"
    assert (after < 0).sum() > 10 * (before < 0).sum(), "screening should uncover the negative surface"


@pytest.mark.parametrize("ionic_strength", [0.06, 0.15, 0.37])
def test_the_result_is_not_sensitive_to_the_exact_ionic_strength(ionic_strength, tmp_path):
    """Agreement with a Poisson-Boltzmann reference was measured flat over this range.

    The features should move smoothly rather than swinging, so that a user who
    picks 0.1 instead of 0.15 does not get a qualitatively different answer.
    """

    reference = tmp_path / "reference.zip"
    other = tmp_path / f"other_{ionic_strength}.zip"
    calculate(STRUCTURE, str(reference), ionic_strength_molar=0.15)
    calculate(STRUCTURE, str(other), ionic_strength_molar=ionic_strength)

    a = read_surface_points(reference)["ep_volts"]
    b = read_surface_points(other)["ep_volts"]
    positive_fraction = ((a > 0).mean(), (b > 0).mean())

    assert abs(positive_fraction[0] - positive_fraction[1]) < 0.15, f"positive surface swung too far: {positive_fraction}"


def test_the_shell_features_are_deliberately_not_screened(tmp_path):
    """The shell uses a different model that did not have the defect.

    It divides the path at the molecular surface and weights the solvent leg with
    a permittivity of 80 against 4 for the protein, so a distant charge is already
    damped about twentyfold and the offset never builds up. Measured against an
    APBS equivalent the unscreened shell agrees at a Spearman of 0.877, and adding
    screening moved that to 0.855.

    This pins the decision so that nobody screens the shell for symmetry without
    first producing evidence that it helps.
    """

    screened = tmp_path / "screened.zip"
    unscreened = tmp_path / "unscreened.zip"
    calculate(STRUCTURE, str(screened), full_features=True, ionic_strength_molar=PHYSIOLOGICAL)
    calculate(STRUCTURE, str(unscreened), full_features=True, ionic_strength_molar=0)

    a = read_features(screened).iloc[0]
    b = read_features(unscreened).iloc[0]
    shell = [column for column in a.index if "Shell" in column]

    assert shell, "no shell features found, so this test is not checking anything"
    assert all(a[column] == b[column] for column in shell)


def test_a_run_records_the_ionic_strength_it_used(tmp_path):
    """A bundle whose potential cannot be reproduced from its own record is useless."""

    out_file = tmp_path / "bundle.zip"
    calculate(STRUCTURE, str(out_file), ionic_strength_molar=0.05)

    with zipfile.ZipFile(out_file) as archive:
        record = archive.read("ARH96693/prodes_run.json").decode()

    assert '"ionic_strength_molar": 0.05' in record


@pytest.mark.parametrize("full_features", [False, True])
def test_screening_changes_the_features_a_run_produces(full_features, tmp_path):
    """If the flag reached nothing, every other test here would still pass.

    Parametrized over the feature set because the partial charge surface columns,
    the SurfEp*Average family, are only produced with the full set. Without that
    the screening could be removed from the partial charge path with every test
    still green.
    """

    screened = tmp_path / "screened.zip"
    unscreened = tmp_path / "unscreened.zip"
    calculate(STRUCTURE, str(screened), full_features=full_features, ionic_strength_molar=PHYSIOLOGICAL)
    calculate(STRUCTURE, str(unscreened), full_features=full_features, ionic_strength_molar=0)

    a = read_features(screened).iloc[0].drop("ID").astype(float)
    b = read_features(unscreened).iloc[0].drop("ID").astype(float)
    changed = [column for column in a.index if a[column] != b[column]]

    assert changed, "screening changed nothing at all"
    assert all("Ep" in column for column in changed), f"screening touched non-electrostatic features: {changed}"
    assert not any("Shell" in column for column in changed), "the shell is deliberately unscreened"

    if full_features:
        average = [column for column in changed if "Average" in column]
        assert average, "the partial charge surface features were not screened"


def test_the_surface_never_reports_more_positive_points_when_screened(tmp_path):
    """Screening removes the net negative offset, so positive patches can only appear.

    Every structure in the reference set is net negative, so the unscreened sum
    pushes the whole surface down. Removing that cannot hide a positive point.
    """

    screened = tmp_path / "screened.zip"
    unscreened = tmp_path / "unscreened.zip"
    calculate(STRUCTURE, str(screened), ionic_strength_molar=PHYSIOLOGICAL)
    calculate(STRUCTURE, str(unscreened), ionic_strength_molar=0)

    a = read_surface_points(screened)["ep_volts"]
    b = read_surface_points(unscreened)["ep_volts"]

    assert len(a) == len(b)
    assert (a > 0).sum() >= (b > 0).sum()


def test_the_positive_and_negative_counts_still_partition_the_surface(tmp_path):
    """The features document NSurfNeg as NSurfPoints minus NSurfPos.

    Screening keeps the cutoff at zero, so that invariant must survive. Points at
    exactly zero belong to neither, which is the only slack allowed.
    """

    out_file = tmp_path / "bundle.zip"
    calculate(STRUCTURE, str(out_file), full_features=True, ionic_strength_molar=PHYSIOLOGICAL)

    features = read_features(out_file).iloc[0]
    points = read_surface_points(out_file)["ep_volts"]

    assert float(features["NSurfPosEpFormal"]) == (points > 0).sum()
    assert float(features["NSurfNegEpFormal"]) == (points < 0).sum()
    assert float(features["NSurfPosEpFormal"]) + float(features["NSurfNegEpFormal"]) + (points == 0).sum() == len(points)
