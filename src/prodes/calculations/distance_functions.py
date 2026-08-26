import math

# Debye length of water at 298 K, in Angstrom per square root of molar ionic
# strength. This is the screening length in bulk water, where the relative
# permittivity is about 78.5.
DEBYE_COEFFICIENT_ANGSTROM = 3.04


def debye_length(ionic_strength_molar):
    """Returns the Debye screening length in Angstrom for a 1:1 salt.

    The length over which mobile buffer ions cancel a charge's field. At 150 mM,
    roughly physiological, it is about 7.9 Angstrom, so a charge on the far side
    of a protein contributes almost nothing to a surface point.

    Zero ionic strength means no mobile ions and therefore no screening, which is
    returned as an infinite length.

    What this is, stated plainly, because it is easy to overclaim. This is the
    Debye length in *water*, where the relative permittivity is about 78.5.
    Prodes evaluates its Coulomb sum at a uniform 4, where the self-consistent
    value would be about 1.8 Angstrom. Worse, the path from a buried charge to a
    surface point runs largely through protein interior, which contains no mobile
    ions to do any screening at all. So exp(-d / debye_length) is better described
    as a locality window, which happens to have close to the right width, than as
    Debye screening. It is used because it was measured to agree with a
    Poisson-Boltzmann reference far better than no damping at all, not because
    the derivation holds.

    The full linearised Debye-Huckel expression also carries an ion-size term,
    exp(kappa*a) / (1 + kappa*a). That is independent of distance, so it scales
    every contribution by the same factor, 1.03 for an ion radius of 2 Angstrom,
    and changes no ranking and no sign. It is omitted deliberately rather than
    overlooked.

    Args:
        ionic_strength_molar: ionic strength in mol/L. 0 disables screening.
    """

    if not math.isfinite(ionic_strength_molar):
        raise ValueError(f"ionic strength must be a finite number, got {ionic_strength_molar}")

    if ionic_strength_molar < 0:
        raise ValueError(f"ionic strength cannot be negative, got {ionic_strength_molar}")

    if ionic_strength_molar == 0:
        return math.inf

    return DEBYE_COEFFICIENT_ANGSTROM / math.sqrt(ionic_strength_molar)


def atom_charge_coulomb(charge):
    return charge * 1.6e-19


def charge_simple(charge, distance, dielectric_constant):
    """takes the charge in coulomb, distance in meters and dielectric constant to calculate the potential"""

    absolute_permittivity = 8.854e-12
    permittivity = dielectric_constant * absolute_permittivity

    return charge / (permittivity * distance * 4 * math.pi)


def potential_multiple_media(charge: float, distances: list, dielectic_constants: list):
    """takes the charge in coulomb, distances in meters and dielectric constants to calculate the potential at a point in space"""

    absolute_permittivity = 8.854e-12
    denominator = 0
    for distance, dielectric_constant in zip(distances, dielectic_constants, strict=True):
        permittivity = dielectric_constant * absolute_permittivity
        denominator += permittivity * distance

    return charge / (denominator * 4 * math.pi)
