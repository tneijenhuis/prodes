"""Finds which cysteines of a structure are joined by a disulfide bond.

A cysteine whose SG is bonded to another cysteine's SG has no thiol proton, so
it does not titrate at all. Prodes used to give every CYS the free-thiol pKa of
8.33 whatever the structure looked like, which put a spurious negative charge on
both halves of every disulfide above about pH 8. On lysozyme that moved the
formal charge at pH 8.5 from +7 to -1 and the isoelectric point by 1.4 units.

Two routes to the answer, in this order:

    SSBOND records   what the depositor says, authoritative for the two
                     cysteines each record names
    geometry         the SG-SG distance, for every cysteine no record claimed

Records are honoured per cysteine rather than per file. A file whose record set
is incomplete keeps the bonds its records forgot, and a file whose only records
are crystallographic symmetry pairs still gets its real bonds. 1ERT is the case
that forces this: its single SSBOND record joins Cys 73 to itself across a
symmetry axis and names nothing in the file's own coordinates.

AlphaFold models carry no SSBOND records at all, so for them the geometric route
is the only one, and it works: they place the SG atoms at bonding distance.

What this cannot see, and what therefore stays wrong, is set out under "Cysteines
and disulfide bonds" in the README: a cysteine bound to a metal, a thioether link
to a haem, and a real bond that a low-confidence model has stretched too far.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

# The longest SG-SG separation still called a bond. PROPKA 3.5.1 and PDB2PQR
# both use 2.5 A, and using their number means prodes cannot disagree with the
# pKa predictor the README tells users to run about which cysteines are bridged.
#
# It is deliberately not wider. Real bonds run to about 2.36 A: 6VXX, a 2.8 A
# cryo-EM structure, models Cys391-Cys525 at 2.361 A in all three protomers,
# which is why the tighter 2.35 A window this started with was abandoned. But the
# AlphaFold model of metallothionein-2, a protein with twenty metal-binding
# cysteines and no disulfides at all, places two of them 2.050 A apart with the
# rest of the cluster at 3.3 A and beyond. No distance cutoff excludes that false
# positive, and widening towards the cluster would turn one into a dozen.
MAX_SG_SG_DISTANCE_ANGSTROM = 2.5

# Below this the two sulfurs are coincident rather than bonded, which happens in
# modelling output that has duplicated a chain or written two conformations
# without alternate location indicators.
MIN_SG_SG_DISTANCE_ANGSTROM = 1.6

# The symmetry operator meaning "this atom, in this asymmetric unit". Any other
# value on an SSBOND record names a crystallographic copy that is not among the
# file's coordinates.
IDENTITY_SYMMETRY_OPERATOR = "1555"


def cysteine_sulfurs(structure):
    """Returns the SG atom of every CYS residue, in parse order.

    A residue with more than one SG is reported and contributes its first: that
    is the visible symptom of the parser merging two residues that differ only
    by an insertion code, which is a separate defect and makes the charge of that
    residue wrong whatever this module does.
    """

    sulfurs = []
    for residue in structure.residues:
        if residue.name != "CYS":
            continue

        found = [atom for atom in residue.atoms if atom.name == "SG"]
        if len(found) > 1:
            logger.warning(
                "%s carries %d SG atoms, so two residues have been read as one; its charge and any disulfide it makes are unreliable",
                residue_label(residue),
                len(found),
            )
        sulfurs.extend(found[:1])

    return sulfurs


def residue_label(residue):
    """Returns a short readable name for a residue, such as 'CYS A30'."""

    chain = residue.chain.name if residue.chain is not None else ""

    return f"{residue.name} {chain}{residue.number}"


def sulfur_distance(first, second):
    """Returns the distance in Angstrom between two atoms."""

    return float(np.linalg.norm(np.array([first.x, first.y, first.z]) - np.array([second.x, second.y, second.z])))


def candidate_pairs(sulfurs):
    """Returns (distance, first, second) for every SG pair close enough to be a bond, shortest first.

    Pairs closer than MIN_SG_SG_DISTANCE_ANGSTROM are reported and left out:
    they are coincident atoms rather than a short bond.
    """

    if len(sulfurs) < 2:
        return []

    coordinates = np.array([[atom.x, atom.y, atom.z] for atom in sulfurs])
    separations = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=-1)

    candidates = []
    for i, j in zip(*np.triu_indices(len(sulfurs), k=1), strict=True):
        distance = float(separations[i, j])
        if distance > MAX_SG_SG_DISTANCE_ANGSTROM:
            continue

        if distance < MIN_SG_SG_DISTANCE_ANGSTROM:
            logger.warning(
                "%s and %s are only %.2f A apart, which is too short for a bond; treating them as coincident atoms rather than a disulfide",
                residue_label(sulfurs[i].residue),
                residue_label(sulfurs[j].residue),
                distance,
            )
            continue

        candidates.append((distance, sulfurs[i], sulfurs[j]))

    return sorted(candidates, key=lambda candidate: candidate[0])


def geometric_disulfides(sulfurs):
    """Pairs SG atoms by distance and returns the residue pairs.

    Each sulfur takes at most one partner. Where two candidate bonds compete for
    the same sulfur the shorter one wins, which is why the candidates are walked
    shortest first.
    """

    paired = set()
    bonds = []
    for _, first, second in candidate_pairs(sulfurs):
        if id(first) in paired or id(second) in paired:
            continue

        paired.update((id(first), id(second)))
        bonds.append((first.residue, second.residue))

    return bonds


def residues_by_chain_and_number(structure):
    """Returns {(chain name, residue number): residue}, which is how SSBOND records name residues."""

    return {(residue.chain.name if residue.chain is not None else "", residue.number): residue for residue in structure.residues}


def record_disulfides(structure, records):
    """Resolves SSBOND records to residue pairs, dropping the ones that name nothing here.

    A record is skipped when it carries a symmetry operator other than the
    identity, when it joins a residue to itself, or when either side is missing
    from the coordinates or is not a cysteine. A record that survives is honoured
    even if the two sulfurs are further apart than a bond should be, because the
    depositor's chemistry outranks our distance window, but the disagreement is
    reported: a stale record on a reduced or mutated model is the one case where
    trusting it silently would be wrong.

    Args:
        structure: the parsed Structure.
        records: tuples of (chain1, number1, chain2, number2, sym1, sym2), as
            read by prodes.io.parser.read_ssbond_line.
    """

    residues = residues_by_chain_and_number(structure)

    bonds = []
    claimed = set()
    for chain1, number1, chain2, number2, symmetry1, symmetry2 in records:
        if symmetry1 != IDENTITY_SYMMETRY_OPERATOR or symmetry2 != IDENTITY_SYMMETRY_OPERATOR:
            # One partner sits in a neighbouring asymmetric unit, so the bond is
            # between two crystallographic copies and not inside this molecule.
            logger.debug("skipping an SSBOND record under symmetry operators %s and %s, which is a crystal contact", symmetry1, symmetry2)
            continue

        if (chain1, number1) == (chain2, number2):
            continue

        first, second = residues.get((chain1, number1)), residues.get((chain2, number2))
        if first is None or second is None or first.name != "CYS" or second.name != "CYS":
            logger.warning(
                "an SSBOND record joins %s%s to %s%s, which are not both cysteines in this file; ignoring it",
                chain1,
                number1,
                chain2,
                number2,
            )
            continue

        if id(first) in claimed or id(second) in claimed:
            logger.warning("two SSBOND records claim the same cysteine, %s or %s; keeping the first", residue_label(first), residue_label(second))
            continue

        report_record_distance(first, second)
        claimed.update((id(first), id(second)))
        bonds.append((first, second))

    return bonds


def report_record_distance(first, second):
    """Warns when an SSBOND record joins two sulfurs that are not at bonding distance."""

    sulfurs = [[atom for atom in residue.atoms if atom.name == "SG"] for residue in (first, second)]
    if not all(sulfurs):
        return

    distance = sulfur_distance(sulfurs[0][0], sulfurs[1][0])
    if not MIN_SG_SG_DISTANCE_ANGSTROM <= distance <= MAX_SG_SG_DISTANCE_ANGSTROM:
        logger.warning(
            "an SSBOND record joins %s to %s, which are %.2f A apart rather than about 2.05 A; honouring the record, but check whether the structure is reduced",
            residue_label(first),
            residue_label(second),
            distance,
        )


def assign_disulfides(structure, records=()):
    """Finds the structure's disulfide bonds and records them on its residues.

    Both members of every bond get a disulfide_partner, which stops
    Residue.pkas offering a titratable side chain and so removes the charge
    everywhere it would otherwise appear. The list of pairs is also stored on
    structure.disulfides so a caller can report how many were found.

    Meant to run on a freshly parsed structure, which is where the parser calls
    it, but it does not depend on being first: any previous assignment is
    cleared and every cysteine's cached pKas go with it, so a second call gives
    the same answer as the first. Note that clearing the cache also discards
    predicted values a pKa file has already supplied for a cysteine, so run this
    before redo_pkas rather than after it.

    Args:
        structure: the parsed Structure to annotate.
        records: SSBOND records read from the same file, if it had any.

    Returns:
        The list of (residue, residue) pairs, records first.
    """

    sulfurs = cysteine_sulfurs(structure)
    for sulfur in sulfurs:
        sulfur.residue.mark_disulfide_partner(None)

    bonds = record_disulfides(structure, records)
    claimed = {id(residue) for bond in bonds for residue in bond}
    bonds.extend(geometric_disulfides([sulfur for sulfur in sulfurs if id(sulfur.residue) not in claimed]))

    for first, second in bonds:
        first.mark_disulfide_partner(second)
        second.mark_disulfide_partner(first)

    structure.disulfides = bonds

    return bonds
