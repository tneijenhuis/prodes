"""Tests for finding disulfide bonds and for the charge they must not carry.

Prodes used to give every CYS the free-thiol pKa of 8.33 whatever the structure
looked like, so both halves of every disulfide titrated. On 1GDW, which is all
disulfides and no free thiol, that put the formal charge at pH 8.5 at -1 when it
should be +7 and moved the isoelectric point by 1.4 units. Issue #6.

Three structures carry most of the load, and between them they cover every case
the detection has to get right:

    1GDW    8 cysteines, all 4 bonds present, no free thiol, no SSBOND records
    1GPB    9 cysteines, none bonded, no SSBOND records
    1AO6    70 cysteines over two chains: 34 bonds and exactly two free thiols,
            with 34 SSBOND records that agree with the geometry

1AO6, serum albumin, is the one that answers the question the others cannot: a
free cysteine and a bonded one in the same structure, so nothing but the bond
differs between them. Its free Cys 34 is the well known one. It is also the
first multi-chain fixture in the repository, so it is what covers resolving a
record by chain as well as by residue number.

tests/data/1AO6.pdb.zip was downloaded from files.rcsb.org and trimmed to its
CRYST1, SSBOND, ATOM, TER and END records, which took it from 780 KB to 187 KB
and changed nothing the parser reads. Its 34 SSBOND records agree exactly with
what the geometry finds, which is the cross-check in
test_geometry_alone_reproduces_the_albumin_records below.

The synthetic cases are written out as files here rather than committed, so that
the doctoring a test depends on is visible in the test itself.
"""

import logging

import pytest

from prodes.calculations.disulfides import (
    MAX_SG_SG_DISTANCE_ANGSTROM,
    assign_disulfides,
    candidate_pairs,
    cysteine_sulfurs,
    geometric_disulfides,
)
from prodes.io.parser import PDBparser, read_ssbond_line

LYSOZYME = "tests/data/1GDW.pdb.zip"
PHOSPHORYLASE = "tests/data/1GPB.pdb.zip"
ALBUMIN = "tests/data/1AO6.pdb.zip"


def parsed(path):
    """Returns the parsed structure, which already has its disulfides assigned."""

    return PDBparser().parse(path)


def residues_by_chain_and_number(structure):
    """Returns {(chain, number): residue}, which is how a PDB file names a residue."""

    return {(residue.chain.name, residue.number): residue for residue in structure.residues}


def bonded_numbers(structure):
    """Returns the detected bonds as a set of residue number pairs."""

    return {frozenset((first.number, second.number)) for first, second in structure.disulfides}


def pdb_line(fields):
    """Returns an 80 column PDB record built from (start column, text) pairs.

    Columns are given as the format's own 1-based numbers, so they can be read
    straight off the specification rather than counted out in an f-string.
    """

    line = [" "] * 80
    for column, text in fields:
        line[column - 1 : column - 1 + len(text)] = text

    return "".join(line).rstrip()


def atom_line(serial, name, residue_name, chain, number, x, y, z, element):
    """Returns one ATOM record, so a test can write the structure it needs."""

    return pdb_line(
        [
            (1, "ATOM"),
            (7, f"{serial:5d}"),
            (13, f"{name:<4s}"),
            (18, f"{residue_name:>3s}"),
            (22, chain),
            (23, f"{number:4d}"),
            (31, f"{x:8.3f}"),
            (39, f"{y:8.3f}"),
            (47, f"{z:8.3f}"),
            (55, "  1.00"),
            (61, "  0.00"),
            (77, f"{element:>2s}"),
        ]
    )


def cysteine_lines(serial, chain, number, x):
    """Returns the six heavy atoms of one cysteine, placed with its SG at x.

    The other atoms only have to be somewhere sensible: nothing under test reads
    them, but a residue with no backbone is not a residue the parser would ever
    produce.
    """

    return [
        atom_line(serial, "N", "CYS", chain, number, x - 3.0, 0.0, 0.0, "N"),
        atom_line(serial + 1, "CA", "CYS", chain, number, x - 2.0, 0.0, 0.0, "C"),
        atom_line(serial + 2, "C", "CYS", chain, number, x - 2.0, 1.5, 0.0, "C"),
        atom_line(serial + 3, "O", "CYS", chain, number, x - 2.0, 2.5, 0.0, "O"),
        atom_line(serial + 4, "CB", "CYS", chain, number, x - 1.0, 0.0, 0.0, "C"),
        atom_line(serial + 5, "SG", "CYS", chain, number, x, 0.0, 0.0, "S"),
    ]


def write_cysteine_structure(path, sulfur_positions, records=()):
    """Writes a PDB file of cysteines whose SG atoms sit at the given x positions.

    Args:
        path: where to write it.
        sulfur_positions: {(chain, residue number): x coordinate of the SG}.
        records: SSBOND lines to put at the top of the file.
    """

    lines = list(records)
    serial = 1
    for (chain, number), x in sulfur_positions.items():
        lines.extend(cysteine_lines(serial, chain, number, x))
        serial += 6

    path.write_text("\n".join(lines) + "\nEND\n")

    return str(path)


def ssbond_line(chain1, number1, chain2, number2, symmetry1="1555", symmetry2="1555"):
    """Returns one SSBOND record in the column layout the format specifies."""

    return pdb_line(
        [
            (1, "SSBOND"),
            (8, "  1"),
            (12, "CYS"),
            (16, chain1),
            (18, f"{number1:4d}"),
            (26, "CYS"),
            (30, chain2),
            (32, f"{number2:4d}"),
            (60, f"{symmetry1:>6s}"),
            (67, f"{symmetry2:>6s}"),
            (74, " 2.03"),
        ]
    )


# Reading the record


def test_read_ssbond_line_reads_a_real_record():
    """A deposited record gives both partners and both symmetry operators."""

    line = "SSBOND   1 CYS A    6    CYS A  128                          1555   1555  2.03  "

    assert read_ssbond_line(line) == ("A", 6, "A", 128, "1555", "1555")


def test_read_ssbond_line_keeps_a_symmetry_operator():
    """A non-identity operator has to survive parsing, since it is what disqualifies the record."""

    line = "SSBOND   1 CYS A   73    CYS A   73                          2655   1555  2.05  "

    assert read_ssbond_line(line) == ("A", 73, "A", 73, "2655", "1555")


def test_read_ssbond_line_treats_a_truncated_line_as_the_identity_operator():
    """Trimmed and hand-written files stop before column 72, and must not lose their bonds."""

    line = "SSBOND   1 CYS A    6    CYS A  128"

    assert read_ssbond_line(line) == ("A", 6, "A", 128, "1555", "1555")


def test_read_ssbond_line_returns_none_for_a_line_it_cannot_read():
    """One malformed record costs its own bond, not the whole file."""

    assert read_ssbond_line("SSBOND   1 CYS A    x    CYS A    y") is None


# Detection on real structures


def test_lysozyme_has_four_disulfides():
    """All eight cysteines of 1GDW are bonded, in four pairs, with no SSBOND records to help."""

    structure = parsed(LYSOZYME)

    assert len(structure.disulfides) == 4
    assert bonded_numbers(structure) == {frozenset(pair) for pair in ((6, 128), (30, 116), (65, 81), (77, 95))}


def test_the_partner_relation_is_symmetric():
    """Both halves of a bond know about it, so neither can be titrated by accident."""

    for first, second in parsed(LYSOZYME).disulfides:
        assert first.disulfide_partner is second
        assert second.disulfide_partner is first


def test_free_cysteines_are_not_paired():
    """1GPB has nine cysteines and no disulfide at all; proximity alone must not pair them."""

    structure = parsed(PHOSPHORYLASE)

    assert sum(1 for residue in structure.residues if residue.name == "CYS") == 9
    assert structure.disulfides == []


def test_albumin_is_detected_from_its_ssbond_records():
    """A structure that carries records: 34 bonds, and the two free thiols left alone."""

    structure = parsed(ALBUMIN)
    free = [residue for residue in structure.residues if residue.name == "CYS" and residue.disulfide_partner is None]

    assert len(structure.disulfides) == 34
    assert [(residue.chain.name, residue.number) for residue in free] == [("A", 34), ("B", 34)]


def test_geometry_alone_reproduces_the_albumin_records():
    """The two routes agree on real data, which is the check that neither is drifting.

    Albumin is the only shipped structure with both SSBOND records and enough
    disulfides for the comparison to mean anything.
    """

    structure = parsed(ALBUMIN)
    from_records = bonded_numbers(structure)

    for residue in structure.residues:
        residue.disulfide_partner = None
    from_geometry = {frozenset((first.number, second.number)) for first, second in geometric_disulfides(cysteine_sulfurs(structure))}

    assert from_geometry == from_records


# The distance criterion


def test_a_pair_just_inside_the_cutoff_is_called(tmp_path):
    """2.45 A is longer than a textbook bond but inside what real structures model."""

    path = write_cysteine_structure(tmp_path / "close.pdb", {("A", 1): 0.0, ("A", 2): 2.45})

    assert len(parsed(path).disulfides) == 1


def test_a_pair_just_outside_the_cutoff_is_not(tmp_path):
    """Past 2.5 A prodes stops guessing, which is where PROPKA stops too."""

    path = write_cysteine_structure(tmp_path / "far.pdb", {("A", 1): 0.0, ("A", 2): 2.55})

    assert parsed(path).disulfides == []


def test_coincident_sulfurs_are_reported_rather_than_bonded(tmp_path, caplog):
    """Two sulfurs half an Angstrom apart are duplicated atoms, not a very short bond."""

    path = write_cysteine_structure(tmp_path / "coincident.pdb", {("A", 1): 0.0, ("A", 2): 0.5})

    with caplog.at_level(logging.WARNING):
        structure = parsed(path)

    assert structure.disulfides == []
    assert "too short for a bond" in caplog.text


def test_the_shorter_bond_wins_a_contested_sulfur(tmp_path):
    """Three sulfurs in a row: the middle one can only bond once, and takes the nearer partner."""

    path = write_cysteine_structure(tmp_path / "contested.pdb", {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 4.45})

    structure = parsed(path)

    assert bonded_numbers(structure) == {frozenset((1, 2))}
    assert residues_by_chain_and_number(structure)[("A", 3)].disulfide_partner is None


def test_candidate_pairs_come_out_shortest_first():
    """The greedy assignment above is only correct because the candidates are ordered."""

    distances = [distance for distance, _, _ in candidate_pairs(cysteine_sulfurs(parsed(ALBUMIN)))]

    assert distances == sorted(distances)
    assert max(distances) <= MAX_SG_SG_DISTANCE_ANGSTROM


# Records against geometry


def test_a_record_is_honoured_even_when_the_sulfurs_are_far_apart(tmp_path, caplog):
    """The depositor's chemistry outranks our window, but the disagreement is reported."""

    records = [ssbond_line("A", 1, "A", 2)]
    path = write_cysteine_structure(tmp_path / "stretched.pdb", {("A", 1): 0.0, ("A", 2): 15.0}, records)

    with caplog.at_level(logging.WARNING):
        structure = parsed(path)

    assert bonded_numbers(structure) == {frozenset((1, 2))}
    assert "check whether the structure is reduced" in caplog.text


def test_a_symmetry_record_does_not_suppress_the_real_bonds(tmp_path):
    """The 1ERT case: its only record joins Cys 73 to itself across a symmetry axis.

    Records are authoritative for the cysteines they name, not for the file. A
    record that names nothing in these coordinates must therefore leave every
    other cysteine to the geometry rather than switching it off.
    """

    records = [ssbond_line("A", 3, "A", 3, symmetry1="2655")]
    path = write_cysteine_structure(tmp_path / "symmetry.pdb", {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 20.0}, records)

    assert bonded_numbers(parsed(path)) == {frozenset((1, 2))}


def test_an_incomplete_record_set_keeps_the_bonds_it_forgot(tmp_path):
    """One record present must not cost the other pair its bond."""

    records = [ssbond_line("A", 1, "A", 2)]
    positions = {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 10.0, ("A", 4): 12.05}
    path = write_cysteine_structure(tmp_path / "partial.pdb", positions, records)

    assert bonded_numbers(parsed(path)) == {frozenset((1, 2)), frozenset((3, 4))}


def test_a_record_wins_the_cysteines_it_names(tmp_path):
    """A cysteine a record has claimed is not re-paired by geometry with someone nearer."""

    records = [ssbond_line("A", 1, "A", 3)]
    positions = {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 20.0}
    path = write_cysteine_structure(tmp_path / "claimed.pdb", positions, records)

    structure = parsed(path)

    assert bonded_numbers(structure) == {frozenset((1, 3))}
    assert residues_by_chain_and_number(structure)[("A", 2)].disulfide_partner is None


def test_a_record_naming_a_missing_residue_is_ignored(tmp_path, caplog):
    """A record can outlive the coordinates it refers to; that must not stop the run."""

    records = [ssbond_line("A", 1, "A", 99)]
    path = write_cysteine_structure(tmp_path / "missing.pdb", {("A", 1): 0.0, ("A", 2): 2.05}, records)

    with caplog.at_level(logging.WARNING):
        structure = parsed(path)

    assert "not both cysteines in this file" in caplog.text
    assert bonded_numbers(structure) == {frozenset((1, 2))}


def test_an_interchain_record_is_honoured(tmp_path):
    """Residue numbers repeat across chains, so a bond has to be resolved by chain as well."""

    records = [ssbond_line("A", 1, "B", 1)]
    path = write_cysteine_structure(tmp_path / "interchain.pdb", {("A", 1): 0.0, ("B", 1): 2.05}, records)

    structure = parsed(path)
    residues = residues_by_chain_and_number(structure)

    assert len(structure.disulfides) == 1
    assert residues[("A", 1)].disulfide_partner is residues[("B", 1)]


# What the bond does to the charge


@pytest.mark.parametrize(
    ("ph", "formal", "expected"),
    [
        # At pH 7 a free thiol at pKa 8.33 is protonated anyway, so only the
        # fractional charge moves. By pH 8.5 the formal charge had the wrong sign.
        (7.0, True, 7.0),
        (7.0, False, 7.092),
        (8.5, True, 7.0),
        (8.5, False, 6.709),
        (9.0, True, 7.0),
        (9.0, False, 6.133),
    ],
)
def test_lysozyme_charge_is_no_longer_dragged_down_by_its_cystines(ph, formal, expected):
    """Before detection these read 7.0, 6.735, -1.0, 1.936, -1.0 and -0.458."""

    assert round(parsed(LYSOZYME).charge(ph, formal=formal), 3) == expected


def test_lysozyme_isoelectric_point():
    """The pI is the number a user is most likely to notice: it was 8.899."""

    assert parsed(LYSOZYME).isoelectric_point() == 10.342


def test_a_bonded_cysteine_at_a_chain_terminus_keeps_its_terminal_charge(tmp_path):
    """The case that crashes if the pKa is read by list position rather than by name.

    With no side-chain entry left, position 0 of the list is the terminus, so
    the SG would take the C-terminal pKa of 2.34 and come out fully charged
    above pH 2.3, or the lookup would fail outright on a residue with no
    entries at all.
    """

    path = write_cysteine_structure(tmp_path / "terminal.pdb", {("A", 1): 0.0, ("A", 2): 2.05})
    structure = parsed(path)
    residues = residues_by_chain_and_number(structure)

    first, last = residues[("A", 1)], residues[("A", 2)]
    assert (first.terminus, last.terminus) == ("N", "C")

    # The side chains are gone but the backbone groups are untouched.
    assert first.pkas == [{"N+": 9.69}]
    assert last.pkas == [{"C-": 2.34}]

    sulfurs = {residue.number: [atom for atom in residue.atoms if atom.name == "SG"][0] for residue in (first, last)}
    assert sulfurs[1].charge(7) == 0
    assert sulfurs[2].charge(7) == 0

    assert round(first.charge(7), 3) == 1
    assert round(last.charge(7), 3) == -1


def test_assign_disulfides_can_be_called_after_the_pkas_have_been_read():
    """Marking a residue clears its cached pKas, so detection is not order dependent."""

    structure = parsed(LYSOZYME)
    cysteine = next(residue for residue in structure.residues if residue.name == "CYS")

    for residue in structure.residues:
        residue.disulfide_partner = None
        residue._pka = None
    assert cysteine.pkas == [{"CYS": 8.33}], "the cysteine should look titratable once its partner is forgotten"

    assign_disulfides(structure)

    assert cysteine.pkas is None


# The record filters, one test each. Written separately on purpose: a single
# record that is both a self-pair and under a symmetry operator is caught twice,
# and each filter then covers for the other being deleted.


def test_a_symmetry_operator_disqualifies_a_record_between_two_residues(tmp_path):
    """Not a self-pair, so only the symmetry field can rule it out.

    The two partners are 20 A apart, which is what the record would have to
    override for the bond to appear.
    """

    records = [ssbond_line("A", 1, "A", 3, symmetry1="2655")]
    path = write_cysteine_structure(tmp_path / "symmetric_pair.pdb", {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 20.0}, records)

    assert bonded_numbers(parsed(path)) == {frozenset((1, 2))}


def test_a_self_pair_is_dropped_even_under_the_identity_operator(tmp_path):
    """A cysteine cannot bond to itself, whatever the record's symmetry fields say."""

    records = [ssbond_line("A", 3, "A", 3)]
    path = write_cysteine_structure(tmp_path / "self_pair.pdb", {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 20.0}, records)

    structure = parsed(path)

    assert bonded_numbers(structure) == {frozenset((1, 2))}
    assert residues_by_chain_and_number(structure)[("A", 3)].disulfide_partner is None


def test_a_second_record_claiming_the_same_cysteine_is_refused(tmp_path, caplog):
    """Contradictory records: the first wins, and the conflict is reported."""

    records = [ssbond_line("A", 1, "A", 2), ssbond_line("A", 1, "A", 3)]
    path = write_cysteine_structure(tmp_path / "contradictory.pdb", {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 20.0}, records)

    with caplog.at_level(logging.WARNING):
        structure = parsed(path)

    assert "claim the same cysteine" in caplog.text
    assert bonded_numbers(structure) == {frozenset((1, 2))}


def test_a_residue_with_two_sg_atoms_is_reported_and_counted_once(tmp_path, caplog):
    """Two residues sharing a number are merged by the parser, which this must survive.

    The merge is a separate defect, of insertion code handling, and it makes the
    residue's charge wrong however this module behaves. What it must not do is
    let one residue take two partners.
    """

    lines = cysteine_lines(1, "A", 1, 0.0) + cysteine_lines(7, "A", 1, 2.05) + cysteine_lines(13, "A", 2, 4.10)
    path = tmp_path / "merged.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")

    with caplog.at_level(logging.WARNING):
        structure = parsed(str(path))

    assert "carries 2 SG atoms" in caplog.text
    assert len(cysteine_sulfurs(structure)) == 2
    assert len(structure.disulfides) <= 1


def test_a_record_works_when_the_file_gives_no_chain_id(tmp_path):
    """Some hand-built and modelling-tool files leave the chain column blank."""

    records = [ssbond_line(" ", 1, " ", 2)]
    path = write_cysteine_structure(tmp_path / "no_chain.pdb", {(" ", 1): 0.0, (" ", 2): 15.0}, records)

    assert bonded_numbers(parsed(path)) == {frozenset((1, 2))}


def test_the_bundle_records_how_many_bonds_were_found(tmp_path):
    """End to end: the count a user checks reaches prodes_run.json from a real run."""

    import json
    import zipfile

    from prodes.run import calculate

    bundle = calculate(LYSOZYME, str(tmp_path / "1GDW.zip"))

    with zipfile.ZipFile(bundle) as archive:
        name = next(member for member in archive.namelist() if member.endswith("prodes_run.json"))
        record = json.loads(archive.read(name))

    assert record["disulfides"] == 4


def test_running_detection_twice_gives_the_same_answer(tmp_path):
    """Re-running must not leave the first run's partners behind.

    The bonds are recomputed from scratch each time, so a stale partner would
    show up as a residue that is still marked bonded while structure.disulfides
    no longer lists it, and that residue would stay untitratable for good.
    """

    records = [ssbond_line("A", 1, "A", 3)]
    path = write_cysteine_structure(tmp_path / "twice.pdb", {("A", 1): 0.0, ("A", 2): 2.05, ("A", 3): 20.0}, records)
    structure = parsed(path)

    assign_disulfides(structure)

    bonded = [residue for residue in structure.residues if residue.disulfide_partner is not None]
    assert len(structure.disulfides) == 1
    assert len(bonded) == 2 * len(structure.disulfides)
    # Without the records this time, so geometry decides and the answer changes.
    assert bonded_numbers(structure) == {frozenset((1, 2))}
