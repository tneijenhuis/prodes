"""Tests for loading externally predicted pKa values into a structure.

prodes carries one default pKa per titratable residue type, and a
structure-based predictor such as PROPKA, H++ or pypka can override them per
residue. Nothing exercised that path: convert_propka and its two siblings had
no test, Structure.redo_pkas had no test, and the only tests that mention a pKa
file pass a made up filename to a mocked calculate(), which pins the argument
plumbing and nothing else.

It is not a corner of the package that can be left to look after itself. On
1GDW a real PROPKA prediction moves 19 of the 54 reduced features, formal
charge included, so a refactor that quietly stopped applying the values would
change more than a third of the output while every other test still passed.
Downstream code already depends on it: biochai runs PROPKA and imports
convert_propka to read the result, so this is a public API with an outside
caller.

The chain has four links and one test each below::

    PROPKA output -> convert_propka -> JSON -> read_pka -> redo_pkas -> features

tests/data/1GDW.pka is a genuine propka3.4.0 run of 1GDW, committed by the
original author, so the whole chain is covered from the repository alone. None
of this needs PROPKA, or any other predictor, to be installed.
"""

import logging

from prodes.core.structure import Structure
from prodes.io.parser import PDBparser, read_pka
from prodes.io.pka_converter import PROPKA_NOT_TITRATABLE, convert_propka, write_json
from tests.pipeline_output import pipeline_output

PDB_PATH = "tests/data/1GDW.pdb.zip"
PROPKA_PATH = "tests/data/1GDW.pka"


def test_convert_propka_reads_a_real_propka_file():
    """The summary table of a propka3 run becomes a residue number to pKa mapping.

    Values are read out of the file by column position, so the assertions name
    residues at both ends of the chain: a shifted slice would still parse and
    still return plausible looking numbers.
    """

    pkas = convert_propka(PROPKA_PATH)

    assert len(pkas) == 45
    assert pkas[18] == [{"ASP": 3.92}]
    assert pkas[35] == [{"GLU": 7.28}]

    # The termini are the two entries propka names differently from the residue
    # types, and the N terminus shares its residue number with a titratable side
    # chain, so both have to land under the same key.
    assert pkas[1] == [{"LYS": 11.34}, {"N+": 8.35}]
    assert pkas[130] == [{"C-": 3.22}]


def test_converted_pkas_survive_the_json_round_trip(tmp_path):
    """What convert_propka writes is what read_pka gives back.

    JSON has no integer keys, so the residue numbers go out as strings and have
    to come back as ints. read_pka is the only thing that converts them, and a
    structure looks up its residues by integer number, so a mapping keyed by
    strings would silently apply to nothing at all.
    """

    pkas = convert_propka(PROPKA_PATH)
    json_path = tmp_path / "1GDW_pka.json"
    write_json(pkas, str(json_path))

    assert read_pka(str(json_path)) == pkas


def test_redo_pkas_overrides_only_the_residues_it_is_given():
    """Predicted values replace the defaults, and untouched residues keep theirs."""

    structure = PDBparser().parse(PDB_PATH)
    residues = {residue.number: residue for residue in structure.residues}

    assert residues[18].pkas == [{"ASP": 3.86}], "default ASP pKa changed, update this test"
    assert residues[1].pkas == [{"LYS": 10.5}, {"N+": 9.69}]

    structure.redo_pkas(convert_propka(PROPKA_PATH))

    assert residues[18].pkas == [{"ASP": 3.92}]
    assert residues[1].pkas == [{"LYS": 11.34}, {"N+": 8.35}]

    # 45 of the 130 residues are in the prediction. The rest are not titratable
    # and must come through untouched rather than being blanked.
    assert residues[2].pkas is None


def test_redo_pkas_is_a_public_method_of_structure():
    """calculate() reaches the override through this name, so it cannot be renamed quietly."""

    assert callable(Structure.redo_pkas)


def test_a_pka_file_changes_the_calculated_features(tmp_path):
    """End to end: predicted pKas reach the output and move the charge features.

    The test is that the values are actually used, not what they come out as.
    Asserting the numbers themselves would pin PROPKA's predictions rather than
    prodes' handling of them, and would have to be rewritten every time the
    charge calculation is legitimately improved.
    """

    json_path = tmp_path / "1GDW_pka.json"
    write_json(convert_propka(PROPKA_PATH), str(json_path))

    default = pipeline_output(PDB_PATH, full_features=False)
    predicted = pipeline_output(PDB_PATH, full_features=False, pkas_file=str(json_path))

    assert list(default.columns) == list(predicted.columns)
    assert default["ID"].iloc[0] == predicted["ID"].iloc[0]

    numeric = default.select_dtypes("number").columns
    changed = [column for column in numeric if default[column].iloc[0] != predicted[column].iloc[0]]

    # Formal charge is the most direct consequence: a predicted pKa either side
    # of pH 7 flips whether a residue counts as charged at all.
    assert "Formal charge" in changed
    assert "Isoelectric point" in changed
    assert len(changed) > 10, f"only {len(changed)} features responded to the pKa file: {changed}"


# The pKa file and the structure can disagree about a residue, and until issue #6
# the file always won by default: redo_pkas replaced a residue's whole list, so a
# group the file did not mention was deleted rather than left at its default.


def test_a_group_the_file_omits_keeps_its_default():
    """A file naming a side chain but not the terminus must not delete the terminus.

    PROPKA always emits both, so this never bit on a PROPKA file, but the
    deletion was silent and the terminal atom then took its charge from the
    side-chain value.
    """

    structure = PDBparser().parse(PDB_PATH)
    first = structure.residues[0]

    structure.redo_pkas({1: [{"LYS": 11.34}]})

    assert first.pkas == [{"LYS": 11.34}, {"N+": 9.69}]
    assert first.group_pka("N+") == 9.69


def test_omitted_groups_are_reported(caplog):
    """Silence used to mean "use the textbook value", which is rarely what was meant."""

    structure = PDBparser().parse(PDB_PATH)

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas({1: [{"LYS": 11.34}]})

    assert "are not in the pKa file" in caplog.text


def test_a_complete_pka_file_is_reported_silently(caplog):
    """PROPKA covers all 46 titratable groups of 1GDW, so nothing should be flagged."""

    structure = PDBparser().parse(PDB_PATH)

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas(convert_propka(PROPKA_PATH))

    assert caplog.text == ""


def test_a_predicted_pka_cannot_titrate_a_bonded_cysteine():
    """The structure outranks the file: the group the file predicts does not exist."""

    structure = PDBparser().parse(PDB_PATH)
    cysteine = next(residue for residue in structure.residues if residue.number == 6)
    assert cysteine.disulfide_partner is not None, "1GDW cysteine 6 is bonded to 128"

    structure.redo_pkas({6: [{"CYS": 9.0}]})

    assert cysteine.pkas is None
    assert round(cysteine.charge(11), 3) == 0


def test_propka_reports_bonded_cysteines_as_not_titratable():
    """PROPKA writes 99.99 for a bridged cysteine rather than leaving it out.

    The issue that prompted this expected an omission, and the safety net was
    designed around that. It is worth pinning what the file really contains, since
    the two mechanisms have to agree rather than fight.
    """

    pkas = convert_propka(PROPKA_PATH)
    bonded = [6, 30, 65, 77, 81, 95, 116, 128]

    assert [pkas[number] for number in bonded] == [[{"CYS": PROPKA_NOT_TITRATABLE}]] * len(bonded)


def test_a_pka_for_a_group_the_residue_does_not_have_is_dropped(caplog):
    """Such an entry used to reach charged_atoms and raise a TypeError there."""

    structure = PDBparser().parse(PDB_PATH)
    second = structure.residues[1]

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas({second.number: [{"SER": 10.0}]})

    assert "which that residue does not have" in caplog.text
    assert second.pkas is None
    assert second.charge(7) == 0
