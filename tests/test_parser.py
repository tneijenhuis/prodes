import pytest

from prodes.io import parser

file_path = "tests/data/1GDW.pdb.zip"
pdb_parser = parser.PDBparser()
pdb_parser.identifier = "ATOM"


def test_main_fail():
    """tests if a incorect file will give an error"""

    with pytest.raises(ValueError):
        pdb_parser.parse(file_path[:-2])


def test_main_accept():
    """tests if a file is parsed correctly"""

    global structure
    structure = pdb_parser.parse(file_path)

    assert len(structure.atoms) == 1022
    assert len(structure.chains) == 1
    assert len(structure.residues) == 130


def test_read(tmp_path):
    """tests if a pdb file read"""

    extracted = parser.extract_pdb(file_path, tmp_path)
    with open(extracted) as file:
        pdb_parser._read_pdb(file, "test")


def test_structure_name_comes_from_the_file_stem():
    """the structure name, which becomes the output ID, is the bare file name

    Guards the Windows case: a backslash path used to leave the whole path in
    the name, which matters because zipped structures are parsed from a
    temporary directory rather than from tests/data.
    """

    assert pdb_parser.parse(file_path).name == "1GDW"


def test_structure_name_keeps_every_extension_but_the_last(tmp_path):
    """a name with more than one dot keeps all of it, since only the suffix is dropped

    1abc.ent.pdb gives 1abc.ent, not 1abc. Pinned because the name becomes the
    output ID, and because splitting on the first dot, which is what the code
    used to do, gave the other answer.
    """

    extracted = parser.extract_pdb(file_path, tmp_path)
    multi_dot = tmp_path / "1abc.ent.pdb"
    multi_dot.write_text(open(extracted).read())

    assert pdb_parser.parse(str(multi_dot)).name == "1abc.ent"


def test_zipped_structure_is_named_after_the_archive(tmp_path):
    """an archive names the structure, not the member file inside it

    A user pointing at bar.pdb.zip expects to get a row called bar, whatever the
    file inside happens to be called.
    """

    import zipfile

    extracted = parser.extract_pdb(file_path, tmp_path)
    archive = tmp_path / "renamed.pdb.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.write(extracted, arcname="something_else.pdb")

    assert pdb_parser.parse(str(archive)).name == "renamed"


def test_zipped_and_plain_parse_identically(tmp_path):
    """a zipped structure parses to the same atoms as the extracted plain file"""

    extracted = parser.extract_pdb(file_path, tmp_path)

    from_archive = pdb_parser.parse(file_path)
    from_plain = pdb_parser.parse(extracted)

    assert len(from_archive.atoms) == len(from_plain.atoms)
    assert [(a.name, a.x, a.y, a.z) for a in from_archive.atoms] == [(a.name, a.x, a.y, a.z) for a in from_plain.atoms]


def test_write_pdb_round_trips(tmp_path):
    """a structure written back out parses to the same atoms it came from

    write_pdb had no test and, once the old scratch scripts were removed, no
    caller either. This keeps its column layout honest: the parser reads by
    fixed column position, so any drift in the writer shows up here.
    """

    structure = pdb_parser.parse(file_path)
    written = tmp_path / "round_trip.pdb"

    parser.write_pdb(structure, str(written))
    reparsed = pdb_parser.parse(str(written))

    assert len(reparsed.atoms) == len(structure.atoms)
    assert [atom.name for atom in reparsed.atoms] == [atom.name for atom in structure.atoms]
    assert [atom.chain_name for atom in reparsed.atoms] == [atom.chain_name for atom in structure.atoms]
    assert [atom.element for atom in reparsed.atoms] == [atom.element for atom in structure.atoms]


def test_write_pdb_rejects_a_non_pdb_name(tmp_path, capsys):
    """asking for a file without a pdb extension writes nothing and says so"""

    target = tmp_path / "structure.txt"

    parser.write_pdb(pdb_parser.parse(file_path), str(target))

    assert "can only make files" in capsys.readouterr().out
    assert not target.exists()


def test_write_pdb_rejects_residue_numbers_that_do_not_fit(tmp_path):
    """a 5 digit residue number is refused rather than silently truncated

    Columns 23-26 are all the format gives the residue number, so a wider one
    would overflow into the chain field and read back as a different residue.
    """

    structure = pdb_parser.parse(file_path)
    structure.atoms[0].residue_number = 12345
    target = tmp_path / "overflow.pdb"

    with pytest.raises(ValueError, match="more than the 4 columns"):
        parser.write_pdb(structure, str(target))

    assert not target.exists()


def test_builder_makes_a_placeholder_atom():
    """the Builder produces a dummy atom at the requested position"""

    atom = parser.Builder().build_dummy_atom(1.5, -2.5, 3.5, chain_name="B")

    assert (atom.x, atom.y, atom.z) == (1.5, -2.5, 3.5)
    assert atom.chain_name == "B"
    assert atom.element == "X"


def test_archive_with_no_pdb_is_rejected(tmp_path):
    """an archive that does not hold exactly one pdb file raises a clear error"""

    import zipfile

    empty = tmp_path / "empty.pdb.zip"
    with zipfile.ZipFile(empty, "w") as archive:
        archive.writestr("readme.txt", "no structure here")

    with pytest.raises(ValueError, match="expected exactly one .pdb file"):
        pdb_parser.parse(str(empty))


def test_parsed_atom():
    """tests parsed atom"""

    atom = structure.atoms[0]

    assert atom.name == "N"
    assert atom.identifier == "ATOM"
    assert atom.chain_name == "A"
    assert atom.residue_number == 1
    assert atom.residue_name == "LYS"
    assert atom.x == 1.134
    assert atom.y == 19.824
    assert atom.z == 22.575
    assert atom.structure == structure


def test_parsed_chain():
    """tests parsed chain"""

    chain = structure.chains[0]
    assert chain.name == "A"


def test_parsed_residue():
    """tests parsed residue"""

    residue = structure.residues[0]

    assert residue.name == "LYS"
    assert residue.number == 1
    assert len(residue.atoms) == 9
