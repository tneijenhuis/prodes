"""Tests for prodes.feature_dictionary.FeatureDictionary.

Ordered by how often the thing being tested is used:

  1. Everyday use    - getting the code lists, translating a code to a label
  2. Staying in sync - the YAML dictionaries vs what calculate() really writes
  3. Details         - field completeness, edge cases, error messages

The YAML files are the single source of truth: the code lists in
prodes.feature_dictionary are derived from them, not written out a second time.
So there is nothing to keep in sync *within* the package, and the only sync worth
testing is against the actual output of calculate().
"""

import importlib.resources

import pandas as pd
import pytest
import yaml

from prodes.feature_dictionary import (
    DROPPED_FEATURE_CODES,
    FULL_FEATURE_CODES,
    FULL_ONLY_DICTIONARY_FILE,
    ID_COLUMN,
    REDUCED_DICTIONARY_FILE,
    REDUCED_FEATURE_CODES,
    FeatureDictionary,
)
from tests.pipeline_output import pipeline_output

# =============================================================================
# 1. Everyday use
# =============================================================================


def test_get_the_reduced_feature_codes():
    """The most common call: what does Prodes write by default?"""
    fd = FeatureDictionary()
    codes = fd.get_feature_codes()

    assert len(codes) == 54
    assert codes[0] == "Molecular weight"
    assert "SurfEpMeanFormal" in codes
    assert "SurfEpMeanAverage" not in codes


def test_get_the_full_feature_codes():
    """The other common call: everything Prodes can write."""
    fd = FeatureDictionary()
    codes = fd.get_feature_codes(full=True)

    assert len(codes) == 105
    assert "SurfEpMeanAverage" in codes


def test_translate_a_feature_code_to_a_plot_label():
    """Labelling an axis from a column heading."""
    fd = FeatureDictionary()

    assert fd.get_plot_name("SurfEpMeanFormal") == "Surface EP mean (formal)"
    assert fd.get_plot_name("Molecular weight") == "Molecular weight (Da)"


def test_translate_a_feature_code_to_a_description():
    """Building a title or caption from a column heading."""
    fd = FeatureDictionary()

    description = fd.get_description("SurfEpMeanFormal")
    assert description == "Mean EP across all surface points, formal charges."
    assert f"Impact of {description} on R2".startswith("Impact of Mean EP")


def test_translate_a_feature_code_to_a_long_description():
    fd = FeatureDictionary()
    long_description = fd.get_long_description("SurfEpMeanFormal")

    assert len(long_description) > len(fd.get_description("SurfEpMeanFormal"))
    assert "surface" in long_description.lower()


def test_ask_why_a_feature_is_not_in_the_default_set():
    fd = FeatureDictionary()

    assert fd.get_reason_dropped("SurfEpMeanFormal") is None
    assert "SurfEpMeanFormal" in fd.get_reason_dropped("SurfEpMeanAverage")


def test_label_a_whole_dataframe_of_results():
    """The realistic plotting loop: every reduced column gets a label."""
    fd = FeatureDictionary()
    labels = {code: fd.get_plot_name(code) for code in fd.get_feature_codes()}

    assert len(labels) == 54
    assert all(labels.values())


def test_line_up_a_full_run_with_a_reduced_run():
    """Mixed datasets: cut the full ones down so the two concatenate cleanly."""
    fd = FeatureDictionary()
    combined = pd.concat([fd.drop_redundant_features(FULL_OUTPUT), REDUCED_OUTPUT], ignore_index=True)

    assert list(combined.columns) == [ID_COLUMN] + fd.get_feature_codes()
    assert not combined.isna().any().any()


# =============================================================================
# 2. Staying in sync with what calculate() actually writes
# =============================================================================

PDB_PATH = "tests/data/ARH96693.pdb.zip"
ORIG_OUTPUT = "tests/data/ARH96693_prodes_orig_output.csv"

# Shared with test_sasa.py, which needs the same two runs of the same structure.
FULL_OUTPUT = pipeline_output(PDB_PATH, full_features=True)
REDUCED_OUTPUT = pipeline_output(PDB_PATH, full_features=False)


def load_yaml(filename):
    with importlib.resources.open_text("prodes.data", filename, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


REDUCED_DOC = load_yaml(REDUCED_DICTIONARY_FILE)
FULL_ONLY_DOC = load_yaml(FULL_ONLY_DICTIONARY_FILE)


def test_full_codes_match_a_real_run():
    """Compared as sets: the full list is reduced-then-dropped, not column order.

    Column order is covered by the regression test in test_sasa.py, which compares
    a full run against the committed original-code reference column by column.
    """
    assert set(FULL_FEATURE_CODES) == {c for c in FULL_OUTPUT.columns if c != ID_COLUMN}
    assert len(FULL_FEATURE_CODES) == len(FULL_OUTPUT.columns) - 1


def test_reduced_codes_match_a_real_run():
    """Order matters here: drop_redundant_features builds frames in this order."""
    assert list(REDUCED_FEATURE_CODES) == [c for c in REDUCED_OUTPUT.columns if c != ID_COLUMN]


def test_full_codes_match_the_original_reference_output():
    """The full feature set must be exactly the columns the original Prodes wrote.

    tests/data/ARH96693_prodes_orig_output.csv was calculated with the original,
    pre-fork version, so it is the reference for what the full feature set is.

    This test is *designed* to fail as soon as anyone adds, removes or renames a
    feature. That is not a nuisance, it is the prompt: a new feature needs a
    decision about whether it belongs in the reduced default set, an entry in the
    YAML dictionaries with a plot name and explanations, and a check of whether it
    is redundant with something already there (see
    docs/redundant_feature_analysis.md, which records how the existing ones were
    checked). Editing the YAML to make this test pass, without making those
    decisions, defeats the point.
    """
    reference = {column for column in pd.read_csv(ORIG_OUTPUT).columns if column != ID_COLUMN}
    ours = set(FULL_FEATURE_CODES)

    added = sorted(ours - reference)
    missing = sorted(reference - ours)

    assert not added, f"in the YAML dictionaries but not in the original Prodes output: {added}. If this is a genuinely new feature, see this test's docstring."
    assert (
        not missing
    ), f"in the original Prodes output but missing from the YAML dictionaries: {missing}. A feature cannot simply be deleted from the dictionaries; it is still written under --full-features."


def test_the_two_yaml_files_partition_the_full_set():
    """Between them the files must cover every feature exactly once.

    A code appearing in both would be swallowed silently: the module merges them
    with ``{**reduced, **full_only}``, so the full-only entry would win and the
    code would also be counted twice in FULL_FEATURE_CODES.
    """
    reduced, full_only = set(REDUCED_DOC["features"]), set(FULL_ONLY_DOC["features"])

    assert reduced & full_only == set(), "a feature code is defined in both files"
    assert reduced | full_only == set(FULL_FEATURE_CODES)
    assert len(FULL_FEATURE_CODES) == len(set(FULL_FEATURE_CODES)), "duplicate code in FULL_FEATURE_CODES"


def test_declared_yaml_counts_match_the_contents():
    assert REDUCED_DOC["meta"]["feature_count"] == len(REDUCED_DOC["features"]) == len(REDUCED_FEATURE_CODES)
    assert FULL_ONLY_DOC["meta"]["feature_count"] == len(FULL_ONLY_DOC["features"]) == len(DROPPED_FEATURE_CODES)


def test_drop_redundant_features_matches_a_real_reduced_run():
    """Reducing a full run must give the same values as calculating reduced directly."""
    reduced_from_full = FeatureDictionary().drop_redundant_features(FULL_OUTPUT)
    pd.testing.assert_frame_equal(reduced_from_full.reset_index(drop=True), REDUCED_OUTPUT, check_dtype=False)


# =============================================================================
# 3. Details
# =============================================================================

DESCRIPTION_FIELDS = ("plot_name", "short_explanation", "long_explanation")
DROP_FIELDS = ("dropped_because", "redundant_with", "drop_note")


@pytest.mark.parametrize("doc_name", ["reduced", "full_only"])
def test_every_entry_has_all_description_fields(doc_name):
    doc = REDUCED_DOC if doc_name == "reduced" else FULL_ONLY_DOC
    for code, entry in doc["features"].items():
        for field in DESCRIPTION_FIELDS:
            assert entry.get(field), f"{code} has no {field}"


def test_full_only_entries_explain_the_drop():
    for code, entry in FULL_ONLY_DOC["features"].items():
        for field in DROP_FIELDS:
            assert entry.get(field), f"{code} has no {field}"
        assert isinstance(entry["redundant_with"], list) and entry["redundant_with"], code


def test_reduced_entries_carry_no_drop_fields():
    """A kept feature was never dropped, so it must not claim a reason."""
    for code, entry in REDUCED_DOC["features"].items():
        for field in DROP_FIELDS:
            assert field not in entry, f"{code} is kept but has {field}"


def test_redundant_with_names_are_real_feature_codes():
    for code, entry in FULL_ONLY_DOC["features"].items():
        for partner in entry["redundant_with"]:
            assert partner in set(FULL_FEATURE_CODES), f"{code} cites unknown feature {partner}"


def test_drop_note_names_the_feature_it_describes():
    """The note must name its own feature, so it cannot be read as describing a neighbour."""
    for code, entry in FULL_ONLY_DOC["features"].items():
        note = entry["drop_note"]
        assert code in note or any(p in note for p in entry["redundant_with"]), f"{code} drop_note names neither itself nor its partner: {note}"


def test_sum_features_do_not_quote_a_pairwise_r_squared():
    """Sum drops rest on an exact identity; a weak pairwise R2 would read as counter-evidence."""
    for code, entry in FULL_ONLY_DOC["features"].items():
        if entry["dropped_because"] == "sum_equals_mean_times_count":
            assert "r_squared_with_first_partner" not in entry, code


def test_getters_cover_every_feature_code():
    fd = FeatureDictionary()
    for code in FULL_FEATURE_CODES:
        assert fd.get_plot_name(code)
        assert fd.get_description(code)
        assert fd.get_long_description(code)


def test_reason_dropped_is_none_for_kept_and_set_for_dropped():
    fd = FeatureDictionary()
    for code in REDUCED_FEATURE_CODES:
        assert fd.get_reason_dropped(code) is None, code
    for code in DROPPED_FEATURE_CODES:
        assert fd.get_reason_dropped(code), code


def test_getters_reject_unknown_codes():
    fd = FeatureDictionary()
    for getter in (fd.get_plot_name, fd.get_description, fd.get_long_description, fd.get_reason_dropped, fd.get_entry):
        with pytest.raises(KeyError, match="not a Prodes feature code"):
            getter("NotAFeature")


def test_id_is_not_a_feature_code():
    assert ID_COLUMN not in FULL_FEATURE_CODES
    assert ID_COLUMN not in REDUCED_FEATURE_CODES


def test_plot_names_are_unique():
    """Two features sharing an axis label would make a plot ambiguous."""
    fd = FeatureDictionary()
    labels = [fd.get_plot_name(code) for code in FULL_FEATURE_CODES]
    duplicates = {label for label in labels if labels.count(label) > 1}
    assert duplicates == set(), f"duplicate plot names: {duplicates}"


def test_returned_lists_and_dicts_are_copies():
    """Callers editing a result must not corrupt the shared dictionary."""
    fd = FeatureDictionary()

    codes = fd.get_feature_codes()
    codes.append("junk")
    assert "junk" not in fd.get_feature_codes()

    entry = fd.get_entry("Area")
    entry["plot_name"] = "junk"
    assert fd.get_plot_name("Area") != "junk"

    everything = fd.get_dictionary()
    everything["Area"]["plot_name"] = "junk"
    assert fd.get_plot_name("Area") != "junk"


def test_instances_share_the_loaded_files():
    """Constructing repeatedly must not re-read the YAML each time."""
    assert FeatureDictionary()._entries is FeatureDictionary()._entries


def test_get_dictionary_covers_the_full_set():
    assert set(FeatureDictionary().get_dictionary()) == set(FULL_FEATURE_CODES)


def test_membership_test():
    fd = FeatureDictionary()
    assert "SurfEpMeanFormal" in fd
    assert "NotAFeature" not in fd


def test_drop_redundant_features_without_id():
    reduced = FeatureDictionary().drop_redundant_features(FULL_OUTPUT, keep_id=False)
    assert list(reduced.columns) == list(REDUCED_FEATURE_CODES)


def test_drop_redundant_features_is_idempotent():
    fd = FeatureDictionary()
    once = fd.drop_redundant_features(FULL_OUTPUT)
    pd.testing.assert_frame_equal(once, fd.drop_redundant_features(once))


def test_drop_redundant_features_reports_missing_columns():
    with pytest.raises(KeyError, match="Area"):
        FeatureDictionary().drop_redundant_features(FULL_OUTPUT.drop(columns=["Area"]))


@pytest.mark.parametrize("doc_name", ["reduced", "full_only"])
def test_unit_and_original_explanation_keys_are_always_present(doc_name):
    """Both keys exist on every entry, holding null when there is nothing to say.

    A missing key reads as an oversight; an explicit null says the feature is
    genuinely dimensionless, or that the paper does not document it.
    """
    doc = REDUCED_DOC if doc_name == "reduced" else FULL_ONLY_DOC
    for code, entry in doc["features"].items():
        assert "unit" in entry, f"{code} has no unit key"
        assert "original_explanation" in entry, f"{code} has no original_explanation key"


def test_get_unit_where_the_paper_gives_one():
    """Units come from Supplemental Table 1 of the original paper."""
    fd = FeatureDictionary()

    assert fd.get_unit("Molecular weight") == "Da"
    assert fd.get_unit("Area") == "Å²"
    assert fd.get_unit("Dipole") == "Å"
    assert fd.get_unit("SurfEpMeanFormal") == "v"


def test_get_unit_is_none_when_dimensionless_or_undocumented():
    fd = FeatureDictionary()

    assert fd.get_unit("ALASurfFrac") is None  # a fraction, dimensionless
    assert fd.get_unit("Shape max") is None
    assert fd.get_unit("SurfPosMhpMedian") is None  # not in the paper at all


def test_every_unit_is_one_of_the_documented_ones():
    fd = FeatureDictionary()
    units = {fd.get_unit(code) for code in FULL_FEATURE_CODES}

    assert units <= {None, "Da", "v", "Å", "Å²"}


def test_get_original_explanation_returns_the_paper_wording():
    fd = FeatureDictionary()

    assert fd.get_original_explanation("Molecular weight") == "Sum of the weight of each amino acid in the protein"
    assert fd.get_original_explanation("SurfMhpSum") == "The sum of all hydrophobicity potentials"


def test_original_explanation_is_none_for_features_the_paper_omits():
    """Two Pos/Neg MHP medians were added to the code after the paper was published."""
    fd = FeatureDictionary()

    assert fd.get_original_explanation("SurfPosMhpMedian") is None
    assert fd.get_original_explanation("SurfNegMhpMedian") is None


def test_the_paper_documents_all_but_two_features():
    fd = FeatureDictionary()
    undocumented = [code for code in FULL_FEATURE_CODES if fd.get_original_explanation(code) is None]

    assert undocumented == ["SurfPosMhpMedian", "SurfNegMhpMedian"]


def test_our_description_and_the_paper_agree_on_the_statistic():
    """Guards against a description being attached to the wrong feature.

    The original paper included feature descriptions in the supplementary material,
    which we have included as original_explanation in the features yaml files.
    Keywords in the original descriptions should match the current short_explanation.
    E.g. if our text says "median" where the paper says "mean", one of the two write-ups
    describes a different feature than the column actually holds.
    """
    import re

    fd = FeatureDictionary()
    # "mean" and "average" are the same statistic, and the paper prefers one word
    # while our descriptions prefer the other. But "average" also names the
    # partial-charge model, so those phrases are stripped first or every *Average
    # feature would look as though it claimed to be a mean.
    charge_model = r"average charges?|partial-charge model|charges ranging between 1 and -1"
    statistics = {
        "trimean": r"trimean",
        "median": r"median",
        "mean": r"\bmean\b|\baverage\b",
        "std": r"standard deviation|\bstd\b",
        "max": r"maximum|\bmax\b",
        "count": r"number of|\bcount\b",
    }

    disagreements = []
    for code in FULL_FEATURE_CODES:
        paper = fd.get_original_explanation(code)
        if not paper:
            continue

        paper = re.sub(charge_model, " ", paper, flags=re.I)
        ours = re.sub(charge_model, " ", fd.get_description(code), flags=re.I)
        for name, pattern in statistics.items():
            in_paper = bool(re.search(pattern, paper, re.I))
            in_ours = bool(re.search(pattern, ours, re.I))
            if in_paper != in_ours:
                disagreements.append(f"{code}: '{name}' in paper={in_paper}, ours={in_ours}")

    assert disagreements == [], "descriptions disagree on which statistic the feature is: " + "; ".join(disagreements)


def test_molecular_weight_is_not_described_as_atom_based():
    """Structure.mw sums residue masses from a table; it never reads atom records.

    The source dictionary originally claimed the opposite, which would have told
    users the value changes with whether hydrogens are in the PDB. It does not.
    """
    fd = FeatureDictionary()
    text = fd.get_long_description("Molecular weight").lower()

    assert "residue" in text
    assert "does not change" in text or "unaffected" in text
