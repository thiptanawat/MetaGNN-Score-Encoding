"""Data-interface contracts with independent small examples, no model optimization."""
import copy
import csv
import json
from pathlib import Path
import numpy as np
import pytest
import xlrd
from src.gpr import evaluate
from src.prepare_data import make_chemistry_map, make_outcomes_and_noise, normalize_context, unique_normalized


def test_conservative_gpr_missing_cells_and_complete_alternatives():
    genes = {"A": np.array([1., 2.]), "B": np.array([np.nan, 3.]), "C": np.array([0., 8.])}
    np.testing.assert_array_equal(evaluate("A and B", genes), [np.nan, 2.])
    np.testing.assert_array_equal(evaluate("(A and B) or C", genes), [0., 8.])
    assert evaluate("A and missing", genes) is None
    np.testing.assert_array_equal(evaluate("(A and missing) or C", genes), [0., 8.])
    np.testing.assert_array_equal(evaluate("A or B and C", genes), [1., 3.])
    np.testing.assert_array_equal(evaluate("(A or B) and C", genes), [0., 3.])


@pytest.mark.parametrize("rule", ["A B", "(A and B", "A or", "A )", "and A", "()"])
def test_malformed_gpr_never_silently_truncates(rule):
    with pytest.raises(ValueError):
        evaluate(rule, {"A": [1], "B": [2]})


def test_context_aliases_and_collision_failure():
    assert normalize_context("786-O") == normalize_context("786-0")
    assert normalize_context("MDA-MB-231/ATCC") == normalize_context("MDA_MB_231")
    assert normalize_context("HL-60(TB)") == normalize_context("HL-60/TB")
    assert normalize_context("SNB-19") != normalize_context("U251")  # group origins, never merge profiles
    with pytest.raises(ValueError):
        unique_normalized(["786-O", "786-0"])


def chemical_fixture():
    item = {"metabolite_id": "m", "source_excel_row": 2, "source_label": "alanine", "method": "HILIC", "calibrated": True, "historical_variable": False, "historical_sd_ratio": 0.1}
    table = {("alanine", "HILIC"): {"source_label": "alanine", "source_kegg": "C00041", "source_pdf_page": 22}}
    model = {"metabolites": [{"id": "ala_e", "name": "L-Alanine", "formula": "C3H7NO2", "compartment": "e", "annotation": {"kegg.compound": ["C00041"]}}],
             "reactions": [{"id": "EX_ala_e", "metabolites": {"ala_e": -1}}]}
    return item, table, model


def test_chemistry_is_independent_of_full_cohort_variable_and_direction(tmp_path):
    (tmp_path / "manifests").mkdir()
    item, table, model = chemical_fixture()
    first = make_chemistry_map(tmp_path, [item], table, model)[0]
    assert first["primary_eligible"] is True and first["sign_factor"] == 1
    item["historical_variable"] = True
    item["historical_sd_ratio"] = 1000
    second = make_chemistry_map(tmp_path, [item], table, model)[0]
    assert (first["primary_eligible"], first["exchange_id"]) == (second["primary_eligible"], second["exchange_id"])
    model["reactions"][0]["metabolites"]["ala_e"] = 2
    reversed_scaled = make_chemistry_map(tmp_path, [item], table, model)[0]
    assert reversed_scaled["sign_factor"] == -2  # positive physical release is -s*v


def test_ambiguous_exchange_and_duplicate_assay_are_excluded(tmp_path):
    (tmp_path / "manifests").mkdir()
    item, table, model = chemical_fixture()
    model["reactions"].append({"id": "EX_ala_second_e", "metabolites": {"ala_e": -1}})
    assert not make_chemistry_map(tmp_path, [item], table, model)[0]["primary_eligible"]
    model["reactions"].pop()
    second = copy.deepcopy(item); second["metabolite_id"] = "m2"
    assert not any(row["primary_eligible"] for row in make_chemistry_map(tmp_path, [item, second], table, model))


def test_noise_uses_only_development_replicates(tmp_path):
    (tmp_path / "data").mkdir()
    item, _, _ = chemical_fixture()
    contexts = [{"context_id": "D", "core_label": "d", "partition": "development"}, {"context_id": "T", "core_label": "t", "partition": "test"}]
    mapping = {"d": "D", "t": "T"}
    values = {"m": {"d": [2., 5.], "t": [-1., 1.]}}
    make_outcomes_and_noise(tmp_path, [item], values, contexts, mapping)
    before = (tmp_path / "data/assay_noise_development.tsv").read_bytes()
    values["m"]["t"] = [-1e100, 1e100]
    make_outcomes_and_noise(tmp_path, [item], values, contexts, mapping)
    assert before == (tmp_path / "data/assay_noise_development.tsv").read_bytes()
    row = next(csv.DictReader((tmp_path / "data/assay_noise_development.tsv").open(), delimiter="\t"))
    assert float(row["tolerance"]) == 3 and int(row["n_dev_pairs"]) == 1


def test_prepared_data_source_contracts():
    root = Path(__file__).resolve().parents[1]
    if not (root / "data/preparation_summary.json").exists():
        pytest.skip("Run preparation first for raw-source integration contracts")
    summary = json.loads((root / "data/preparation_summary.json").read_text())
    assert len(summary["development_contexts"]) == 11
    assert len(summary["test_contexts"]) == 47
    assert len(summary["test_origin_groups"]) == 44
    assert not summary["historical_variable_used_for_primary"]
    scores = np.load(root / "data/scores_conservative.npz")
    assert scores["A"].shape == (5468, 59) and np.isfinite(scores["A"]).all()
    assert len(set(scores["rxn"])) == 5468
    rows = list(csv.DictReader((root / "manifests/metabolite_map.tsv").open(), delimiter="\t"))
    panel = [row for row in rows if row["primary_eligible"] == "True"]
    assert len({row["exchange_id"] for row in panel}) == len(panel)
    assert all(row["exchange_id"].startswith("EX_") and row["model_metabolite_id"].endswith("_e") for row in panel)
    assert not any(row["source_label"] in {"creatinine", "3-OH-kynurenate", "homocystine", "glycerol_1", "glycerol_2"} for row in panel)


def test_every_raw_outcome_is_preserved_with_original_orientation():
    root = Path(__file__).resolve().parents[1]
    if not (root / "data/outcomes_long.tsv").exists():
        pytest.skip("Preparation output is not present")
    sheet = xlrd.open_workbook(root / "data/raw/core/NIHMS419088-supplement-Database_S1.xls").sheet_by_name("CORE data")
    headers = sheet.row_values(0)
    exported = list(csv.DictReader((root / "data/outcomes_long.tsv").open(), delimiter="\t"))
    assert len(exported) == 16800
    seen = set()
    for record in exported:
        source_row = int(record["source_excel_row"]) - 1
        columns = [col for col in range(3, sheet.ncols) if str(headers[col]).strip() == record["core_label"]]
        source_column = columns[int(record["replicate"]) - 1]
        assert float(record["value"]) == sheet.cell_value(source_row, source_column)
        assert record["unit"] == ("fmol/cell/h" if sheet.cell_value(source_row, 2) == 1 else "arbitrary")
        seen.add((source_row, source_column))
    assert len(seen) == 16800
