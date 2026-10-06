"""Tests für die Design-Erkennung. Das Beispiel ist der echte Aufbau von GSE12220."""
import pandas as pd
import pytest

from src import design as d
from src.metadata import MetadataError


def make_samples() -> pd.DataFrame:
    rows = []
    gsm = 307217
    for stress, source_tag, desc in (
            ("oxidative stress", "H2O2", "0.3mM of H2O2"),
            ("DNA damage stress", "MMS", "0.1% of MMS")):
        for t in (0, 30, 60, 100, 140, 180):
            title = f"{stress}, t={t}" + ("" if t == 0 else "min")
            source = (f"time 0 following {source_tag} treatment" if t == 0
                      else f"{t} min following {source_tag} treatment")
            rows.append({
                "geo_accession": f"GSM{gsm}", "title": title, "source_name_ch1": source,
                "organism_ch1": "Saccharomyces cerevisiae", "characteristics_1": "strain Y262",
                "treatment_protocol_ch1": "At time 0, cells were treated with 0.3mM of H2O2 or 0.1% of MMS.",
                "description": f"Gene expression following treatment with {desc}"})
            gsm += 1
    return pd.DataFrame(rows).set_index("geo_accession")


@pytest.mark.parametrize("text, expected", [
    ("30 min", 30.0), ("2 h", 120.0), ("1.5 hours", 90.0), ("t=30min", 30.0),
    ("t=0", 0.0), ("time 0 following H2O2 treatment", 0.0),
    ("0.3mM of H2O2", None), ("no time here", None),
])
def test_extract_time_minutes(text, expected):
    assert d.extract_time_minutes(text) == expected


def test_build_design_gse12220_structure():
    design = d.build_design(make_samples())
    assert design["treatment"].value_counts().to_dict() == {"H2O2": 6, "MMS": 6}
    assert design["is_baseline"].sum() == 2
    assert sorted(design["time_min"].unique()) == [0, 30, 60, 100, 140, 180]
    assert not design["needs_review"].any()
    assert design.loc["GSM307218", "concentration"] == "0.3 mM"
    assert design.loc["GSM307224", "concentration"] == "0.1 %"
    assert (design["strain"] == "Y262").all()
    assert (design["replicate"] == 1).all()  # keine Replikate


def test_protocol_text_naming_both_treatments_is_ignored():
    # Würde treatment_protocol_ch1 ausgewertet, wäre jedes Sample 'unklar'.
    design = d.build_design(make_samples())
    assert "unknown" not in set(design["treatment"])


def test_ambiguous_sample_needs_review():
    samples = make_samples()
    samples.loc["GSM307218", "source_name_ch1"] = "30 min following H2O2 and MMS"
    design = d.build_design(samples)
    assert design.loc["GSM307218", "needs_review"]


def test_conflicting_times_not_guessed():
    samples = make_samples()
    samples.loc["GSM307218", "title"] = "oxidative stress, t=60min"  # Quelle sagt 30
    design = d.build_design(samples)
    assert pd.isna(design.loc["GSM307218", "time_min"])
    assert design.loc["GSM307218", "needs_review"]


def test_apply_overrides():
    samples = make_samples()
    samples.loc["GSM307218", "source_name_ch1"] = "unklar"
    samples.loc["GSM307218", "description"] = "unklar"
    design = d.build_design(samples)
    assert design.loc["GSM307218", "treatment"] == "unknown"
    fixed = d.apply_overrides(design, pd.DataFrame(
        {"sample_id": ["GSM307218"], "treatment": ["H2O2"], "time_min": [30.0]}))
    assert fixed.loc["GSM307218", "treatment"] == "H2O2"
    assert not fixed["needs_review"].any()


def test_confirm_saves_and_reuses(tmp_path):
    design = d.build_design(make_samples())
    d.confirm_design(design, tmp_path, "GSE12220", input_func=lambda _: "j")
    assert (tmp_path / "GSE12220_design_confirmed.csv").exists()
    # zweiter Aufruf darf nicht mehr fragen
    again = d.confirm_design(design, tmp_path, "GSE12220",
                             input_func=lambda _: pytest.fail("Es sollte nicht gefragt werden"))
    assert list(again.index) == list(design.index)
    assert again["is_baseline"].sum() == 2
    assert (again["strain"] == "Y262").all()


def test_reject_writes_override_template(tmp_path):
    design = d.build_design(make_samples())
    with pytest.raises(d.DesignNotConfirmedError):
        d.confirm_design(design, tmp_path, "X", input_func=lambda _: "n")
    assert (tmp_path / "X_design_overrides.csv").exists()


def test_unclear_design_refuses_even_with_yes(tmp_path):
    samples = make_samples()
    samples.loc["GSM307218", "source_name_ch1"] = "30 min following H2O2 and MMS"
    design = d.build_design(samples)
    with pytest.raises(d.DesignNotConfirmedError):
        d.confirm_design(design, tmp_path, "X", assume_yes=True)


def test_non_interactive_without_yes_refuses(tmp_path):
    with pytest.raises(d.DesignNotConfirmedError):
        d.confirm_design(d.build_design(make_samples()), tmp_path, "X", interactive=False)


def test_unknown_override_sample_raises():
    design = d.build_design(make_samples())
    with pytest.raises(MetadataError):
        d.apply_overrides(design, pd.DataFrame({"sample_id": ["GSM1"], "treatment": ["H2O2"]}))
