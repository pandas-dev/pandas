import pandas as pd
from pandas._diagnostics import report


def test_report_returns_wrapper_and_report_object():
    left = pd.DataFrame({"k": [1, 2, 3], "v": [10, 20, 30]})
    right = pd.DataFrame({"k": [1, 2, 4], "w": [100, 200, 400]})

    out = report(left.merge, right, on="k", how="left", profile=False, memory="off")
    assert hasattr(out, "value")
    assert hasattr(out, "report")

    df = out.value
    rep = out.report

    assert df.shape == (3, 3)
    assert rep is not None
    assert rep.total_seconds >= 0.0
    assert isinstance(rep.metadata, dict)
    assert isinstance(rep.counters, dict)
    assert isinstance(rep.phases, list)


def test_merge_report_has_metadata_and_memory_fields_when_enabled():
    left = pd.DataFrame({"k": [1, 2, 3], "v": [10, 20, 30]})
    right = pd.DataFrame({"k": [1, 2, 4], "w": [100, 200, 400]})

    out = report(
        left.merge,
        right,
        on="k",
        how="left",
        memory="delta",
        allocations=True,
        profile=True,
        bucket_phases=True,
        include_diagnostics_phases=False,
    )
    rep = out.report

    assert rep.metadata.get("operation") in ("merge", "DataFrame.merge")
    assert rep.metadata.get("left_rows") == 3
    assert rep.metadata.get("right_rows") == 3
    assert rep.metadata.get("how") == "left"

    assert rep.memory_bytes_net is not None
    assert rep.memory_current_bytes is not None
    assert rep.memory_peak_bytes is not None

    assert isinstance(rep.auto_phases, list)
    assert len(rep.auto_phases) > 0

    assert all(not p.name.startswith("diag.") for p in rep.phases)


def test_include_diagnostics_phases_flag_controls_diag_phases():
    left = pd.DataFrame({"k": [1, 2, 3], "v": [10, 20, 30]})
    right = pd.DataFrame({"k": [1, 2, 4], "w": [100, 200, 400]})

    out_hidden = report(left.merge, right, on="k", how="left", memory="peak")
    names_hidden = [p.name for p in out_hidden.report.phases]
    assert all(not n.startswith("diag.") for n in names_hidden)

    out_shown = report(
        left.merge,
        right,
        on="k",
        how="left",
        memory="peak",
        include_diagnostics_phases=True,
    )
    names_shown = [p.name for p in out_shown.report.phases]
    assert any(n.startswith("diag.") for n in names_shown)


def test_merge_counters_present_when_merge_instrumentation_available():
    left = pd.DataFrame({"k": [1, 2, 3], "v": [10, 20, 30]})
    right = pd.DataFrame({"k": [1, 2, 4], "w": [100, 200, 400]})

    out = report(left.merge, right, on="k", how="left", profile=True, memory="peak")
    rep = out.report

    if "join_len" in rep.counters:
        assert rep.counters["join_len"] == 3


def test_report_to_text_smoke():
    left = pd.DataFrame({"k": [1, 2, 3], "v": [10, 20, 30]})
    right = pd.DataFrame({"k": [1, 2, 4], "w": [100, 200, 400]})

    out = report(left.merge, right, on="k", how="left", profile=True, memory="peak")
    txt = out.report.to_text()

    assert isinstance(txt, str)
    assert "pandas diagnostics report" in txt
