import pandas as pd
import pandas.diagnostics as diag


def test_merge_report_attaches_and_has_phases():
    left = pd.DataFrame({"k": [1, 2, 3], "v": [10, 20, 30]})
    right = pd.DataFrame({"k": [1, 2, 4], "w": [100, 200, 400]})

    out = diag.report(lambda: left.merge(right, on="k", how="left"), memory=True)
    assert out.shape == (3, 3)

    rep = getattr(out, "_diagnostics_report", None)
    assert rep is not None

    phase_names = [p.name for p in rep.phases]
    assert "merge.get_join_info" in phase_names
    assert "merge.reindex_and_concat" in phase_names

    assert rep.counters.get("join_len") == 3
