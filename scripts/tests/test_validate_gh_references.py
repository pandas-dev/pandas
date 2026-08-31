import pytest

from scripts import validate_gh_references
from scripts.validate_gh_references import (
    main,
    update_baseline,
)

MSG = "https://github.com/pandas-dev/pandas/issues/1234"


@pytest.mark.parametrize(
    "src",
    [
        "# GH#1234\n",
        "# GH 1234\n",
        "# GH-1234\n",
        "# GH1234\n",
        "# gh-1234\n",
        "# GH #1234\n",
        "# GH: 1234\n",
        "# GH: #1234\n",
        "# GH##1234\n",
        "# Gh 1234\n",
        "x = 1  # GH#1234 trailing prose\n",
        '"""See GH#1234."""\n',
    ],
)
def test_short_forms_are_rejected(capsys, src):
    ret = main(src, "t.py", {})
    out, _ = capsys.readouterr()
    assert ret == 1
    assert MSG in out


@pytest.mark.parametrize(
    "src",
    [
        # the form we want
        "# https://github.com/pandas-dev/pandas/issues/1234\n",
        "# https://github.com/pandas-dev/pandas/pull/1234\n",
        # another project's tracker
        "# https://github.com/apache/arrow/issues/49410\n",
        # identifiers, not references
        "gh_13141_expected = {}\n",
        "ghost1234 = 1\n",
        # csv fixture data, and numbers too long to be issue numbers
        'data = """GH,100102040,jkl,0205"""\n',
        # "gh" inside a word
        "# ran through 1234 rows\n",
        "# a high 1234 water mark\n",
        # too short to be an issue number
        "# GH 1\n",
    ],
)
def test_non_references_are_ignored(capsys, src):
    ret = main(src, "t.py", {})
    out, _ = capsys.readouterr()
    assert ret == 0
    assert out == ""


def test_baseline_grandfathers_existing_references(capsys):
    ret = main("# GH#1234\n", "t.py", {"t.py": {"1234"}})
    out, _ = capsys.readouterr()
    assert ret == 0
    assert out == ""


def test_baseline_is_per_file(capsys):
    ret = main("# GH#1234\n", "t.py", {"other.py": {"1234"}})
    out, _ = capsys.readouterr()
    assert ret == 1
    assert out == (f"t.py:1:3 found 'GH#1234'; reference issues as {MSG}\n")


def test_baseline_is_per_issue_number(capsys):
    ret = main("# GH#1234\n# GH#5678\n", "t.py", {"t.py": {"1234"}})
    out, _ = capsys.readouterr()
    assert ret == 1
    assert "issues/5678" in out
    assert "issues/1234" not in out


def test_multiple_references_on_one_line(capsys):
    ret = main("# GH#1234, GH#5678\n", "t.py", {})
    out, _ = capsys.readouterr()
    assert ret == 1
    assert len(out.splitlines()) == 2


def test_update_baseline_refuses_to_grow(capsys, monkeypatch, tmp_path):
    path = tmp_path / "baseline.txt"
    path.write_text("t.py 1234\n", encoding="utf-8")
    monkeypatch.setattr(
        validate_gh_references, "collect_baseline", lambda: {"t.py": {"1234", "5678"}}
    )

    ret = update_baseline(path)
    out, _ = capsys.readouterr()
    assert ret == 1
    assert "refusing to update" in out
    # the refusal leaves the baseline alone
    assert path.read_text(encoding="utf-8") == "t.py 1234\n"


def test_update_baseline_allow_growth(capsys, monkeypatch, tmp_path):
    path = tmp_path / "baseline.txt"
    path.write_text("t.py 1234\n", encoding="utf-8")
    monkeypatch.setattr(
        validate_gh_references, "collect_baseline", lambda: {"t.py": {"1234", "5678"}}
    )

    ret = update_baseline(path, allow_growth=True)
    out, _ = capsys.readouterr()
    assert ret == 0
    assert "1 -> 2 references" in out
    assert path.read_text(encoding="utf-8").endswith("t.py 1234 5678\n")


def test_update_baseline_shrinks_without_the_flag(capsys, monkeypatch, tmp_path):
    path = tmp_path / "baseline.txt"
    path.write_text("t.py 1234 5678\n", encoding="utf-8")
    monkeypatch.setattr(
        validate_gh_references, "collect_baseline", lambda: {"t.py": {"1234"}}
    )

    ret = update_baseline(path)
    out, _ = capsys.readouterr()
    assert ret == 0
    assert "2 -> 1 references" in out
    assert path.read_text(encoding="utf-8").endswith("t.py 1234\n")
