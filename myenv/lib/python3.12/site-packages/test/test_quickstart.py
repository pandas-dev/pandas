# Licensed under a 3-clause BSD style license - see LICENSE.rst

from io import StringIO
from pathlib import Path

from asv import util

from . import tools


def test_quickstart(tmpdir, monkeypatch):
    dest = Path(tmpdir, "separate")
    dest.mkdir()

    tools.run_asv("quickstart", "--no-top-level", "--dest", str(dest))

    asv_conf_path = dest / "asv.conf.json"
    assert asv_conf_path.exists()
    benchmark_path = dest / "benchmarks" / "benchmarks.py"
    assert benchmark_path.exists()

    conf = util.load_json(str(asv_conf_path), js_comments=True)
    assert "env_dir" not in conf
    assert "html_dir" not in conf
    assert "results_dir" not in conf

    dest = Path(tmpdir, "same")
    dest.mkdir()

    monkeypatch.setattr("sys.stdin", StringIO("1"))
    tools.run_asv("quickstart", "--dest", str(dest))

    asv_conf_path = dest / "asv.conf.json"
    assert asv_conf_path.exists()
    benchmark_path = dest / "benchmarks" / "benchmarks.py"
    assert benchmark_path.exists()

    conf = util.load_json(str(asv_conf_path), js_comments=True)
    assert conf["env_dir"] != "env"
    assert conf["html_dir"] != "html"
    assert conf["results_dir"] != "results"
