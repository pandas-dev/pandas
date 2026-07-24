"""tests of pyodide-lock customization"""

from __future__ import annotations

import os
import abc
import json
import textwrap
import shutil

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

from .conftest import WHEELS

if TYPE_CHECKING:
    from collections.abc import Callable

    from pytest_console_scripts import ScriptRunner, RunResult

    TConfig = dict[str, dict[str, Any]]
    TLockRunner = Callable[[list[str], int], RunResult]

#: JSON-compatible encoding
UTF8 = {"encoding": "utf-8"}

#: configuration path
JLC_JSON = "jupyter_lite_config.json"

#: a valid timestamp for IPython 9.11.0
IPY911_EPOCH = 1_772_814_229

#: a set of de-normalized specs for IPython 9.11.0 to override the pyodide defaults
IPY911_SPECS = [
    "IPython ==9.11.0",
    "jedi >=0.18.2",
    "matplotlib_inline >=0.1.6",
    "Pygments >=2.14.0",
]

#: A URL for a federated extension/kernel package
BQ_URL = "https://pypi.org/packages/py2.py3/b/bqplot/bqplot-0.13.1-py2.py3-none-any.whl"

#: a set of configurations, keyed use as fixture values; `lock_enabled` set to True
CONFIGS: dict[str, TConfig] = dict(
    defaults={
        # enabled, no other configuration
        "PyodideAddon": {},
    },
    wheel={
        # a single local wheel with no dependencies
        "PyodideLockAddon": {"wheels": [f"{WHEELS[0]}"]},
    },
    wheels={
        # a single folder containing a wheel
        "PyodideLockAddon": {"wheels": [f"{WHEELS[0].parent}"]},
    },
    wheel_well_known={
        # a wheel in `lite_dir/static/pyodide-lock`
        "PyodideLockAddon": {},
    },
    ipy911_specs={
        # specs for a specific IPython, not in any Pyodide distribution
        "PyodideLockAddon": {"specs": IPY911_SPECS},
        "LiteBuildConfig": {"source_date_epoch": IPY911_EPOCH},
    },
    ipy911_constraints={
        # constraints for a specific version of IPython, not in any Pyodide distribution
        "PyodideLockAddon": {
            "constraints": IPY911_SPECS,
            "lite_constraints_file": "lite-constraints.txt",
        },
    },
    widgets={
        # ipywidgets
        "PyodideLockAddon": {"specs": ["ipywidgets"]},
    },
    all_remote={
        # use any publicly-available URL instead of locally-cached copies
        "PyodideLockAddon": {
            "pyodide_lock_uv_options": {"preserve_url_prefixes": ["https://"]},
        }
    },
    fed_ext={
        # a federated extension wheel that should constrain the solved version
        "LiteBuildConfig": {"federated_extensions": [BQ_URL, f"{WHEELS[0]}"]},
        "PyodideLockAddon": {"specs": ["bqplot", "the-smallest-extension"]},
    },
)

#: keys of CONFIGS that should not use the test fixture tarball
CONFIG_NO_TARBALL: set[str] = {"all_remote"}

#: keys of CONFIGS that should get a wheel in the well-known folder
CONFIG_ADD_WELL_KNOWN: set[str] = {"wheel_well_known"}

#: the PyPI dependencies of ``pyodide_kernel`` not included in the current Pyodide distribution
DEFAULT_WHEELS = {"comm", "ipython_pygments_lexers"}

#: additional wheels expected for a specific version of IPython
IPY911_WHEELS = {"ipython", "pygments", "jedi", "matplotlib_inline"}

#: the wheels that come with widgets
WIDGETS_WHEELS = {"ipywidgets", "jupyterlab_widgets", "widgetsnbextension"}

#: wheels to expect in the ``static/pyodide-lock`` folder after a build
CONFIG_EXPECT_WHEEL_STEMS: dict[str, set[str]] = dict(
    defaults=DEFAULT_WHEELS,
    wheel={"the_smallest_extension", *DEFAULT_WHEELS},
    wheels={"the_smallest_extension", *DEFAULT_WHEELS},
    wheel_well_known={"the_smallest_extension", *DEFAULT_WHEELS},
    ipy911_specs={*IPY911_WHEELS, *DEFAULT_WHEELS},
    ipy911_constraints={*IPY911_WHEELS, *DEFAULT_WHEELS},
    widgets={*WIDGETS_WHEELS, *DEFAULT_WHEELS},
    all_remote=set(),
    fed_ext={
        "bqplot",
        "bqscales",
        "the_smallest_extension",
        "traittypes",
        *WIDGETS_WHEELS,
        *DEFAULT_WHEELS,
    },
)

#: hide expected warnings about favicon, translation, etc.
CLI_ENV = {
    **os.environ,
    "JUPYTERLITE_NO_JUPYTER_SERVER": "true",
    "JUPYTERLITE_NO_JUPYTERLAB_SERVER": "true",
}


#: extra checks to perform
CONFIG_POST: dict[str, Callable[[], list[type[PostCheck]]]] = dict(
    ipy911_specs=lambda: [CheckReqFiles, CheckBreakLock],
    ipy911_constraints=lambda: [CheckLiteConstraints],
    fed_ext=lambda: [CheckFederated],
)


def test_pyodide_lock(
    some_post_run_checks: list[type[PostCheck]],
    an_empty_lite_dir: Path,
    run_with_lock: TLockRunner,
    a_lock_config: str,
    the_pyodide_lock_version: str,
) -> None:
    """can we generate a custom Pyodide lockfile?"""
    run_with_lock(["status"], 0)
    run_with_lock(["build"], 0)
    run_with_lock(["check"], 0)

    for post in some_post_run_checks:
        post(
            an_empty_lite_dir, a_lock_config, the_pyodide_lock_version, run_with_lock
        ).check()


# fixtures
@pytest.fixture(params=sorted(CONFIGS))
def a_lock_config(request: pytest.FixtureRequest, the_pyodide_lock_version: str) -> str:
    """Provide a key from ``CONFIGS`` above."""
    return f"{request.param}"


@pytest.fixture
def some_post_run_checks(a_lock_config: str) -> list[type[PostCheck]]:
    get_custom = CONFIG_POST.get(a_lock_config)
    return [CheckPaths, *(get_custom() if get_custom else [])]


@pytest.fixture
def run_with_lock(
    a_lock_config: str,
    an_empty_lite_dir: Path,
    a_pyodide_tarball: str,
    script_runner: ScriptRunner,
) -> TLockRunner:
    """Provide a pre-configured script runner."""
    conf = deepcopy(CONFIGS[a_lock_config])
    pyodide_config = conf.setdefault("PyodideAddon", {})
    lock_config = conf.setdefault("PyodideLockAddon", {})

    pyodide_config.update(
        pyodide_url=(
            None if a_lock_config in CONFIG_NO_TARBALL else f"{a_pyodide_tarball}"
        ),
    )
    lock_config.update(enabled=True)
    (an_empty_lite_dir / JLC_JSON).write_text(json.dumps(conf), **UTF8)

    if a_lock_config in CONFIG_ADD_WELL_KNOWN:
        well_known = an_empty_lite_dir / "static/pyodide-lock"
        well_known.mkdir(parents=True)
        shutil.copy2(WHEELS[0], well_known / WHEELS[0].name)

    def run(args: list[str], expect_rc: int = 0) -> RunResult:
        res = script_runner.run(
            ["jupyter", "lite", *args],
            cwd=an_empty_lite_dir,
            env=CLI_ENV,
        )
        rc = res.returncode
        ok = rc == expect_rc
        if not ok:  # pragma: no cover
            cache_dir = an_empty_lite_dir / ".cache/pyodide-lock"
            paths = map(str, sorted(cache_dir.rglob("*.*")))
            print("\n".join(["cached paths", *paths]))
        assert ok, f"did not return {expect_rc} from {args}"
        return res

    return run


# post-run checks
@dataclass
class PostCheck:
    an_empty_lite_dir: Path
    a_lock_config: str
    the_pyodide_lock_version: str
    run: TLockRunner

    @property
    def out(self) -> Path:
        return self.an_empty_lite_dir / "_output"

    @property
    def out_pyodide(self) -> Path:
        return self.out / "static/pyodide/pyodide.mjs"

    @property
    def out_lock(self) -> Path:
        return self.out / "static/pyodide-lock/pyodide-lock.json"

    @abc.abstractmethod
    def check(self) -> None:
        """Subclasses implement this and add their own ``assert`."""

    def clean_out(self) -> None:
        shutil.rmtree(self.out)


class CheckPaths(PostCheck):
    """Check whether key paths exist."""

    def check(self) -> None:
        if self.a_lock_config in CONFIG_NO_TARBALL:
            assert not self.out_pyodide.exists()
        else:
            assert self.out_pyodide.exists()
        assert self.out_lock.exists()

        expect_wheels = CONFIG_EXPECT_WHEEL_STEMS[self.a_lock_config]
        wheels = sorted(self.out_lock.parent.glob("*.whl"))
        wheel_names = {w.name.split("-")[0] for w in wheels}
        extra_wheels = wheel_names - expect_wheels
        missing_wheels = expect_wheels - wheel_names
        assert not missing_wheels, "not enough wheels after build"
        assert not extra_wheels, "too many wheels after build"


class CheckBreakLock(PostCheck):
    """Break a bunch of things."""

    def check(self) -> None:
        lock = json.loads(self.out_lock.read_text(**UTF8))
        lock["packages"].pop("jedi")
        next(self.out_lock.parent.glob("matplotlib*.whl")).unlink()
        self.out_lock.write_text(json.dumps(lock, indent=2, sort_keys=True), **UTF8)
        res = self.run(["doit", "--", "-s", "check"], 1)
        assert "[ipython] missing dependency: jedi" in res.stderr
        assert "[matplotlib-inline] missing wheel" in res.stderr
        self.out_lock.unlink()
        res = self.run(["doit", "--", "-s", "check"], 1)
        lock.pop("info")
        self.out_lock.write_text(json.dumps(lock, indent=2, sort_keys=True), **UTF8)
        res = self.run(["doit", "--", "-s", "check"], 1)
        msg = f"Failed to parse lock with pyodide-lock v{self.the_pyodide_lock_version}"
        assert msg in res.stderr


class CheckFederated(PostCheck):
    """Check whether a federated extension constrains the solve."""

    def check(self) -> None:
        for path in sorted(
            (self.an_empty_lite_dir / ".cache/pyodide-lock/_work").glob("*")
        ):
            print(path.name)
            print(textwrap.indent(path.read_text(**UTF8), "\t"))
        wheels = sorted(w.name for w in self.out_lock.parent.glob("*.whl"))
        assert Path(BQ_URL).name in wheels, (
            "federated extension did not constrain solve"
        )


class CheckReqFiles(PostCheck):
    """Re-run with specs in a nested ``requirements.txt`` file."""

    def check(self) -> None:
        lite_dir = self.an_empty_lite_dir
        jlc_json = lite_dir / JLC_JSON
        reqs = lite_dir / "reqs"
        reqs_txt = reqs / "requirements.txt"
        reqs_txt2 = reqs / "real/requirements.txt"
        old_wheels = {p.name for p in self.out_lock.parent.glob("*.whl")}

        jlc = json.loads(jlc_json.read_text(**UTF8))
        pla = jlc["PyodideLockAddon"]
        real_specs = pla.pop("specs")

        self.clean_out()

        spec_lines = [
            "# a comment",
            "# a misleading comment about -r requirements.txt",
            "-r ./real/requirements.txt # another comment",
            "--requirements ../reqs/real/requirements.txt # a second copy",
            "",
            "# the above line is intentionally blank",
        ]
        spec_lines2 = [
            *real_specs,
            "-r ../requirements.txt  # intentional circular reference",
        ]
        reqs_txt2.parent.mkdir(parents=True)
        reqs_txt2.write_text("\n".join(spec_lines2), **UTF8)
        reqs_txt.write_text("\n".join(spec_lines), **UTF8)

        pla["specs"] = [f"-r reqs/{reqs_txt.name}"]
        jlc_json_text = json.dumps(jlc, indent=2, sort_keys=True)
        print("building with -r", jlc_json_text)
        jlc_json.write_text(jlc_json_text, **UTF8)

        self.run(["build"], 0)
        req_wheels = {p.name for p in self.out_lock.parent.glob("*.whl")}
        assert old_wheels == req_wheels

        self.clean_out()
        shutil.rmtree(reqs)

        ppt = lite_dir / "pyproject.toml"
        ppt_lines = [
            "[dependency-groups]",
            f"demo = {json.dumps(real_specs)}",
            "demo2 = [{include-group='demo'}]",
        ]
        ppt.write_text("\n".join(ppt_lines))

        pla["specs"] = ["-g demo", f"--group={ppt.name}:demo2"]
        jlc_json_text = json.dumps(jlc, indent=2, sort_keys=True)
        print("building with -g", jlc_json_text)
        jlc_json.write_text(jlc_json_text, **UTF8)

        self.run(["build"], 0)
        ppt_wheels = {p.name for p in self.out_lock.parent.glob("*.whl")}
        assert old_wheels == ppt_wheels


class CheckLiteConstraints(PostCheck):
    """Check that the lite constraints file was written."""

    def check(self) -> None:
        lcp = self.an_empty_lite_dir / "lite-constraints.txt"
        old_text = lcp.read_text(**UTF8)
        old_lines = old_text.splitlines()
        old_wheels = {*self.out_lock.parent.glob("*.whl")}
        non_empty_lines = [line for line in old_lines if line.split("#")[0].strip()]
        assert non_empty_lines
        bad_names = ("ipykernel", "pyodide-kernel", "piplite")
        bad_lines = {
            line for line in old_lines if line.split(" ")[0].strip() in bad_names
        }
        assert not bad_lines, "unexpected bad names"
        self.clean_out()
        self.run(["build"], 0)
        new_text = lcp.read_text(**UTF8)
        assert new_text == old_text, "the constraints should not change"
        new_wheels = {*self.out_lock.parent.glob("*.whl")}
        assert new_wheels == old_wheels, "unexpected wheels"
