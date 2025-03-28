"""integration tests for overall CLI functionality"""

import platform
import time

from pytest import mark

from jupyterlite_core import __version__
from jupyterlite_core.constants import HOOKS

PY_IMPL = platform.python_implementation()
IS_PYPY = "pypy" in PY_IMPL.lower()
IS_WIN = platform.system() == "Windows"


# TODO: others?
LITE_INVOCATIONS = [
    ["jupyter-lite"],
    ["jupyter", "lite"],
    ["python", "-m", "jupyterlite_core"],
]

# nothing we can do about this, at present
TRASH = [".jupyterlite.doit.db"]

# some files we expect to exist after a full build
A_GOOD_BUILD = [
    *TRASH,
    "_output/package.json",
    "_output/jupyter-lite.json",
    "_output/index.html",
]

# these hooks will not generate a build
FAST_HOOKS = ["list", "status"]

# serve is handled separately
NOT_SERVE_HOOK = [h for h in HOOKS if h != "serve"]

# a simple overrides.json
AN_OVERRIDES = """{
  "@jupyterlab/apputils-plugin:themes": {
    "themeScrollbars": true
  }
}
"""

# a simple jupyter-lite.json describing a remote entry
A_SIMPLE_JUPYTERLITE_JSON = """{ "jupyter-config-data": {
    "federated_extensions": [
        {
            "extension": "./extension",
            "load": "static/remoteEntry.abc123.js",
            "name": "@org/pkg"
        }
    ],
    "disabledExtensions": ["@org/pkg"],
    "settingsOverrides": {}
} }"""


@mark.parametrize("lite_args", LITE_INVOCATIONS)
def test_cli_version(lite_args, script_runner):
    """do various invocations work"""
    returned_version = script_runner.run([*lite_args, "--version"])
    assert returned_version.success
    assert __version__ in returned_version.stdout
    assert returned_version.stderr == ""


@mark.parametrize("lite_args", LITE_INVOCATIONS)
@mark.parametrize("help", ["-h", "--help"])
def test_cli_help(lite_args, help, script_runner):  # noqa: A002
    """does help work"""
    returned_version = script_runner.run([*lite_args, help])
    assert returned_version.success
    assert returned_version.stderr == ""


@mark.parametrize("lite_args", LITE_INVOCATIONS)
def test_nonzero_rc(lite_args, script_runner):
    a_step = script_runner.run([*lite_args, "doit", "this-is-not-a-step"])
    assert not a_step.success


@mark.parametrize("lite_hook", ["list", "status"])
def test_cli_status_null(lite_hook, an_empty_lite_dir, script_runner):
    """do the "side-effect-free" commands create exactly one file?"""
    returned_status = script_runner.run(["jupyter", "lite", lite_hook], cwd=str(an_empty_lite_dir))
    assert returned_status.success
    files = set(an_empty_lite_dir.rglob("*"))
    # we would expect to see our build cruft sqlite
    assert len(files) == 1
    dododb = an_empty_lite_dir / ".jupyterlite.doit.db"
    assert files == {dododb}


@mark.parametrize("lite_hook", NOT_SERVE_HOOK)
def test_cli_any_hook(  # noqa: PLR0915
    lite_hook, an_empty_lite_dir, script_runner, a_simple_lite_ipynb
):
    """does all the hooks basically work

    TODO: this should be broken up into a hypothesis state machine, perhaps
    """
    expected_files = TRASH if lite_hook in FAST_HOOKS else A_GOOD_BUILD
    started = time.time()
    returned_status = script_runner.run(["jupyter", "lite", lite_hook], cwd=str(an_empty_lite_dir))
    duration_1 = time.time() - started
    assert returned_status.success
    files = set(an_empty_lite_dir.rglob("*"))
    assert len(files) >= 1
    for expected in expected_files:
        assert (an_empty_lite_dir / expected).exists()

    if lite_hook in FAST_HOOKS:
        return

    # re-run, be faster
    restarted = time.time()
    rereturned_status = script_runner.run(
        ["jupyter", "lite", lite_hook],
        cwd=str(an_empty_lite_dir),
    )
    duration_2 = time.time() - restarted
    assert rereturned_status.success

    if not (IS_PYPY or IS_WIN):
        # some caching doesn't seep to work reliably
        assert duration_1 > duration_2

    # force, detect a root file with ``=`` in the name
    readme = an_empty_lite_dir / "== README ==.md"
    readme.write_text("# hello world", encoding="utf-8")

    # ... and a nested folder
    more = an_empty_lite_dir / "more"
    details = more / "details"
    details_readme = details / "README.md"
    details_readme.parent.mkdir(parents=True)
    details_readme.write_text("# more details", encoding="utf-8")

    # some federated stuff
    lite_json = an_empty_lite_dir / "jupyter-lite.json"
    lite_json.write_text(A_SIMPLE_JUPYTERLITE_JSON, encoding="utf-8")

    lite_ipynb = an_empty_lite_dir / "jupyter-lite.ipynb"
    lite_ipynb.write_text(a_simple_lite_ipynb, encoding="utf-8")

    # ... and app overrides
    app_overrides = an_empty_lite_dir / "lab/overrides.json"
    app_overrides.parent.mkdir()
    app_overrides.write_text(AN_OVERRIDES, encoding="utf-8")

    forced_status = script_runner.run(
        [
            "jupyter",
            "lite",
            lite_hook,
            "--force",
            "--contents",
            str(readme),
            "--contents",
            str(more),
        ],
        cwd=str(an_empty_lite_dir),
    )

    if A_GOOD_BUILD[-1] in expected_files and lite_hook not in ["init"]:
        out = an_empty_lite_dir / "_output"

        # did the files make it...
        expected_readme = out / "files/== README ==.md"
        assert expected_readme.exists()
        assert "world" in expected_readme.read_text(encoding="utf-8")
        expected_details = out / "files/details/README.md"
        assert expected_details.exists()
        assert "details" in expected_details.read_text(encoding="utf-8")

        # ...and get indexed
        missed = 0
        for path in ["", "details"]:
            contents = (out / f"api/contents/{path}/all.json").read_text(encoding="utf-8")
            print("contents of", path, contents)
            if "README" not in contents:  # pragma: no cover
                missed += 1
        assert not missed, "some contents were not indexed"

        # default translation files should also be created
        all_packs_file = out / "api/translations/all.json"
        assert all_packs_file.exists()
        all_packs = all_packs_file.read_text(encoding="utf-8")
        assert "English" in all_packs

        en_pack_file = out / "api/translations/en.json"
        assert en_pack_file.exists()

    assert forced_status.success


def test_cli_raw_doit(an_empty_lite_dir, script_runner):
    """does raw doit work"""
    returned_status = script_runner.run(
        ["jupyter", "lite", "doit", "--", "--help"], cwd=str(an_empty_lite_dir)
    )
    assert returned_status.success
    assert "http://pydoit.org" in returned_status.stdout


def test_build_repl_no_sourcemaps(an_empty_lite_dir, script_runner):
    """does (re-)building create a predictable pattern of file counts"""
    out = an_empty_lite_dir / "_output"

    # disable the federated_extensions addon as it may output settings for federated extensions
    args = original_args = "jupyter", "lite", "build", "--disable-addons", "federated_extensions"
    status = script_runner.run(args, cwd=str(an_empty_lite_dir))
    norm_files = sorted(out.rglob("*"))
    assert status.success
    assert [f for f in norm_files if f.name.endswith(".map")], "expected maps"

    args = [*args, "--apps", "repl", "--apps", "foobarbaz"]
    status = script_runner.run(args, cwd=str(an_empty_lite_dir))
    repl_files = sorted(out.rglob("*"))
    repl_bundles = sorted(out.glob("build/*/bundle.js"))
    assert status.success

    assert len(repl_files) < len(norm_files), "expected fewer files"
    assert len(repl_bundles) == 1, "only expected one bundle"
    assert "'foobarbaz' is not one of" in status.stderr

    args = [*args, "--no-unused-shared-packages"]
    status = script_runner.run(args, cwd=str(an_empty_lite_dir))
    no_chunk_files = sorted(out.rglob("*"))
    # assert "pruning unused shared package" in status.stderr

    unexpected = sorted(set(map(str, no_chunk_files)) - set(map(str, repl_files)))

    assert len(no_chunk_files) < len(repl_files), f"unexpected {unexpected}"

    args = [*args, "--no-sourcemaps"]
    status = script_runner.run(args, cwd=str(an_empty_lite_dir))
    min_files = sorted(out.rglob("*"))
    assert status.success

    assert not [f for f in min_files if f.name.endswith(".map")], "expected no maps"

    assert len(min_files) < len(no_chunk_files), "expected fewer files still"

    status = script_runner.run(original_args, cwd=str(an_empty_lite_dir))
    rebuild_files = sorted(out.rglob("*"))
    assert status.success

    assert len(norm_files) == len(rebuild_files), "expected the same files"
