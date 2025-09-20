"""feature tests of generated artifacts"""

import pprint
import shutil
import tarfile
import tempfile
from hashlib import sha256
from pathlib import Path

# use the generally-documented invocation
LITE_ARGS = "jupyter", "lite"


def test_archive_input(a_lite_app_archive, the_npm_source_date_epoch):
    """validate some assumptions about the structure of our baseline tarball"""
    expected = dict(gid=0, uid=0, mtime=the_npm_source_date_epoch)
    unexpected = []

    with tarfile.open(a_lite_app_archive) as tar:
        for member in tar.getmembers():
            for expected_key, expected_value in expected.items():
                observed_value = getattr(member, expected_key)
                if observed_value == expected_value:
                    continue
                unexpected += [  # pragma: no cover
                    [expected_key, expected_value, observed_value, member.name]
                ]

    pprint.pprint(unexpected)
    assert not unexpected, f"{a_lite_app_archive} does not work the way we expect"


def test_archive_is_reproducible(an_empty_lite_dir, script_runner, source_date_epoch):
    """do we build reproducible artifacts?"""
    # TODO: handle macro-scale reproducibility in CI
    extra_args = "--debug", "--source-date-epoch", source_date_epoch
    archive_args = (*LITE_ARGS, "archive", *extra_args)
    cwd = dict(cwd=str(an_empty_lite_dir))

    # put some content in
    readme = an_empty_lite_dir / "files/nested/README.md"
    readme.parent.mkdir(parents=True)
    readme.write_text("# Hello world\n", encoding="utf-8")

    # build once for initial tarball
    before = an_empty_lite_dir / "v1.tgz"
    initial = script_runner.run(
        [*archive_args, "--output-archive", str(before), "--no-libarchive"], **cwd
    )
    assert initial.success, "failed to build the first tarball"

    # reset
    _reset_a_lite_dir(an_empty_lite_dir, before, readme)

    # build another tarball
    after = an_empty_lite_dir / "v2.tgz"
    subsequent = script_runner.run([*archive_args, "--output-archive", str(after)], **cwd)
    assert subsequent.success, "failed to build the second tarball"

    # check them
    _assert_same_tarball("two successive builds should be the same", script_runner, before, after)


def test_archive_is_idempotent(an_empty_lite_dir, script_runner, source_date_epoch):
    extra_args = "--debug", "--source-date-epoch", source_date_epoch
    archive_args = (*LITE_ARGS, "archive", *extra_args)
    cwd = dict(cwd=str(an_empty_lite_dir))

    # build once for initial tarball
    before = an_empty_lite_dir / "v1.tgz"
    initial = script_runner.run([*archive_args, "--output-archive", str(before)], **cwd)
    assert initial.success, "failed to build the first tarball"

    # build another tarball
    after = an_empty_lite_dir / "v2.tgz"
    subsequent = script_runner.run(
        [*archive_args, "--app-archive", str(before), "--output-archive", str(after)],
        **cwd,
    )
    assert subsequent.success, "failed to build the second tarball"

    # check them
    _assert_same_tarball("a build repeated should be the same", script_runner, before, after)


def _reset_a_lite_dir(lite_dir, *skip):
    """clean out a lite dir, except for the named files"""
    for path in lite_dir.glob("*"):
        if path in skip or path.is_dir():
            continue
        else:
            path.unlink()


def _assert_same_tarball(message, script_runner, before, after):  # pragma: no cover
    """helper function to compare two tarballs.

    TODO: the `diffoscope` HTML output is _definitely_ good enough for end users
          but is linux only...
    """
    fails = []

    hashes = {p: sha256(p.read_bytes()).hexdigest() for p in [before, after]}

    if len(set(hashes.values())) > 1:
        fails += [hashes]

    if fails:
        print("FAILS...", fails, flush=True)

    diffoscope = shutil.which("diffoscope")

    if diffoscope:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            print(tdp)

        diffoscope_result = script_runner.run([diffoscope, str(before), str(after)])

        if not diffoscope_result.success:
            fails += [diffoscope_result.stdout, diffoscope_result.stderr]
            print(
                f"{message} DIFFOSCOPE",
                diffoscope_result.stdout,
                diffoscope_result.stderr,
            )
    else:
        print("on unix platforms, install diffoscope for **MUCH** clearer output")

    assert not fails, f"tarballs were not the same: {message} {fails}"
