"""
Try to identify which asvs can be ignored for a PR.
"""
import os
import shlex
import subprocess

here = os.path.dirname(__file__)  # assumes repo directory
bmark_dir = os.path.join("asv_bench", "benchmarks")


def find_changed():
    """
    Find the names of files changes (compared to master) in this branch.
    """
    cmd = "git diff upstream/master --name-only"

    p = subprocess.Popen(
        shlex.split(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        close_fds=True,
    )
    stdout, stderr = p.communicate()

    flist = stdout.decode("utf-8").splitlines()
    return flist


def find_asv_paths():
    """
    Find all of our asv benchmark files.
    """
    bmarks = os.path.join(here, bmark_dir)
    walk = os.scandir(bmarks)
    walk = [
        x
        for x in walk
        if x.is_file()
        and x.name.endswith(".py")
        and x.name not in ["__init__.py", "pandas_vb_common.py"]
    ]

    paths = [x.path.replace(bmarks, "").lstrip(os.sep) for x in walk]
    return paths


def trim_asv_paths():
    """
    Exclude benchmark files we can be confident are not affected
    by this PR.
    """
    flist = find_changed()

    changed = [x for x in flist if "pandas/" in x and "pandas/tests/" not in x]

    asv_paths = find_asv_paths()

    if "setup.py" in flist:
        # We can't exclude anything
        return asv_paths
    if any("_libs/src" in x for x in changed):
        # We can't exclude anything
        return asv_paths
    if not any("_libs/tslibs" in x for x in changed):
        # If nothing in tslibs changed, then none of the tslibs asvs should
        #  be affected
        asv_paths = [x for x in asv_paths if not x.startswith("tslibs/")]
    if not any("_libs/" in x for x in changed):
        asv_paths = [
            x for x in asv_paths if x not in ["libs.py", "indexing_engines.py"]
        ]
        # TODO:
        #  inference.MaybeConvertNumeric
        #  dtypes.InferDtypes
        #  algorithms.MaybeConvertObjects  # (not 100% disjoint)
    return asv_paths


def path_to_module(path: str) -> str:
    """
    >>> path = "asv_bench/benchmarks/ctors.py"
    >>> path_to_module(path)
    'ctors'

    >>> path = "asv_bench/benchmarks/tslibs/fields.py"
    >>> path_to_module(path)
    'tslibs.fields'
    """
    assert path.startswith(bmark_dir)
    assert path.endswith(".py")

    name = path.replace(bmark_dir, "").strip(os.sep)
    name = name[:-3]
    name = name.replace(os.sep, ".")
    return name


def run_relevant_asvs():
    paths = trim_asv_paths()

    names = [path_to_module(x) for x in paths]

    pattern = "-b ".join(names)
    run_asvs(pattern)


def run_asvs(pattern: str):

    cmd = (
        "asv continuous "
        "-E virtualenv "
        "-f 1.05 "
        "--record-samples --append-samples "
        f"master HEAD -b {pattern}"
    )

    os.chdir("asv_bench")
    try:
        p = subprocess.Popen(
            shlex.split(cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
        )
        stdout, stderr = p.communicate()
    finally:
        os.chdir(here)

    stdout = stdout.decode("utf-8")
    stderr = stderr.decode("utf-8")
    return stdout, stderr
