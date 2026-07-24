"""
Helper functions for testing.
"""
import itertools
import locale
import logging
import os
from pathlib import Path
import string
import subprocess
import sys
from tempfile import TemporaryDirectory

import matplotlib as mpl
from matplotlib import _api

_log = logging.getLogger(__name__)


def set_font_settings_for_testing():
    mpl.rcParams['font.family'] = 'DejaVu Sans'
    mpl.rcParams['text.hinting'] = 'default'


def set_reproducibility_for_testing():
    mpl.rcParams['svg.hashsalt'] = 'matplotlib'


def setup():
    # The baseline images are created in this locale, so we should use
    # it during all of the tests.

    try:
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    except locale.Error:
        try:
            locale.setlocale(locale.LC_ALL, 'English_United States.1252')
        except locale.Error:
            _log.warning(
                "Could not set locale to English/United States. "
                "Some date-related tests may fail.")

    mpl.use('Agg')

    with _api.suppress_matplotlib_deprecation_warning():
        mpl.rcdefaults()  # Start with all defaults

    # These settings *must* be hardcoded for running the comparison tests and
    # are not necessarily the default values as specified in rcsetup.py.
    set_font_settings_for_testing()
    set_reproducibility_for_testing()


def subprocess_run_for_testing(command, env=None, timeout=60, stdout=None,
                               stderr=None, check=False, text=True,
                               capture_output=False, **kwargs):
    """
    Create and run a subprocess.

    Thin wrapper around `subprocess.run`, intended for testing.  Will
    mark fork() failures on Cygwin as expected failures: not a
    success, but not indicating a problem with the code either.

    Parameters
    ----------
    args : list of str
    env : dict[str, str]
    timeout : float
    stdout, stderr
    check : bool
    text : bool
        Also called ``universal_newlines`` in subprocess.  I chose this
        name since the main effect is returning bytes (`False`) vs. str
        (`True`), though it also tries to normalize newlines across
        platforms.
    capture_output : bool
        Set stdout and stderr to subprocess.PIPE

    Returns
    -------
    proc : subprocess.Popen

    See Also
    --------
    subprocess.run

    Raises
    ------
    pytest.skip
        If running on emscripten, which does not support subprocesses.
    pytest.xfail
        If platform is Cygwin and subprocess reports a fork() failure.
    """
    if sys.platform == 'emscripten':
        import pytest
        pytest.skip('emscripten does not support subprocesses')
    if capture_output:
        stdout = stderr = subprocess.PIPE
    # Add CREATE_NO_WINDOW flag on Windows to prevent console window overhead
    # This is added in an attempt to fix flaky timeouts of subprocesses on Windows
    if sys.platform == 'win32':
        if 'creationflags' not in kwargs:
            kwargs['creationflags'] = subprocess.CREATE_NO_WINDOW
        else:
            kwargs['creationflags'] |= subprocess.CREATE_NO_WINDOW
    try:
        proc = subprocess.run(
            command, env=env,
            timeout=timeout, check=check,
            stdout=stdout, stderr=stderr,
            text=text, **kwargs
        )
    except BlockingIOError:
        if sys.platform == "cygwin":
            # Might want to make this more specific
            import pytest
            pytest.xfail("Fork failure")
        raise
    except subprocess.CalledProcessError as e:
        if e.stdout:
            _log.error(f"Subprocess output:\n{e.stdout}")
        if e.stderr:
            _log.error(f"Subprocess error:\n{e.stderr}")
        raise e
    if proc.stdout:
        _log.debug(f"Subprocess output:\n{proc.stdout}")
    if proc.stderr:
        _log.debug(f"Subprocess error:\n{proc.stderr}")
    return proc


def subprocess_run_helper(func, *args, timeout, extra_env=None):
    """
    Run a function in a sub-process.

    Parameters
    ----------
    func : function
        The function to be run.  It must be in a module that is importable.
    *args : str
        Any additional command line arguments to be passed in
        the first argument to ``subprocess.run``.
    extra_env : dict[str, str]
        Any additional environment variables to be set for the subprocess.
    """
    target = func.__name__
    module = func.__module__
    file = func.__code__.co_filename
    proc = subprocess_run_for_testing(
        [
            sys.executable,
            "-c",
            f"import importlib.util;"
            f"_spec = importlib.util.spec_from_file_location({module!r}, {file!r});"
            f"_module = importlib.util.module_from_spec(_spec);"
            f"_spec.loader.exec_module(_module);"
            f"_module.{target}()",
            *args,
        ],
        env={
            **os.environ,
            "SOURCE_DATE_EPOCH": "0",
            # subprocess_run_helper sets SOURCE_DATE_EPOCH=0 above, so for a dirty tree,
            # the version will have the date 19700101 which breaks pickle tests with a
            # warning if the working tree is dirty.
            #
            # This will also avoid at least one additional subprocess call for
            # setuptools-scm query git, so we tell the subprocess what version
            # to report as the test process.
            "SETUPTOOLS_SCM_PRETEND_VERSION_FOR_MATPLOTLIB": mpl.__version__,
            **(extra_env or {}),
        },
        timeout=timeout,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return proc


def _check_for_pgf(texsystem):
    """
    Check if a given TeX system + pgf is available

    Parameters
    ----------
    texsystem : str
        The executable name to check
    """
    with TemporaryDirectory() as tmpdir:
        tex_path = Path(tmpdir, "test.tex")
        tex_path.write_text(r"""
            \documentclass{article}
            \usepackage{pgf}
            \begin{document}
            \typeout{pgfversion=\pgfversion}
            \makeatletter
            \@@end
        """, encoding="utf-8")
        try:
            subprocess.check_call(
                [texsystem, "-halt-on-error", "-no-shell-escape",
                 str(tex_path)], cwd=tmpdir,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            return False
        return True


def _has_tex_package(package):
    try:
        mpl.dviread.find_tex_file(f"{package}.sty")
        return True
    except (FileNotFoundError, OSError):
        return False


def ipython_in_subprocess(requested_backend_or_gui_framework, all_expected_backends):
    import pytest
    IPython = pytest.importorskip("IPython")

    if sys.platform == "win32":
        pytest.skip("Cannot change backend running IPython in subprocess on Windows")

    if (IPython.version_info[:3] == (8, 24, 0) and
            requested_backend_or_gui_framework == "osx"):
        pytest.skip("Bug using macosx backend in IPython 8.24.0 fixed in 8.24.1")

    # This code can be removed when Python 3.12, the latest version supported
    # by IPython < 8.24, reaches end-of-life in late 2028.
    for min_version, backend in all_expected_backends.items():
        if IPython.version_info[:2] >= min_version:
            expected_backend = backend
            break

    code = ("import matplotlib as mpl, matplotlib.pyplot as plt;"
            "fig, ax=plt.subplots(); ax.plot([1, 3, 2]); mpl.get_backend()")
    proc = subprocess_run_for_testing(
        [
            "ipython",
            "--no-simple-prompt",
            f"--matplotlib={requested_backend_or_gui_framework}",
            "-c", code,
        ],
        check=True,
        capture_output=True,
    )

    assert proc.stdout.strip().endswith(f"'{expected_backend}'")


def is_ci_environment():
    # Common CI variables
    ci_environment_variables = [
        'CI',        # Generic CI environment variable
        'CONTINUOUS_INTEGRATION',  # Generic CI environment variable
        'TRAVIS',    # Travis CI
        'CIRCLECI',  # CircleCI
        'JENKINS',   # Jenkins
        'GITLAB_CI',  # GitLab CI
        'GITHUB_ACTIONS',  # GitHub Actions
        'TEAMCITY_VERSION'  # TeamCity
        # Add other CI environment variables as needed
    ]

    for env_var in ci_environment_variables:
        if os.getenv(env_var):
            return True

    return False


def _gen_multi_font_text():
    """
    Generate text intended for use with multiple fonts to exercise font fallbacks.

    Returns
    -------
    fonts : list of str
        The names of the fonts used to render the test string, sorted by intended
        priority. This should be set as the font family for the Figure or Text artist.
    text : str
        The test string.
    """
    # These fonts are serif and sans-serif, and would not normally be combined, but that
    # should make it easier to see which glyph is from which font.
    fonts = ['cmr10', 'DejaVu Sans']
    # cmr10 does not contain accented characters, so they should fall back to DejaVu
    # Sans. However, some accented capital A versions *are* in cmr10 with non-standard
    # glyph shapes, so don't test those (otherwise this Latin1 supplement group would
    # start at 0xA0.)
    start = 0xC5
    latin1_supplement = [chr(x) for x in range(start, 0xFF+1)]
    latin_extended_A = [chr(x) for x in range(0x100, 0x17F+1)]
    latin_extended_B = [chr(x) for x in range(0x180, 0x24F+1)]
    non_basic_multilingual_plane = [chr(x) for x in range(0x1F600, 0x1F610)]
    count = itertools.count(start - 0xA0)
    non_basic_characters = '\n'.join(
        ''.join(line)
        for _, line in itertools.groupby(  # Replace with itertools.batched for Py3.12+.
            [*latin1_supplement, *latin_extended_A, *latin_extended_B,
             *non_basic_multilingual_plane],
            key=lambda x: next(count) // 32)  # 32 characters per line.
    )
    test_str = f"""There are basic characters
{string.ascii_uppercase} {string.ascii_lowercase}
{string.digits} {string.punctuation}
and accented characters
{non_basic_characters}
in between!"""
    # The resulting string contains 491 unique characters. Some file formats use 8-bit
    # tables, which the large number of characters exercises twice over.
    return fonts, test_str
