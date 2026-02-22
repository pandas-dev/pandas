from __future__ import annotations

import os
import sys
from typing import Literal
import uuid
import warnings

import pytest


_sys_path = list(sys.path)  # freeze a copy of sys.path at interpreter startup


@pytest.hookimpl
def pytest_xdist_auto_num_workers(config: pytest.Config) -> int:
    env_var = os.environ.get("PYTEST_XDIST_AUTO_NUM_WORKERS")
    if env_var:
        try:
            return int(env_var)
        except ValueError:
            warnings.warn(
                "PYTEST_XDIST_AUTO_NUM_WORKERS is not a number: {env_var!r}. Ignoring it."
            )

    try:
        import psutil
    except ImportError:
        pass
    else:
        use_logical: bool = config.option.numprocesses == "logical"
        count = psutil.cpu_count(logical=use_logical) or psutil.cpu_count()
        if count:
            return count
    try:
        from os import sched_getaffinity

        def cpu_count() -> int:
            return len(sched_getaffinity(0))

    except ImportError:
        if os.environ.get("TRAVIS") == "true":
            # workaround https://bitbucket.org/pypy/pypy/issues/2375
            return 2
        try:
            from os import cpu_count  # type: ignore[assignment]
        except ImportError:
            from multiprocessing import cpu_count
    try:
        n = cpu_count()
    except NotImplementedError:
        return 1
    return n if n else 1


def parse_numprocesses(s: str) -> int | Literal["auto", "logical"]:
    if s in ("auto", "logical"):
        return s  # type: ignore[return-value]
    elif s is not None:
        return int(s)


@pytest.hookimpl
def pytest_addoption(parser: pytest.Parser) -> None:
    # 'Help' formatting (same rules as pytest's):
    # Start with capitalized letters.
    # If a single phrase, do not end with period. If more than one phrase, all phrases end with periods.
    # Use \n to separate logical lines.
    group = parser.getgroup("xdist", "distributed and subprocess testing")
    group._addoption(
        "-n",
        "--numprocesses",
        dest="numprocesses",
        metavar="numprocesses",
        action="store",
        type=parse_numprocesses,
        help="Shortcut for '--dist=load --tx=NUM*popen'.\n"
        "With 'logical', attempt to detect logical CPU count (requires psutil, falls back to 'auto').\n"
        "With 'auto', attempt to detect physical CPU count. If physical CPU count cannot be determined, "
        "falls back to 1.\n"
        "Forced to 0 (disabled) when used with --pdb.",
    )
    group.addoption(
        "--maxprocesses",
        dest="maxprocesses",
        metavar="maxprocesses",
        action="store",
        type=int,
        help="Limit the maximum number of workers to process the tests when using --numprocesses "
        "with 'auto' or 'logical'",
    )
    group.addoption(
        "--max-worker-restart",
        action="store",
        default=None,
        dest="maxworkerrestart",
        help="Maximum number of workers that can be restarted "
        "when crashed (set to zero to disable this feature)",
    )
    group.addoption(
        "--dist",
        metavar="distmode",
        action="store",
        choices=[
            "each",
            "load",
            "loadscope",
            "loadfile",
            "loadgroup",
            "worksteal",
            "no",
        ],
        dest="dist",
        default="no",
        help=(
            "Set mode for distributing tests to exec environments.\n\n"
            "each: Send each test to all available environments.\n\n"
            "load: Load balance by sending any pending test to any"
            " available environment.\n\n"
            "loadscope: Load balance by sending pending groups of tests in"
            " the same scope to any available environment.\n\n"
            "loadfile: Load balance by sending test grouped by file"
            " to any available environment.\n\n"
            "loadgroup: Like 'load', but sends tests marked with 'xdist_group' to the same worker.\n\n"
            "worksteal: Split the test suite between available environments,"
            " then re-balance when any worker runs out of tests.\n\n"
            "(default) no: Run tests inprocess, don't distribute."
        ),
    )
    group.addoption(
        "--loadscope-reorder",
        dest="loadscopereorder",
        action="store_true",
        default=True,
        help=(
            "Pytest-xdist will default reorder tests by number of tests per scope "
            "when used in conjunction with loadscope.\n"
            "This option will enable loadscope reorder which will improve the "
            "parallelism of the test suite.\n"
            "However, the partial order of tests might not be retained.\n"
        ),
    )
    group.addoption(
        "--no-loadscope-reorder",
        dest="loadscopereorder",
        action="store_false",
        help=(
            "Pytest-xdist will default reorder tests by number of tests per scope "
            "when used in conjunction with loadscope.\n"
            "This option will disable loadscope reorder, "
            "and the partial order of tests can be retained.\n"
            "This is useful when pytest-xdist is used together with "
            "other plugins that specify tests in a specific order."
        ),
    )
    group.addoption(
        "--tx",
        dest="tx",
        action="append",
        default=[],
        metavar="xspec",
        help=(
            "Add a test execution environment. Some examples:\n"
            "--tx popen//python=python2.5 --tx socket=192.168.1.102:8888\n"
            "--tx ssh=user@codespeak.net//chdir=testcache"
        ),
    )
    group.addoption(
        "--px",
        dest="px",
        action="append",
        default=[],
        metavar="xspec",
        help=(
            "Add a proxy gateway to pass to test execution environments using `via`. Example:\n"
            "--px id=my_proxy//socket=192.168.1.102:8888 --tx 5*popen//via=my_proxy"
        ),
    )
    group._addoption(
        "-d",
        action="store_true",
        dest="distload",
        default=False,
        help="Load-balance tests. Shortcut for '--dist=load'.",
    )
    group.addoption(
        "--rsyncdir",
        action="append",
        default=[],
        metavar="DIR",
        help="Add directory for rsyncing to remote tx nodes",
    )
    group.addoption(
        "--rsyncignore",
        action="append",
        default=[],
        metavar="GLOB",
        help="Add expression for ignores when rsyncing to remote tx nodes",
    )
    group.addoption(
        "--testrunuid",
        action="store",
        help=(
            "Provide an identifier shared amongst all workers as the value of "
            "the 'testrun_uid' fixture.\n"
            "If not provided, 'testrun_uid' is filled with a new unique string "
            "on every test run."
        ),
    )
    group.addoption(
        "--maxschedchunk",
        action="store",
        type=int,
        help=(
            "Maximum number of tests scheduled in one step for --dist=load.\n"
            "Setting it to 1 will force pytest to send tests to workers one by "
            "one - might be useful for a small number of slow tests.\n"
            "Larger numbers will allow the scheduler to submit consecutive "
            "chunks of tests to workers - allows reusing fixtures.\n"
            "Due to implementation reasons, at least 2 tests are scheduled per "
            "worker at the start. Only later tests can be scheduled one by one.\n"
            "Unlimited if not set."
        ),
    )

    parser.addini(
        "rsyncdirs",
        "list of (relative) paths to be rsynced for remote distributed testing.",
        type="paths",
    )
    parser.addini(
        "rsyncignore",
        "list of (relative) glob-style paths to be ignored for rsyncing.",
        type="paths",
    )
    parser.addini(
        "looponfailroots",
        type="paths",
        help="directories to check for changes. Default: current directory.",
    )


# -------------------------------------------------------------------------
# distributed testing hooks
# -------------------------------------------------------------------------


@pytest.hookimpl
def pytest_addhooks(pluginmanager: pytest.PytestPluginManager) -> None:
    from xdist import newhooks

    pluginmanager.add_hookspecs(newhooks)


# -------------------------------------------------------------------------
# distributed testing initialization
# -------------------------------------------------------------------------


@pytest.hookimpl(trylast=True)
def pytest_configure(config: pytest.Config) -> None:
    config_line = (
        "xdist_group: specify group for tests should run in same session."
        "in relation to one another. Provided by pytest-xdist."
    )
    config.addinivalue_line("markers", config_line)

    # Skip this plugin entirely when only doing collection.
    if config.getvalue("collectonly"):
        return

    # Create the distributed session in case we have a valid distribution
    # mode and test environments.
    if _is_distribution_mode(config):
        from xdist.dsession import DSession

        session = DSession(config)
        config.pluginmanager.register(session, "dsession")
        tr = config.pluginmanager.getplugin("terminalreporter")
        if tr:
            tr.showfspath = False

    # Deprecation warnings for deprecated command-line/configuration options.
    if config.getoption("looponfail", None) or config.getini("looponfailroots"):
        warning = DeprecationWarning(
            "The --looponfail command line argument and looponfailroots config variable are deprecated.\n"
            "The loop-on-fail feature will be removed in pytest-xdist 4.0."
        )
        config.issue_config_time_warning(warning, 2)

    if config.getoption("rsyncdir", None) or config.getini("rsyncdirs"):
        warning = DeprecationWarning(
            "The --rsyncdir command line argument and rsyncdirs config variable are deprecated.\n"
            "The rsync feature will be removed in pytest-xdist 4.0."
        )
        config.issue_config_time_warning(warning, 2)


def _is_distribution_mode(config: pytest.Config) -> bool:
    """Whether distribution mode is on."""
    return config.getoption("dist") != "no" and bool(config.getoption("tx"))


@pytest.hookimpl(tryfirst=True)
def pytest_cmdline_main(config: pytest.Config) -> None:
    if config.option.distload:
        config.option.dist = "load"

    usepdb = config.getoption("usepdb", False)  # a core option
    if config.option.numprocesses in ("auto", "logical"):
        if usepdb:
            config.option.numprocesses = 0
            config.option.dist = "no"
        else:
            auto_num_cpus = config.hook.pytest_xdist_auto_num_workers(config=config)
            config.option.numprocesses = auto_num_cpus

    if config.option.numprocesses:
        if config.option.dist == "no":
            config.option.dist = "load"
        numprocesses = config.option.numprocesses
        if config.option.maxprocesses:
            numprocesses = min(numprocesses, config.option.maxprocesses)
        config.option.tx = ["popen"] * numprocesses

    if config.option.numprocesses == 0:
        config.option.dist = "no"
        config.option.tx = []

    val = config.getvalue
    if not val("collectonly") and _is_distribution_mode(config) and usepdb:
        raise pytest.UsageError(
            "--pdb is incompatible with distributing tests; try using -n0 or -nauto."
        )


# -------------------------------------------------------------------------
# fixtures and API to easily know the role of current node
# -------------------------------------------------------------------------


def is_xdist_worker(
    request_or_session: pytest.FixtureRequest | pytest.Session,
) -> bool:
    """Return `True` if this is an xdist worker, `False` otherwise.

    :param request_or_session: the `pytest` `request` or `session` object
    """
    return hasattr(request_or_session.config, "workerinput")


def is_xdist_controller(
    request_or_session: pytest.FixtureRequest | pytest.Session,
) -> bool:
    """Return `True` if this is the xdist controller, `False` otherwise.

    Note: this method also returns `False` when distribution has not been
    activated at all.

    :param request_or_session: the `pytest` `request` or `session` object
    """
    return (
        not is_xdist_worker(request_or_session)
        and request_or_session.config.option.dist != "no"
    )


# ALIAS: TODO, deprecate (#592)
is_xdist_master = is_xdist_controller


def get_xdist_worker_id(
    request_or_session: pytest.FixtureRequest | pytest.Session,
) -> str:
    """Return the id of the current worker ('gw0', 'gw1', etc) or 'master'
    if running on the controller node.

    If not distributing tests (for example passing `-n0` or not passing `-n` at all)
    also return 'master'.

    :param request_or_session: the `pytest` `request` or `session` object
    """
    if hasattr(request_or_session.config, "workerinput"):
        workerid: str = request_or_session.config.workerinput["workerid"]
        return workerid
    else:
        # TODO: remove "master", ideally for a None
        return "master"


@pytest.fixture(scope="session")
def worker_id(request: pytest.FixtureRequest) -> str:
    """Return the id of the current worker ('gw0', 'gw1', etc) or 'master'
    if running on the master node.
    """
    # TODO: remove "master", ideally for a None
    return get_xdist_worker_id(request)


@pytest.fixture(scope="session")
def testrun_uid(request: pytest.FixtureRequest) -> str:
    """Return the unique id of the current test."""
    if hasattr(request.config, "workerinput"):
        testrunid: str = request.config.workerinput["testrunuid"]
        return testrunid
    else:
        return uuid.uuid4().hex
