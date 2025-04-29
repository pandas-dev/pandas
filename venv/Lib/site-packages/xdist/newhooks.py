"""
xdist hooks.

Additionally, pytest-xdist will also decorate a few other hooks
with the worker instance that executed the hook originally:

``pytest_runtest_logreport``: ``rep`` parameter has a ``node`` attribute.

You can use this hooks just as you would use normal pytest hooks, but some care
must be taken in plugins in case ``xdist`` is not installed. Please see:

    http://pytest.org/en/latest/writing_plugins.html#optionally-using-hooks-from-3rd-party-plugins
"""

from __future__ import annotations

import os
from typing import Any
from typing import Sequence
from typing import TYPE_CHECKING

import execnet
import pytest


if TYPE_CHECKING:
    from xdist.remote import Producer
    from xdist.scheduler.protocol import Scheduling
    from xdist.workermanage import WorkerController


@pytest.hookspec()
def pytest_xdist_setupnodes(
    config: pytest.Config, specs: Sequence[execnet.XSpec]
) -> None:
    """Called before any remote node is set up."""


@pytest.hookspec()
def pytest_xdist_newgateway(gateway: execnet.Gateway) -> None:
    """Called on new raw gateway creation."""


@pytest.hookspec(
    warn_on_impl=DeprecationWarning(
        "rsync feature is deprecated and will be removed in pytest-xdist 4.0"
    )
)
def pytest_xdist_rsyncstart(
    source: str | os.PathLike[str],
    gateways: Sequence[execnet.Gateway],
) -> None:
    """Called before rsyncing a directory to remote gateways takes place."""


@pytest.hookspec(
    warn_on_impl=DeprecationWarning(
        "rsync feature is deprecated and will be removed in pytest-xdist 4.0"
    )
)
def pytest_xdist_rsyncfinish(
    source: str | os.PathLike[str],
    gateways: Sequence[execnet.Gateway],
) -> None:
    """Called after rsyncing a directory to remote gateways takes place."""


@pytest.hookspec(firstresult=True)
def pytest_xdist_getremotemodule() -> Any:
    """Called when creating remote node."""


@pytest.hookspec()
def pytest_configure_node(node: WorkerController) -> None:
    """Configure node information before it gets instantiated."""


@pytest.hookspec()
def pytest_testnodeready(node: WorkerController) -> None:
    """Test Node is ready to operate."""


@pytest.hookspec()
def pytest_testnodedown(node: WorkerController, error: object | None) -> None:
    """Test Node is down."""


@pytest.hookspec()
def pytest_xdist_node_collection_finished(
    node: WorkerController, ids: Sequence[str]
) -> None:
    """Called by the controller node when a worker node finishes collecting."""


@pytest.hookspec(firstresult=True)
def pytest_xdist_make_scheduler(
    config: pytest.Config, log: Producer
) -> Scheduling | None:
    """Return a node scheduler implementation."""


@pytest.hookspec(firstresult=True)
def pytest_xdist_auto_num_workers(config: pytest.Config) -> int:
    """
    Return the number of workers to spawn when ``--numprocesses=auto`` is given in the
    command-line.

    .. versionadded:: 2.1
    """
    raise NotImplementedError()


@pytest.hookspec(firstresult=True)
def pytest_handlecrashitem(
    crashitem: str, report: pytest.TestReport, sched: Scheduling
) -> None:
    """
    Handle a crashitem, modifying the report if necessary.

    The scheduler is provided as a parameter to reschedule the test if desired with
    `sched.mark_test_pending`.

    def pytest_handlecrashitem(crashitem, report, sched):
        if should_rerun(crashitem):
            sched.mark_test_pending(crashitem)
            report.outcome = "rerun"

    .. versionadded:: 2.2.1
    """
