import os
import threading
from pathlib import Path

import pytest
from unittest import mock

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.testing import subprocess_run_helper


_test_timeout = 60


def _test_cached_renderer():
    # Make sure that figures have an associated renderer after
    # a fig.canvas.draw() call
    fig = plt.figure(1)
    fig.canvas.draw()
    assert fig.canvas.get_renderer()._renderer is not None

    fig = plt.figure(2)
    fig.draw_without_rendering()
    assert fig.canvas.get_renderer()._renderer is not None


@pytest.mark.backend('macosx', skip_on_importerror=True)
def test_cached_renderer():
    subprocess_run_helper(_test_cached_renderer, timeout=_test_timeout,
                          extra_env={"MPLBACKEND": "macosx"})


def _test_savefig_rcparam():
    tmp_path = Path(os.environ["TEST_SAVEFIG_PATH"])

    def new_choose_save_file(title, directory, filename):
        # Replacement function instead of opening a GUI window
        # Make a new directory for testing the update of the rcParams
        assert directory == str(tmp_path)
        os.makedirs(f"{directory}/test")
        return f"{directory}/test/{filename}"

    fig = plt.figure()
    with (mock.patch("matplotlib.backends._macosx.choose_save_file",
                     new_choose_save_file),
          mpl.rc_context({"savefig.directory": tmp_path})):
        fig.canvas.toolbar.save_figure()
        # Check the saved location got created
        save_file = f"{tmp_path}/test/{fig.canvas.get_default_filename()}"
        assert os.path.exists(save_file)

        # Check the savefig.directory rcParam got updated because
        # we added a subdirectory "test"
        assert mpl.rcParams["savefig.directory"] == f"{tmp_path}/test"


@pytest.mark.backend('macosx', skip_on_importerror=True)
def test_savefig_rcparam(tmp_path):
    subprocess_run_helper(
        _test_savefig_rcparam, timeout=_test_timeout,
        extra_env={"MPLBACKEND": "macosx", "TEST_SAVEFIG_PATH": tmp_path})


@pytest.mark.backend('macosx', skip_on_importerror=True)
def test_ipython():
    from matplotlib.testing import ipython_in_subprocess
    ipython_in_subprocess("osx", {(8, 24): "macosx", (7, 0): "MacOSX"})


def _test_save_figure_return():
    fig, ax = plt.subplots()
    ax.imshow([[1]])
    prop = "matplotlib.backends._macosx.choose_save_file"
    with mock.patch(prop, return_value="foobar.png"):
        fname = fig.canvas.manager.toolbar.save_figure()
        os.remove("foobar.png")
        assert fname == "foobar.png"
    with mock.patch(prop, return_value=None):
        fname = fig.canvas.manager.toolbar.save_figure()
        assert fname is None


@pytest.mark.backend('macosx', skip_on_importerror=True)
def test_save_figure_return():
    subprocess_run_helper(_test_save_figure_return, timeout=_test_timeout,
                          extra_env={"MPLBACKEND": "macosx"})


def _test_create_figure_on_worker_thread_fails():
    def create_figure():
        warn_msg = "Matplotlib GUI outside of the main thread will likely fail."
        err_msg = "Cannot create a GUI FigureManager outside the main thread"
        with pytest.warns(UserWarning, match=warn_msg):
            with pytest.raises(RuntimeError, match=err_msg):
                plt.gcf()

    worker = threading.Thread(target=create_figure)
    worker.start()
    worker.join()


@pytest.mark.backend('macosx', skip_on_importerror=True)
def test_create_figure_on_worker_thread_fails():
    subprocess_run_helper(
        _test_create_figure_on_worker_thread_fails,
        timeout=_test_timeout,
        extra_env={"MPLBACKEND": "macosx"}
    )
