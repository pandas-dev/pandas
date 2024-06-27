# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import pytest

from asv.util import check_output, which, git_default_branch

from . import tools
from .conftest import generate_basic_conf

WIN = (os.name == 'nt')


@pytest.mark.skipif(tools.HAS_PYPY, reason="Times out randomly on pypy")
def test_find(capfd, tmpdir):
    values = [
        (None, None),
        (1, 1),
        (3, 1),
        (None, 1),
        (6, None),
        (5, 1),
        (6, 1),
        (6, 1),
        (6, 6),
        (6, 6),
    ]

    tmpdir, local, conf, machine_file = generate_basic_conf(tmpdir,
                                                            values=values,
                                                            dummy_packages=False)

    # Test find at least runs
    tools.run_asv_with_conf(conf, 'find',
                            f"{git_default_branch()}~5..{git_default_branch()}",
                            "params_examples.track_find_test",
                            _machine_file=machine_file)

    def get_result_files():
        results_dir = os.path.join(conf.results_dir, 'orangutan')
        if not os.path.isdir(results_dir):
            return ['machine.json']
        return sorted(os.listdir(results_dir))

    # Test find at least runs, not saving results
    prev_result_files = get_result_files()
    tools.run_asv_with_conf(conf, 'find', "--skip-save",
                            f"{git_default_branch()}~5..{git_default_branch()}",
                            "params_examples.track_find_test",
                            _machine_file=machine_file)
    assert get_result_files() == prev_result_files

    # Check it found the first commit after the initially tested one
    output, err = capfd.readouterr()

    regression_hash = check_output(
        [which('git'), 'rev-parse', f'{git_default_branch()}^'], cwd=conf.repo)

    assert f"Greatest regression found: {regression_hash[:8]}" in output

    # Test find at least runs, saving results --- should give same outcome
    tools.run_asv_with_conf(conf, 'find',
                            f"{git_default_branch()}~5..{git_default_branch()}",
                            "params_examples.track_find_test",
                            _machine_file=machine_file)
    output2, err2 = capfd.readouterr()
    assert get_result_files() == prev_result_files


@pytest.mark.skipif(tools.HAS_PYPY, reason="Times out randomly on pypy")
@pytest.mark.flaky(reruns=1, reruns_delay=5)  # depends on a timeout
def test_find_timeout(capfd, tmpdir):
    values = [
        (1, 0),
        (1, 0),
        (1, -1)
    ]

    tmpdir, local, conf, machine_file = generate_basic_conf(tmpdir,
                                                            values=values,
                                                            dummy_packages=False)

    # Test find at least runs
    tools.run_asv_with_conf(conf, 'find', "-e",
                            f"{git_default_branch()}",
                            "params_examples.time_find_test_timeout",
                            _machine_file=machine_file)

    # Check it found the first commit after the initially tested one
    output, err = capfd.readouterr()

    regression_hash = check_output(
        [which('git'), 'rev-parse', f'{git_default_branch()}'], cwd=conf.repo)

    assert f"Greatest regression found: {regression_hash[:8]}" in output
    assert "asv: benchmark timed out (timeout 1.0s)" in output


@pytest.mark.skipif(tools.HAS_PYPY, reason="Times out randomly on pypy")
def test_find_inverted(capfd, tmpdir):
    values = [
        (5, 6),
        (6, 6),
        (6, 6),
        (6, 1),
        (6, 1),
    ]

    tmpdir, local, conf, machine_file = generate_basic_conf(tmpdir,
                                                            values=values,
                                                            dummy_packages=False)
    tools.run_asv_with_conf(*[conf, 'find',
                              "-i", f"{git_default_branch()}~4..{git_default_branch()}",
                              "params_examples.track_find_test"],
                            _machine_file=machine_file)

    output, err = capfd.readouterr()

    regression_hash = check_output(
        [which('git'), 'rev-parse', f'{git_default_branch()}^'], cwd=conf.repo)

    formatted = f"Greatest improvement found: {regression_hash[:8]}"
    assert formatted in output
