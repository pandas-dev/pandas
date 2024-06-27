# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import re

from asv.results import iter_results_for_machine
from asv import util

from . import tools
from .tools import get_default_environment_type


def test_continuous(capfd, basic_conf_2):
    tmpdir, local, conf, machine_file = basic_conf_2

    python = "{0[0]}.{0[1]}".format(sys.version_info)
    env_type = get_default_environment_type(conf, python)
    env_spec = ("-E", env_type + ":" + python)

    # Check that asv continuous runs
    tools.run_asv_with_conf(conf, 'continuous',
                            f"{util.git_default_branch()}^", '--show-stderr',
                            '--bench=params_examples.track_find_test',
                            '--bench=params_examples.track_param',
                            '--bench=time_examples.TimeSuite.time_example_benchmark_1',
                            '--attribute=repeat=1', '--attribute=number=1',
                            '--attribute=warmup_time=0',
                            *env_spec, _machine_file=machine_file)

    text, err = capfd.readouterr()
    assert "SOME BENCHMARKS HAVE CHANGED SIGNIFICANTLY" in text
    assert "PERFORMANCE INCREASED" in text or "PERFORMANCE DECREASED" in text
    # Check output, but not the whole row
    pattern = r"params_examples\.track_find_test\(2\) \[orangutan\/\w+-py"
    assert re.search(pattern, text) is not None

    assert "params_examples.ClassOne" in text

    # Check rounds were interleaved (timing benchmark was run twice)
    assert re.search(r"For.*commit [a-f0-9]+ (<[a-z0-9~^]+> )?\(round 1/2\)", text, re.M), text

    result_found = False
    for results in iter_results_for_machine(conf.results_dir, "orangutan"):
        result_found = True
        stats = results.get_result_stats('time_examples.TimeSuite.time_example_benchmark_1', [])
        assert stats[0]['repeat'] == 2
    assert result_found
