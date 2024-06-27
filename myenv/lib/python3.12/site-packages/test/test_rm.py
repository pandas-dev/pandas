# Licensed under a 3-clause BSD style license - see LICENSE.rst
import shutil
from os.path import join

from asv import config, results

from . import tools


def test_rm(tmpdir, example_results):
    tmpdir = str(tmpdir)

    shutil.copytree(
        example_results,
        join(tmpdir, 'example_results'))

    conf = config.Config.from_json({
        'results_dir': join(tmpdir, 'example_results'),
        'repo': "### IGNORED, BUT REQUIRED ###"
    })

    tools.run_asv_with_conf(conf, 'rm', '-y', 'benchmark=time_quantity*')

    results_a = list(results.iter_results(tmpdir))
    for result in results_a:
        for key in result.get_all_result_keys():
            assert not key.startswith('time_quantity')
        for key in result.started_at.keys():
            assert not key.startswith('time_quantity')
        for key in result.duration.keys():
            assert not key.startswith('time_quantity')

    tools.run_asv_with_conf(conf, 'rm', '-y', 'commit_hash=05d283b9')

    results_b = list(results.iter_results(tmpdir))
    assert len(results_b) == len(results_a) - 1
