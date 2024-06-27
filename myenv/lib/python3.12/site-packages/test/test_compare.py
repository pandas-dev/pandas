# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import re
import shutil
from os.path import abspath, dirname, join

import pytest

from asv import config, util
from asv.commands.compare import Compare

from . import tools

try:
    import hglib
except ImportError:
    hglib = None

from .tools import WIN

MACHINE_FILE = abspath(join(dirname(__file__), 'asv-machine.json'))

REFERENCE = """
All benchmarks:

| Change   | Before [22b920c6]    | After [fcf8c079]    | Ratio   | Benchmark (Parameter)                         |
|----------|----------------------|---------------------|---------|-----------------------------------------------|
| !        | n/a                  | failed              | n/a     | params_examples.ParamSuite.track_value        |
|          | failed               | failed              | n/a     | time_AAA_failure                              |
|          | n/a                  | n/a                 | n/a     | time_AAA_skip                                 |
|          | 1.00±1s              | 3.00±1s             | ~3.00   | time_ci_big                                   |
| +        | 1.00±0s              | 3.00±0s             | 3.00    | time_ci_small                                 |
| !        | 454μs                | failed              | n/a     | time_coordinates.time_latitude                |
|          | 1.00s                | 1.00s               | 1.00    | time_other.time_parameterized(1)              |
|          | 2.00s                | 4.00s               | 2.00    | time_other.time_parameterized(2)              |
| !        | 3.00s                | failed              | n/a     | time_other.time_parameterized(3)              |
| +        | 1.75ms               | 153ms               | 87.28   | time_quantity.time_quantity_array_conversion  |
| +        | 934μs                | 108ms               | 115.90  | time_quantity.time_quantity_init_array        |
|          | 83.6μs               | 55.4μs              | 0.66    | time_quantity.time_quantity_init_scalar       |
|          | 282μs                | 147μs               | 0.52    | time_quantity.time_quantity_scalar_conversion |
| +        | 1.31ms               | 7.75ms              | 5.91    | time_quantity.time_quantity_ufunc_sin         |
|          | 42                   | 42                  | 1.00    | time_secondary.track_value                    |
|          | 5.73m                | 5.73m               | 1.00    | time_units.mem_unit                           |
| +        | 125μs                | 3.81ms              | 30.42   | time_units.time_simple_unit_parse             |
|          | 1.64ms               | 1.53ms              | 0.93    | time_units.time_unit_compose                  |
| +        | 372μs                | 11.5ms              | 30.81   | time_units.time_unit_parse                    |
| -        | 69.1μs               | 18.3μs              | 0.27    | time_units.time_unit_to                       |
|          | 11.9μs               | 13.1μs              | 1.10    | time_units.time_very_simple_unit_parse        |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_match                       |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_mismatch_bench              |
| x        | 1.00s                | 3.00s               | 3.00    | time_with_version_mismatch_other              |
"""  # noqa: E501 line too long

REFERENCE_SPLIT = """
Benchmarks that have improved:

| Change   | Before [22b920c6]    | After [fcf8c079]    |   Ratio | Benchmark (Parameter)   |
|----------|----------------------|---------------------|---------|-------------------------|
| -        | 69.1μs               | 18.3μs              |    0.27 | time_units.time_unit_to |

Benchmarks that have stayed the same:

| Change   | Before [22b920c6]    | After [fcf8c079]    | Ratio   | Benchmark (Parameter)                         |
|----------|----------------------|---------------------|---------|-----------------------------------------------|
|          | failed               | failed              | n/a     | time_AAA_failure                              |
|          | n/a                  | n/a                 | n/a     | time_AAA_skip                                 |
|          | 1.00±1s              | 3.00±1s             | ~3.00   | time_ci_big                                   |
|          | 1.00s                | 1.00s               | 1.00    | time_other.time_parameterized(1)              |
|          | 2.00s                | 4.00s               | 2.00    | time_other.time_parameterized(2)              |
|          | 83.6μs               | 55.4μs              | 0.66    | time_quantity.time_quantity_init_scalar       |
|          | 282μs                | 147μs               | 0.52    | time_quantity.time_quantity_scalar_conversion |
|          | 42                   | 42                  | 1.00    | time_secondary.track_value                    |
|          | 5.73m                | 5.73m               | 1.00    | time_units.mem_unit                           |
|          | 1.64ms               | 1.53ms              | 0.93    | time_units.time_unit_compose                  |
|          | 11.9μs               | 13.1μs              | 1.10    | time_units.time_very_simple_unit_parse        |

Benchmarks that have got worse:

| Change   | Before [22b920c6]    | After [fcf8c079]    | Ratio   | Benchmark (Parameter)                        |
|----------|----------------------|---------------------|---------|----------------------------------------------|
| !        | n/a                  | failed              | n/a     | params_examples.ParamSuite.track_value       |
| +        | 1.00±0s              | 3.00±0s             | 3.00    | time_ci_small                                |
| !        | 454μs                | failed              | n/a     | time_coordinates.time_latitude               |
| !        | 3.00s                | failed              | n/a     | time_other.time_parameterized(3)             |
| +        | 1.75ms               | 153ms               | 87.28   | time_quantity.time_quantity_array_conversion |
| +        | 934μs                | 108ms               | 115.90  | time_quantity.time_quantity_init_array       |
| +        | 1.31ms               | 7.75ms              | 5.91    | time_quantity.time_quantity_ufunc_sin        |
| +        | 125μs                | 3.81ms              | 30.42   | time_units.time_simple_unit_parse            |
| +        | 372μs                | 11.5ms              | 30.81   | time_units.time_unit_parse                   |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_match                      |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_mismatch_bench             |

Benchmarks that are not comparable:

| Change   | Before [22b920c6]    | After [fcf8c079]    |   Ratio | Benchmark (Parameter)            |
|----------|----------------------|---------------------|---------|----------------------------------|
| x        | 1.00s                | 3.00s               |       3 | time_with_version_mismatch_other |
"""  # noqa: E501 line too long

REFERENCE_ONLY_CHANGED = """
| Change   | Before [22b920c6] <name1>   | After [fcf8c079] <name2>   | Ratio   | Benchmark (Parameter)                        |
|----------|-----------------------------|----------------------------|---------|----------------------------------------------|
| !        | n/a                         | failed                     | n/a     | params_examples.ParamSuite.track_value       |
| !        | 454μs                       | failed                     | n/a     | time_coordinates.time_latitude               |
| !        | 3.00s                       | failed                     | n/a     | time_other.time_parameterized(3)             |
| +        | 1.75ms                      | 153ms                      | 87.28   | time_quantity.time_quantity_array_conversion |
| +        | 1.31ms                      | 7.75ms                     | 5.91    | time_quantity.time_quantity_ufunc_sin        |
| +        | 372μs                       | 11.5ms                     | 30.81   | time_units.time_unit_parse                   |
| +        | 125μs                       | 3.81ms                     | 30.42   | time_units.time_simple_unit_parse            |
| +        | 1.00±0s                     | 3.00±0s                    | 3.00    | time_ci_small                                |
| +        | 1.00s                       | 3.00s                      | 3.00    | time_with_version_match                      |
| +        | 1.00s                       | 3.00s                      | 3.00    | time_with_version_mismatch_bench             |
| +        | 934μs                       | 108ms                      | 115.90  | time_quantity.time_quantity_init_array       |
| -        | 69.1μs                      | 18.3μs                     | 0.27    | time_units.time_unit_to                      |
"""  # noqa: E501 line too long

REFERENCE_ONLY_CHANGED_MULTIENV = """
| Change   | Before [22b920c6]    | After [fcf8c079]    | Ratio   | Benchmark (Parameter)                                                 |
|----------|----------------------|---------------------|---------|-----------------------------------------------------------------------|
| !        | n/a                  | failed              | n/a     | params_examples.ParamSuite.track_value [cheetah/py2.7-numpy1.8]       |
| !        | 454μs                | failed              | n/a     | time_coordinates.time_latitude [cheetah/py2.7-numpy1.8]               |
| !        | 3.00s                | failed              | n/a     | time_other.time_parameterized(3) [cheetah/py2.7-numpy1.8]             |
| +        | 1.75ms               | 153ms               | 87.28   | time_quantity.time_quantity_array_conversion [cheetah/py2.7-numpy1.8] |
| +        | 1.31ms               | 7.75ms              | 5.91    | time_quantity.time_quantity_ufunc_sin [cheetah/py2.7-numpy1.8]        |
| +        | 372μs                | 11.5ms              | 30.81   | time_units.time_unit_parse [cheetah/py2.7-numpy1.8]                   |
| +        | 125μs                | 3.81ms              | 30.42   | time_units.time_simple_unit_parse [cheetah/py2.7-numpy1.8]            |
| +        | 1.00±0s              | 3.00±0s             | 3.00    | time_ci_small [cheetah/py2.7-numpy1.8]                                |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_match [cheetah/py2.7-numpy1.8]                      |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_mismatch_bench [cheetah/py2.7-numpy1.8]             |
| +        | 934μs                | 108ms               | 115.90  | time_quantity.time_quantity_init_array [cheetah/py2.7-numpy1.8]       |
| -        | 69.1μs               | 18.3μs              | 0.27    | time_units.time_unit_to [cheetah/py2.7-numpy1.8]                      |
"""  # noqa: E501 line too long

REFERENCE_ONLY_CHANGED_NOSTATS = """
| Change   | Before [22b920c6]    | After [fcf8c079]    | Ratio   | Benchmark (Parameter)                        |
|----------|----------------------|---------------------|---------|----------------------------------------------|
| !        | n/a                  | failed              | n/a     | params_examples.ParamSuite.track_value       |
| !        | 454μs                | failed              | n/a     | time_coordinates.time_latitude               |
| !        | 3.00s                | failed              | n/a     | time_other.time_parameterized(3)             |
| +        | 1.75ms               | 153ms               | 87.28   | time_quantity.time_quantity_array_conversion |
| +        | 1.31ms               | 7.75ms              | 5.91    | time_quantity.time_quantity_ufunc_sin        |
| +        | 372μs                | 11.5ms              | 30.81   | time_units.time_unit_parse                   |
| +        | 125μs                | 3.81ms              | 30.42   | time_units.time_simple_unit_parse            |
| +        | 1.00±1s              | 3.00±1s             | 3.00    | time_ci_big                                  |
| +        | 1.00±0s              | 3.00±0s             | 3.00    | time_ci_small                                |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_match                      |
| +        | 1.00s                | 3.00s               | 3.00    | time_with_version_mismatch_bench             |
| +        | 934μs                | 108ms               | 115.90  | time_quantity.time_quantity_init_array       |
| -        | 69.1μs               | 18.3μs              | 0.27    | time_units.time_unit_to                      |
"""  # noqa: E501 line too long


def test_compare(capsys, tmpdir, example_results):
    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    conf = config.Config.from_json(
        {'results_dir': example_results,
         'repo': tools.generate_test_repo(tmpdir).path,
         'project': 'asv',
         'environment_type': "shouldn't matter what"})

    tools.run_asv_with_conf(conf, 'compare', '22b920c6', 'fcf8c079', '--machine=cheetah',
                            '--factor=2', '--environment=py2.7-numpy1.8')

    text, err = capsys.readouterr()
    assert text.strip() == REFERENCE.strip()

    tools.run_asv_with_conf(conf, 'compare', '22b920c6', 'fcf8c079', '--machine=cheetah',
                            '--factor=2', '--split', '--environment=py2.7-numpy1.8')
    text, err = capsys.readouterr()
    assert text.strip() == REFERENCE_SPLIT.strip()

    # Check print_table output as called from Continuous
    status = Compare.print_table(conf, '22b920c6', 'fcf8c079', factor=2, machine='cheetah',
                                 split=False, only_changed=True, sort='ratio',
                                 env_names=["py2.7-numpy1.8"],
                                 commit_names={'22b920c6': 'name1', 'fcf8c079': 'name2'})
    worsened, improved = status
    assert worsened
    assert improved
    text, err = capsys.readouterr()
    # Removing tailing spaces, so they don't need to be in `REFERENCE_ONLY_CHANGED`
    text = re.sub(r' *\n', r'\n', text)
    assert text.strip() == REFERENCE_ONLY_CHANGED.strip()

    # Check table with multiple environments
    status = Compare.print_table(conf, '22b920c6', 'fcf8c079', factor=2, machine='cheetah',
                                 split=False, only_changed=True, sort='ratio')
    text, err = capsys.readouterr()
    assert text.strip() == REFERENCE_ONLY_CHANGED_MULTIENV.strip()

    # Check results with no stats
    tools.run_asv_with_conf(conf, 'compare', '22b920c6', 'fcf8c079', '--machine=cheetah',
                            '--factor=2', '--sort=ratio', '--environment=py2.7-numpy1.8',
                            '--no-stats', '--only-changed')
    text, err = capsys.readouterr()
    assert text.strip() == REFERENCE_ONLY_CHANGED_NOSTATS.strip()
    assert "time_ci_big" in text.strip()


@pytest.mark.parametrize("dvcs_type", [
    "git",
    pytest.param("hg", marks=pytest.mark.skipif(hglib is None or WIN,
                                                reason="needs hglib"))
])
def test_compare_name_lookup(dvcs_type, capsys, tmpdir, example_results):
    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    repo = tools.generate_test_repo(tmpdir, dvcs_type=dvcs_type)

    branch_name = util.git_default_branch() if dvcs_type == 'git' else 'default'
    commit_hash = repo.get_branch_hashes(branch_name)[0]

    result_dir = os.path.join(tmpdir, 'results')

    src = os.path.join(example_results, 'cheetah')
    dst = os.path.join(result_dir, 'cheetah')
    os.makedirs(dst)

    for fn in ['feea15ca-py2.7-Cython-numpy1.8.json', 'machine.json']:
        shutil.copyfile(os.path.join(src, fn), os.path.join(dst, fn))

    shutil.copyfile(os.path.join(example_results, 'benchmarks.json'),
                    os.path.join(result_dir, 'benchmarks.json'))

    # Copy to different commit
    fn_1 = os.path.join(dst, 'feea15ca-py2.7-Cython-numpy1.8.json')
    fn_2 = os.path.join(dst, commit_hash[:8] + '-py2.7-Cython-numpy1.8.json')
    data = util.load_json(fn_1)
    data['commit_hash'] = commit_hash
    util.write_json(fn_2, data)

    conf = config.Config.from_json(
        {'results_dir': result_dir,
         'repo': repo.path,
         'project': 'asv',
         'environment_type': "shouldn't matter what"})

    # Lookup with symbolic name
    tools.run_asv_with_conf(conf, 'compare', branch_name, 'feea15ca', '--machine=cheetah',
                            '--factor=2', '--environment=py2.7-Cython-numpy1.8',
                            '--only-changed')

    # Nothing should be printed since no results were changed
    text, err = capsys.readouterr()
    assert text.strip() == ''
