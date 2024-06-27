# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import datetime
import shutil
from os.path import join

import pytest

from asv import results, runner, util


def _truncate_floats(item, digits=5):
    if isinstance(item, float):
        fmt = f'{{:.{digits - 1}e}}'
        return float(fmt.format(item))
    elif isinstance(item, list):
        return [_truncate_floats(x, digits) for x in item]
    elif isinstance(item, dict):
        return dict((k, _truncate_floats(v)) for k, v in item.items())
    else:
        return item


def test_results(tmpdir):
    tmpdir = str(tmpdir)

    timestamp1 = datetime.datetime.now(datetime.timezone.utc)
    duration = 1.5

    resultsdir = join(tmpdir, "results")
    for i in range(10):
        r = results.Results(
            {'machine': 'foo',
             'arch': 'x86_64'},
            {},
            hex(i),
            i * 1000000,
            '2.7',
            'some-environment-name',
            {})

        x1 = float(i * 0.001)
        x2 = float(i * 0.001)
        x3 = float((i + 1) ** -1)

        values = {
            'suite1.benchmark1': {'result': [x1], 'number': [1],
                                  'samples': [[x1, x1]], 'params': [['a']],
                                  'version': "1", 'profile': b'\x00\xff'},
            'suite1.benchmark2': {'result': [x2], 'number': [1],
                                  'samples': [[x2, x2, x2]], 'params': [],
                                  'version': "1", 'profile': b'\x00\xff'},
            'suite2.benchmark1': {'result': [x3], 'number': [None],
                                  'samples': [None], 'params': [['c']],
                                  'version': None, 'profile': b'\x00\xff'},
            'suite3.benchmark1': {'result': [x1, x2], 'number': [1, 1],
                                  'samples': [[x1, x1], [x2, x2, x3]],
                                  'params': [['c', 'd']],
                                  'version': None, 'profile': b'\x00\xff'}
        }

        for key, val in values.items():
            v = runner.BenchmarkResult(result=val['result'],
                                       samples=val['samples'],
                                       number=val['number'],
                                       profile=val['profile'],
                                       errcode=0,
                                       stderr='')
            benchmark = {'name': key, 'version': val['version'], 'params': val['params']}
            r.add_result(benchmark, v, record_samples=True,
                         started_at=timestamp1, duration=duration)

        # Save / add_existing_results roundtrip
        r.save(resultsdir)

        r2 = results.Results.load(join(resultsdir, r._filename))
        assert r2.date == r.date
        assert r2.commit_hash == r.commit_hash
        assert r2._filename == r._filename

        r3 = results.Results(r.params,
                             r._requirements,
                             r.commit_hash,
                             r.date,
                             r._python,
                             r.env_name,
                             {})
        r3.load_data(resultsdir)

        for rr in [r2, r3]:
            assert rr._results == r._results
            assert rr._stats == _truncate_floats(r._stats)
            assert rr._samples == r._samples
            assert rr._profiles == r._profiles
            assert rr.started_at == r._started_at
            assert rr.duration == _truncate_floats(r._duration)
            assert rr.benchmark_version == r._benchmark_version

        # Check the get_* methods
        assert sorted(r2.get_all_result_keys()) == sorted(values.keys())
        for bench in r2.get_all_result_keys():
            # Get with same parameters as stored
            params = r2.get_result_params(bench)
            assert params == values[bench]['params']
            assert r2.get_result_value(bench, params) == values[bench]['result']
            assert r2.get_result_samples(bench, params) == values[bench]['samples']
            stats = r2.get_result_stats(bench, params)
            if values[bench]['number'][0] is None:
                assert stats == [None]
            else:
                assert stats[0]['number'] == values[bench]['number'][0]

            # Get with different parameters than stored (should return n/a)
            bad_params = [['foo', 'bar']]
            assert r2.get_result_value(bench, bad_params) == [None, None]
            assert r2.get_result_stats(bench, bad_params) == [None, None]
            assert r2.get_result_samples(bench, bad_params) == [None, None]

            # Get profile
            assert r2.get_profile(bench) == b'\x00\xff'

        # Check get_result_keys
        mock_benchmarks = {
            'suite1.benchmark1': {'version': '1'},
            'suite1.benchmark2': {'version': '2'},
            'suite2.benchmark1': {'version': '2'},
        }
        assert sorted(r2.get_result_keys(mock_benchmarks)) == ['suite1.benchmark1',
                                                               'suite2.benchmark1']


def test_get_result_hash_from_prefix(tmpdir):
    results_dir = tmpdir.mkdir('results')
    machine_dir = results_dir.mkdir('cheetah')

    machine_json = join(os.path.dirname(__file__), 'example_results', 'cheetah', 'machine.json')
    shutil.copyfile(machine_json, join(str(machine_dir), 'machine.json'))

    for f in ['e5b6cdbc', 'e5bfoo12']:
        open(join(str(machine_dir), f'{f}-py2.7-Cython-numpy1.8.json'), 'a').close()

    # check unique, valid case
    full_commit = results.get_result_hash_from_prefix(str(results_dir), 'cheetah', 'e5b6')
    assert full_commit == 'e5b6cdbc'

    # check invalid hash case
    bad_commit = results.get_result_hash_from_prefix(str(results_dir), 'cheetah', 'foobar')
    assert bad_commit is None

    # check non-unique case
    with pytest.raises(util.UserError) as excinfo:
        results.get_result_hash_from_prefix(str(results_dir), 'cheetah', 'e')

    assert 'one of multiple commits' in str(excinfo.value)


def test_backward_compat_load(example_results): # noqa F811 redefinition of the imported fixture,it can be removed when the fixture is moved to conftest.py file and the import deleted
    resultsdir = example_results
    filename = join('cheetah', '624da0aa-py2.7-Cython-numpy1.8.json')

    r = results.Results.load(join(resultsdir, filename))
    assert r._filename == filename
    assert r._env_name == 'py2.7-Cython-numpy1.8'


def test_json_timestamp(tmpdir):
    # Check that per-benchmark timestamps are saved as JS timestamps in the result file
    tmpdir = str(tmpdir)

    stamp0 = datetime.datetime(
        1970, 1, 1,
        tzinfo = datetime.timezone.utc
    )
    stamp1 = datetime.datetime(
        1971, 1, 1,
        tzinfo = datetime.timezone.utc
    )
    duration = 1.5

    r = results.Results({'machine': 'mach'},
                        {},
                        'aaaa',
                        util.datetime_to_timestamp(stamp0),
                        'py',
                        'env',
                        {})
    value = runner.BenchmarkResult(
        result=[42],
        samples=[None],
        number=[None],
        profile=None,
        errcode=0,
        stderr=''
    )
    benchmark = {'name': 'some_benchmark', 'version': 'some version', 'params': []}
    r.add_result(benchmark, value, started_at=stamp1, duration=duration)
    r.save(tmpdir)

    r = util.load_json(join(tmpdir, 'mach', 'aaaa-env.json'))
    keys = r['result_columns']
    values = dict(zip(keys, r['results']['some_benchmark']))
    assert values['started_at'] == util.datetime_to_js_timestamp(stamp1)
    assert values['duration'] == duration


def test_iter_results(capsys, tmpdir, example_results): # noqa F811 noqa F811 redefinition of the imported fixture,it can be removed when the fixture is moved to conftest.py file and the import deleted
    dst = os.path.join(str(tmpdir), 'example_results')
    shutil.copytree(example_results, dst)

    path = os.path.join(dst, 'cheetah')

    skip_list = [
        'machine.json',
        'aaaaaaaa-py2.7-Cython-numpy1.8.json',  # malformed file
        'bbbbbbbb-py2.7-Cython-numpy1.8.json',  # malformed file
        'cccccccc-py2.7-Cython-numpy1.8.json',  # malformed file
    ]

    files = [f for f in os.listdir(path) if f.endswith('.json') and f not in skip_list]
    res = list(results.iter_results(path))
    assert len(res) == len(files)
    out, err = capsys.readouterr()
    assert skip_list[1] in out
    assert skip_list[2] in out
    assert skip_list[3] in out
    assert skip_list[0] not in out

    # The directory should be ignored without machine.json
    os.unlink(os.path.join(path, 'machine.json'))
    res = list(results.iter_results(path))
    assert len(res) == 0
    out, err = capsys.readouterr()
    assert "machine.json" in out


def test_filename_format():
    r = results.Results({'machine': 'foo'}, [], "commit", 0, "", "env", {})
    assert r._filename == join("foo", "commit-env.json")

    r = results.Results({'machine': 'foo'}, [], "hash", 0, "", "a" * 128, {})
    assert r._filename == join("foo", "hash-env-e510683b3f5ffe4093d021808bc6ff70.json")


def test_remove_samples():
    benchmark1 = {'name': 'a', 'version': '1', 'params': []}
    benchmark2 = {'name': 'b', 'version': '1', 'params': [['1', '2', '3']]}

    r = results.Results.unnamed()

    v1 = runner.BenchmarkResult(result=[True], samples=[[1]], number=[1],
                                profile=None, errcode=0, stderr='')
    v2 = runner.BenchmarkResult(result=[True] * 3, samples=[[1], [2], [3]], number=[1, 1, 1],
                                profile=None, errcode=0, stderr='')

    r.add_result(benchmark1, v1, record_samples=True)
    r.add_result(benchmark2, v2, record_samples=True)

    assert r.get_result_samples(benchmark1['name'], benchmark1['params']) == v1.samples
    assert r.get_result_samples(benchmark2['name'], benchmark2['params']) == v2.samples

    r.remove_samples(benchmark1['name'])
    assert r.get_result_samples(benchmark1['name'], benchmark1['params']) == [None]

    r.remove_samples(benchmark2['name'], selected_idx=[1])
    assert r.get_result_samples(benchmark2['name'], benchmark2['params']) == [[1], None, [3]]

    r.remove_samples(benchmark2['name'])
    assert r.get_result_samples(benchmark2['name'], benchmark2['params']) == [None, None, None]


def test_table_formatting():
    benchmark = {'params': [], 'param_names': [], 'unit': 's'}
    result = []
    expected = ["[]"]
    assert results._format_benchmark_result(result, benchmark) == expected

    benchmark = {'params': [['a', 'b', 'c']], 'param_names': ['param1'], "unit": "seconds"}
    result = list(zip([1e-6, 2e-6, 3e-6], [3e-6, 2e-6, 1e-6]))
    expected = ("======== ==========\n"
                " param1            \n"
                "-------- ----------\n"
                "   a      1.00\u00b13\u03bcs \n"
                "   b      2.00\u00b12\u03bcs \n"
                "   c      3.00\u00b11\u03bcs \n"
                "======== ==========")
    table = "\n".join(results._format_benchmark_result(result, benchmark, max_width=80))
    assert table == expected

    benchmark = {'params': [["'a'", "'b'", "'c'"], ["[1]", "[2]"]],
                 'param_names': ['param1', 'param2'], "unit": "seconds"}
    result = list(zip([1, 2, None, 4, 5, float('nan')], [None] * 6))
    expected = ("======== ======== =======\n"
                "--            param2     \n"
                "-------- ----------------\n"
                " param1    [1]      [2]  \n"
                "======== ======== =======\n"
                "   a      1.00s    2.00s \n"
                "   b      failed   4.00s \n"
                "   c      5.00s     n/a  \n"
                "======== ======== =======")
    table = "\n".join(results._format_benchmark_result(result, benchmark, max_width=80))
    assert table == expected

    expected = ("======== ======== ========\n"
                " param1   param2          \n"
                "-------- -------- --------\n"
                "   a       [1]     1.00s  \n"
                "   a       [2]     2.00s  \n"
                "   b       [1]     failed \n"
                "   b       [2]     4.00s  \n"
                "   c       [1]     5.00s  \n"
                "   c       [2]      n/a   \n"
                "======== ======== ========")
    table = "\n".join(results._format_benchmark_result(result, benchmark, max_width=0))
    assert table == expected
