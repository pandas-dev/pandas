# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import shutil
import textwrap
from os.path import join, dirname
from hashlib import sha256

import pytest

from asv import benchmarks, config, environment, util
from asv.repo import get_repo

from . import tools

BENCHMARK_DIR = join(dirname(__file__), 'benchmark')

INVALID_BENCHMARK_DIR = join(
    dirname(__file__), 'benchmark.invalid')

ASV_CONF_JSON = {
    'project': 'asv'
}

if util.ON_PYPY:
    ASV_CONF_JSON['pythons'] = ["pypy{0[0]}.{0[1]}".format(sys.version_info)]


def test_discover_benchmarks(benchmarks_fixture):
    conf, repo, envs, commit_hash = benchmarks_fixture

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex='secondary')
    assert len(b) == 6

    old_branches = conf.branches
    conf.branches = [f"{util.git_default_branch()}",
                     "some-missing-branch"]  # missing branches ignored
    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex='example')
    conf.branches = old_branches
    if util.ON_PYPY:
        assert len(b) == 34
    else:
        assert len(b) == 36

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex='time_example_benchmark_1')
    assert len(b) == 2

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex=['time_example_benchmark_1',
                                       'some regexp that does not match anything'])
    assert len(b) == 2

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash], regex='custom')
    assert sorted(b.keys()) == ['custom.time_function', 'custom.track_method',
                                'named.track_custom_pretty_name']
    assert 'pretty_name' not in b['custom.track_method']
    assert b['custom.time_function']['pretty_name'] == 'My Custom Function'
    assert b['named.track_custom_pretty_name']['pretty_name'] == 'this.is/the.answer'

    # benchmark param selection with regex
    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex=r'track_param_selection\(.*, 3\)')
    assert list(b.keys()) == ['params_examples.track_param_selection']
    assert b._benchmark_selection['params_examples.track_param_selection'] == [0, 2]
    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex=r'track_param_selection\(1, ')
    assert list(b.keys()) == ['params_examples.track_param_selection']
    assert b._benchmark_selection['params_examples.track_param_selection'] == [0, 1]
    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex='track_param_selection')
    assert list(b.keys()) == ['params_examples.track_param_selection']
    assert b._benchmark_selection['params_examples.track_param_selection'] == [0, 1, 2, 3]

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash])
    if not util.ON_PYPY:
        assert len(b) == 50
    else:
        assert len(b) == 48

    assert 'named.OtherSuite.track_some_func' in b

    params = b['params_examples.FunctionParamSuite.time_func']['params']
    assert len(params) == 1
    assert len(params[0]) == 2
    assert params[0][0] == '<function track_param>'
    # repr is a bit different on py2 vs py3 here
    assert params[0][1] in ['<function FunctionParamSuite.<lambda>>', '<function <lambda>>']

    # Raw timing benchmarks
    assert b['timeraw_examples.TimerawSuite.timeraw_count']['repeat'] == 7
    assert b['timeraw_examples.TimerawSuite.timeraw_count']['number'] == 3
    assert b['timeraw_examples.TimerawSuite.timeraw_setup']['number'] == 1


def test_invalid_benchmark_tree(tmpdir):
    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    d = {}
    d.update(ASV_CONF_JSON)
    d['benchmark_dir'] = INVALID_BENCHMARK_DIR
    d['env_dir'] = "env"
    d['repo'] = tools.generate_test_repo(tmpdir, [0]).path
    conf = config.Config.from_json(d)

    repo = get_repo(conf)
    envs = list(environment.get_environments(conf, None))
    commit_hash = repo.get_hash_from_name(repo.get_branch_name())

    with pytest.raises(util.UserError):
        benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash])


def test_find_benchmarks_cwd_imports(tmpdir):
    # Test that files in the directory above the benchmark suite are
    # not importable

    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    os.makedirs('benchmark')
    with open(os.path.join('benchmark', '__init__.py'), 'w') as f:
        pass

    with open(os.path.join('benchmark', 'test.py'), 'w') as f:
        f.write("""
try:
    import this_should_really_not_be_here
    raise AssertionError('This should not happen!')
except ImportError:
    pass

def track_this():
    return 0
""")

    with open(os.path.join('this_should_really_not_be_here.py'), 'w') as f:
        f.write("raise AssertionError('Should not be imported!')")

    d = {}
    d.update(ASV_CONF_JSON)
    d['env_dir'] = "env"
    d['benchmark_dir'] = 'benchmark'
    d['repo'] = tools.generate_test_repo(tmpdir, [[0, 1]]).path
    conf = config.Config.from_json(d)

    repo = get_repo(conf)
    envs = list(environment.get_environments(conf, None))
    commit_hash = repo.get_hash_from_name(repo.get_branch_name())

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex='track_this')
    assert len(b) == 1


def test_import_failure_retry(tmpdir):
    # Test that a different commit is tried on import failure

    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    os.makedirs('benchmark')
    with open(os.path.join('benchmark', '__init__.py'), 'w') as f:
        f.write(textwrap.dedent("""
        import asv_test_repo

        def time_foo():
            pass

        time_foo.number = asv_test_repo.dummy_value

        if asv_test_repo.dummy_value == 0:
            raise RuntimeError("fail discovery")
        """))

    dvcs = tools.generate_test_repo(tmpdir, [2, 1, 0])

    d = {}
    d.update(ASV_CONF_JSON)
    d['env_dir'] = "env"
    d['benchmark_dir'] = 'benchmark'
    d['repo'] = dvcs.path
    conf = config.Config.from_json(d)

    repo = get_repo(conf)
    envs = list(environment.get_environments(conf, None))
    commit_hashes = dvcs.get_branch_hashes()

    b = benchmarks.Benchmarks.discover(conf, repo, envs, commit_hashes)
    assert len(b) == 1
    assert b['time_foo']['number'] == 1


def test_conf_inside_benchmarks_dir(tmpdir):
    # Test that the configuration file can be inside the benchmark suite

    tmpdir = str(tmpdir)
    benchmark_dir = os.path.join(tmpdir, 'benchmark')

    os.makedirs(benchmark_dir)
    with open(os.path.join(benchmark_dir, '__init__.py'), 'w') as f:
        # Test also benchmark in top-level __init__.py
        f.write("def track_this(): pass")

    with open(os.path.join(benchmark_dir, 'bench.py'), 'w') as f:
        f.write("def track_this(): pass")

    d = {}
    d.update(ASV_CONF_JSON)
    d['env_dir'] = "env"
    d['benchmark_dir'] = '.'
    d['repo'] = tools.generate_test_repo(tmpdir, [[0, 1]]).path
    conf = config.Config.from_json(d)

    # NB. conf_dir == getcwd()
    os.chdir(benchmark_dir)

    repo = get_repo(conf)
    envs = list(environment.get_environments(conf, None))
    commit_hash = repo.get_hash_from_name(repo.get_branch_name())

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex='track_this')
    assert set(b.keys()) == {'track_this', 'bench.track_this'}


def test_code_extraction(tmpdir):
    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    shutil.copytree(BENCHMARK_DIR, 'benchmark')

    d = {}
    d.update(ASV_CONF_JSON)
    d['env_dir'] = "env"
    d['benchmark_dir'] = 'benchmark'
    d['repo'] = tools.generate_test_repo(tmpdir, [0]).path
    conf = config.Config.from_json(d)

    repo = get_repo(conf)
    envs = list(environment.get_environments(conf, None))
    commit_hash = repo.get_hash_from_name(repo.get_branch_name())

    b = benchmarks.Benchmarks.discover(conf, repo, envs, [commit_hash],
                                       regex=r'^code_extraction\.')

    expected_code = textwrap.dedent("""
    def track_test():
        # module-level 難
        return 0

    def setup():
        # module-level
        pass

    def setup_cache():
        # module-level
        pass
    """).strip()

    bench = b['code_extraction.track_test']
    assert bench['version'] == sha256(bench['code'].encode('utf-8')).hexdigest()
    assert bench['code'] == expected_code

    expected_code = textwrap.dedent("""
    int track_pretty_source_test() {
        return 0;
    }

    def setup():
        # module-level
        pass

    def setup_cache():
        # module-level
        pass
    """).strip()

    bench = b['code_extraction.track_pretty_source_test']
    assert bench['version'] == sha256(bench['code'].encode('utf-8')).hexdigest()
    assert bench['code'] == expected_code

    expected_code = textwrap.dedent("""
    class MyClass:
        def track_test(self):
            # class-level 難
            return 0

    def setup():
        # module-level
        pass

    class MyClass:
        def setup(self):
            # class-level
            pass

        def setup_cache(self):
            # class-level
            pass
    """).strip()

    bench = b['code_extraction.MyClass.track_test']
    assert bench['version'] == sha256(bench['code'].encode('utf-8')).hexdigest()

    if sys.version_info[:2] != (3, 2):
        # Python 3.2 doesn't have __qualname__
        assert bench['code'] == expected_code


def test_asv_benchmark_timings():
    # Check the benchmark runner runs
    util.check_call([sys.executable, '-masv.benchmark', 'timing',
                     '--setup=import time',
                     'time.sleep(0)'],
                    cwd=os.path.join(os.path.dirname(__file__), '..'))
