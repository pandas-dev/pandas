# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import os
import re
import shutil
import glob
import datetime
import textwrap
from os.path import join

import pytest

from asv import results, environment, repo, util
from asv.commands.run import Run
from asv.commands import make_argparser

from . import tools
from .tools import WIN


def test_set_commit_hash(capsys, existing_env_conf):
    tmpdir, local, conf, machine_file = existing_env_conf

    r = repo.get_repo(conf)
    commit_hash = r.get_hash_from_name(r.get_branch_name())

    tools.run_asv_with_conf(conf, 'run', '--set-commit-hash=' + r.get_branch_name(),
                            _machine_file=join(tmpdir, 'asv-machine.json'))

    env_name = list(environment.get_environments(conf, None))[0].name
    result_filename = commit_hash[:conf.hash_length] + '-' + env_name + '.json'
    assert result_filename in os.listdir(join('results_workflow', 'orangutan'))

    result_path = join('results_workflow', 'orangutan', result_filename)
    times = results.Results.load(result_path)
    assert times.commit_hash == commit_hash


def test_run_spec(basic_conf_2):
    tmpdir, local, conf, machine_file = basic_conf_2
    conf.build_cache_size = 5

    extra_branches = [(f'{util.git_default_branch()}~1', 'some-branch', [12])]
    dvcs_path = os.path.join(tmpdir, 'test_repo2')
    dvcs = tools.generate_test_repo(dvcs_path, [1, 2],
                                    extra_branches=extra_branches)
    conf.repo = dvcs.path

    initial_commit = dvcs.get_hash(f"{util.git_default_branch()}~1")
    main_commit = dvcs.get_hash(f"{util.git_default_branch()}")
    branch_commit = dvcs.get_hash("some-branch")
    template_dir = os.path.join(tmpdir, "results_workflow_template")
    results_dir = os.path.join(tmpdir, 'results_workflow')
    tools.run_asv_with_conf(conf, 'run', initial_commit + "^!",
                            '--bench=time_secondary.track_value',
                            '--quick',
                            _machine_file=join(tmpdir, 'asv-machine.json'))
    shutil.copytree(results_dir, template_dir)

    def _test_run(range_spec, branches, expected_commits):
        # Rollback initial results
        shutil.rmtree(results_dir)
        shutil.copytree(template_dir, results_dir)

        args = ["run", "--quick", "--skip-existing-successful",
                "--bench=time_secondary.track_value",
                "-s", "1000"]  # large number of steps should be noop
        if range_spec is not None:
            args.append(range_spec)
        conf.branches = branches
        tools.run_asv_with_conf(conf, *args, _machine_file=machine_file)

        # Check that files for all commits expected were generated
        envs = list(environment.get_environments(conf, None))
        tool_name = envs[0].tool_name

        pyver = conf.pythons[0]
        if pyver.startswith('pypy'):
            pyver = pyver[2:]

        expected = set(['machine.json'])
        for commit in expected_commits:
            for psver in tools.DUMMY2_VERSIONS:
                expected.add(f'{commit[:8]}-{tool_name}-py{pyver}-asv_dummy_'
                             f'test_package_1-asv_dummy_test_package_2{psver}.json')

        result_files = os.listdir(join(tmpdir, 'results_workflow', 'orangutan'))

        assert set(result_files) == expected

    for branches, expected_commits in (
        # Without branches in config, shoud just use main
        ([None], [initial_commit, main_commit]),

        # With one branch in config, should just use that branch
        (["some-branch"], [initial_commit, branch_commit]),

        # With two branch in config, should apply to specified branches
        ([f"{util.git_default_branch()}", "some-branch"],
         [initial_commit, main_commit, branch_commit]),
    ):
        for range_spec in (None, "NEW", "ALL"):
            _test_run(range_spec, branches, expected_commits)

    # test the HASHFILE version of range_spec'ing
    expected_commits = (initial_commit, branch_commit)
    with open(os.path.join(tmpdir, 'hashes_to_benchmark'), 'w') as f:
        for commit in expected_commits:
            f.write(commit + "\n")
        f.write(f"{util.git_default_branch()}~1\n")
        f.write("some-bad-hash-that-will-be-ignored\n")
        expected_commits += (dvcs.get_hash(f"{util.git_default_branch()}~1"),)
    _test_run('HASHFILE:hashes_to_benchmark', [None], expected_commits)


def test_run_build_failure(basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    conf.matrix = {}

    # Add a commit that fails to build
    dvcs = tools.Git(conf.repo)
    setup_py = join(dvcs.path, 'setup.py')
    with open(setup_py, 'r') as f:
        setup_py_content = f.read()
    with open(setup_py, 'w') as f:
        f.write("assert False")
    dvcs.add(join(dvcs.path, 'setup.py'))
    dvcs.commit("Break setup.py")
    with open(setup_py, 'w') as f:
        f.write(setup_py_content)
    dvcs.add(join(dvcs.path, 'setup.py'))
    dvcs.commit("Fix setup.py")

    # Test running it
    timestamp = util.datetime_to_js_timestamp(datetime.datetime.now(datetime.timezone.utc))

    bench_name = 'time_secondary.track_value'
    for commit in [f'{util.git_default_branch()}^!',
                   f'{util.git_default_branch()}~1^!']:
        tools.run_asv_with_conf(conf, 'run', commit,
                                '--quick', '--show-stderr',
                                '--bench', bench_name,
                                _machine_file=machine_file)

    # Check results
    hashes = dvcs.get_branch_hashes()
    fn_broken, = glob.glob(join(tmpdir, 'results_workflow', 'orangutan',
                                hashes[1][:8] + '-*.json'))
    fn_ok, = glob.glob(join(tmpdir, 'results_workflow', 'orangutan',
                            hashes[0][:8] + '-*.json'))

    data_broken = util.load_json(fn_broken)
    data_ok = util.load_json(fn_ok)

    for data in (data_broken, data_ok):
        value = dict(zip(data['result_columns'], data['results'][bench_name]))
        assert value['started_at'] >= timestamp
        if data is data_broken:
            assert 'duration' not in value
        else:
            assert value['duration'] >= 0

    assert len(data_broken['results']) == 1
    assert len(data_ok['results']) == 1
    assert data_broken['result_columns'][0] == 'result'
    assert data_ok['result_columns'][0] == 'result'
    assert data_broken['results'][bench_name][0] is None
    assert data_ok['results'][bench_name][0] == [42.0]

    # Check that parameters were also saved
    assert data_broken['params'] == data_ok['params']


def test_run_with_repo_subdir(basic_conf_with_subdir):
    """
    Check 'asv run' with the Python project inside a subdirectory.
    """
    tmpdir, local, conf, machine_file = basic_conf_with_subdir

    conf.matrix = {}

    # This benchmark imports the project under test (asv_test_repo)
    bench_name = 'params_examples.track_find_test'
    # Test with a single changeset
    tools.run_asv_with_conf(conf, 'run', f'{util.git_default_branch()}^!',
                            '--quick', '--show-stderr',
                            '--bench', bench_name,
                            _machine_file=machine_file)

    # Check it ran ok
    fn_results, = glob.glob(join(tmpdir, 'results_workflow', 'orangutan',
                                 '*-*.json'))  # avoid machine.json
    data = util.load_json(fn_results)

    value = dict(zip(data['result_columns'], data['results'][bench_name]))
    assert value['params'] == [['1', '2']]
    assert value['result'] == [6, 6]


def test_benchmark_param_selection(basic_conf):
    tmpdir, local, conf, machine_file = basic_conf
    conf.matrix = {}
    tools.generate_test_repo(tmpdir, values=[(1, 2, 3)])
    tools.run_asv_with_conf(conf, 'run', f'{util.git_default_branch()}^!',
                            '--quick', '--show-stderr',
                            '--bench', r'track_param_selection\(.*, 3\)',
                            _machine_file=machine_file)

    def get_results():
        results = util.load_json(glob.glob(join(
            tmpdir, 'results_workflow', 'orangutan', '*-*.json'))[0])
        # replacing NaN by 'n/a' make assertions easier
        keys = results['result_columns']
        value = dict(zip(keys, results['results']['params_examples.track_param_selection']))
        return ['n/a' if util.is_nan(item) else item
                for item in value['result']]

    assert get_results() == [4, 'n/a', 5, 'n/a']
    tools.run_asv_with_conf(conf, 'run', '--show-stderr',
                            '--bench', r'track_param_selection\(1, ',
                            _machine_file=machine_file)
    assert get_results() == [4, 6, 5, 'n/a']
    tools.run_asv_with_conf(conf, 'run', '--show-stderr',
                            '--bench', 'track_param_selection',
                            _machine_file=machine_file)


def test_run_append_samples(basic_conf_2):
    tmpdir, local, conf, machine_file = basic_conf_2

    # Only one environment
    conf.matrix['asv_dummy_test_package_2'] = conf.matrix['asv_dummy_test_package_2'][:1]

    # Tests multiple calls to "asv run --append-samples"
    def run_it():
        tools.run_asv_with_conf(conf, 'run', f"{util.git_default_branch()}^!",
                                '--bench', 'time_examples.TimeSuite.time_example_benchmark_1',
                                '--append-samples', '-a', 'repeat=(1, 1, 10.0)', '-a', 'rounds=1',
                                '-a', 'number=1', '-a', 'warmup_time=0',
                                _machine_file=machine_file)

    run_it()

    result_dir = join(tmpdir, 'results_workflow', 'orangutan')
    result_fn, = [join(result_dir, fn) for fn in os.listdir(result_dir)
                  if fn != 'machine.json']

    data = util.load_json(result_fn)
    value = dict(zip(
        data['result_columns'],
        data['results']['time_examples.TimeSuite.time_example_benchmark_1']))
    assert value['stats_q_25'][0] is not None
    assert len(value['samples'][0]) == 1

    run_it()
    data = util.load_json(result_fn)
    value = dict(zip(
        data['result_columns'],
        data['results']['time_examples.TimeSuite.time_example_benchmark_1']))
    assert len(value['samples'][0]) == 2


def test_cpu_affinity(basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    # Only one environment
    conf.matrix = {}

    tools.run_asv_with_conf(conf, 'run', f"{util.git_default_branch()}^!",
                            '--bench', 'time_examples.TimeSuite.time_example_benchmark_1',
                            '--cpu-affinity=0', '-a', 'repeat=(1, 1, 10.0)', '-a', 'rounds=1',
                            '-a', 'number=1', '-a', 'warmup_time=0',
                            _machine_file=machine_file)
    # Check run produced a result
    result_dir = join(tmpdir, 'results_workflow', 'orangutan')
    result_fn, = [join(result_dir, fn) for fn in os.listdir(result_dir)
                  if fn != 'machine.json']
    data = util.load_json(result_fn)
    assert data['results']['time_examples.TimeSuite.time_example_benchmark_1']


def test_env_matrix_value(basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    conf.matrix = {}

    def check_env_matrix(env_build, env_nobuild):
        conf.matrix = {"env": env_build, "env_nobuild": env_nobuild}

        tools.run_asv_with_conf(conf, 'run', f"{util.git_default_branch()}^!",
                                '--bench', 'time_secondary.track_environment_value',
                                _machine_file=machine_file)

        # Check run produced a result
        result_dir = join(tmpdir, 'results_workflow', 'orangutan')

        result_fn1, = glob.glob(result_dir + '/*-SOME_TEST_VAR1.json')
        result_fn2, = glob.glob(result_dir + '/*-SOME_TEST_VAR2.json')

        data = util.load_json(result_fn1)
        assert data['result_columns'][0] == 'result'
        assert data['results']['time_secondary.track_environment_value'][0] == [1]

        data = util.load_json(result_fn2)
        assert data['results']['time_secondary.track_environment_value'][0] == [2]

    check_env_matrix({}, {'SOME_TEST_VAR': ['1', '2']})
    check_env_matrix({'SOME_TEST_VAR': ['1', '2']}, {})


@pytest.mark.skipif(tools.HAS_PYPY, reason="Times out randomly on pypy")
def test_parallel(basic_conf_2, dummy_packages):
    tmpdir, local, conf, machine_file = basic_conf_2

    if WIN and os.path.basename(sys.argv[0]).lower().startswith('py.test'):
        # Multiprocessing in spawn mode can result to problems with py.test
        # Find.run calls Setup.run in parallel mode by default
        pytest.skip("Multiprocessing spawn mode on Windows not safe to run "
                    "from py.test runner.")

    conf.matrix = {
        "req": dict(conf.matrix),
        "env": {"SOME_TEST_VAR": ["1", "2"]},
        "env_nobuild": {"SOME_OTHER_TEST_VAR": ["1", "2"]}
    }

    tools.run_asv_with_conf(conf, 'run', f"{util.git_default_branch()}^!",
                            '--bench', 'time_secondary.track_environment_value',
                            '--parallel=2', _machine_file=machine_file)


def test_filter_date_period(tmpdir, basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    dates = [
        datetime.datetime(2001, 1, 1, tzinfo=datetime.timezone.utc),
        datetime.datetime(2001, 1, 2, tzinfo=datetime.timezone.utc),
        datetime.datetime(2001, 1, 8, tzinfo=datetime.timezone.utc)
    ]

    dvcs = tools.generate_repo_from_ops(
        tmpdir, 'git',
        [("commit", j, dates[j]) for j in range(len(dates))])
    commits = dvcs.get_branch_hashes()[::-1]

    conf.repo = dvcs.path
    conf.matrix = {}

    tools.run_asv_with_conf(conf, 'run', f'{util.git_default_branch()}',
                            '--date-period=1w',
                            '--quick', '--show-stderr',
                            '--bench=time_secondary.track_value',
                            _machine_file=machine_file)

    expected_commits = [commits[0], commits[2]]

    fns = glob.glob(join(tmpdir, 'results_workflow', 'orangutan', '*-*.json'))

    for commit in expected_commits:
        assert any(os.path.basename(c).startswith(commit[:8]) for c in fns)

    assert len(fns) == len(expected_commits)


def test_format_durations():
    durations = {'foo': 1, 'bar': 2, 'quux': 3}

    msg = Run.format_durations(durations, 2)
    # Removing tailing spaces, so they don't need to be in `expected`
    msg = re.sub(r' *\n', r'\n', msg)
    expected = textwrap.dedent("""\
    =========== ================
     benchmark   total duration
    ----------- ----------------
        quux         3.00s
        bar          2.00s
        ...           ...
       total         6.00s
    =========== ================""")
    assert msg == expected


def test_return_code(tmpdir, basic_conf_2):
    tmpdir, local, conf, machine_file = basic_conf_2

    res = tools.run_asv_with_conf(conf, 'run', f'{util.git_default_branch()}^!', '--quick',
                                  '--bench', 'TimeSecondary',
                                  _machine_file=machine_file)
    assert res == 2


def test_run_python_same(capsys, basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    # Test Run runs with python=same
    tools.run_asv_with_conf(conf, 'run', '--python=same',
                            '--bench=time_secondary.TimeSecondary.time_exception',
                            '--bench=time_secondary.track_value',
                            _machine_file=machine_file)
    text, err = capsys.readouterr()

    assert re.search("time_exception.*failed", text, re.S)
    assert re.search(r"time_secondary.track_value\s+42.0", text)

    # Check that it did not clone or install
    assert "Cloning" not in text
    assert "Installing" not in text


def test_run_python_arg():
    parser, subparsers = make_argparser()

    argv = ['run', 'ALL']
    args = parser.parse_args(argv)
    assert args.env_spec == []


def test_run_steps_arg():
    parser, subparsers = make_argparser()

    argv = ['run', '--steps=20', 'ALL']
    args = parser.parse_args(argv)
    assert args.steps == 20
