# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import sys
import json

import pytest

from asv import config, environment, util
from asv.repo import get_repo
from asv.util import shlex_quote as quote

from . import tools
from .tools import (PYTHON_VER1, PYTHON_VER2, DUMMY1_VERSION, DUMMY2_VERSIONS, WIN, HAS_PYPY,
                    HAS_CONDA, HAS_VIRTUALENV, HAS_PYTHON_VER2, generate_test_repo)


@pytest.mark.skipif(not (HAS_PYTHON_VER2 or HAS_CONDA),
                    reason="Requires two usable Python versions")
def test_matrix_environments(tmpdir, dummy_packages):
    conf = config.Config()

    conf.env_dir = str(tmpdir.join("env"))

    conf.pythons = [PYTHON_VER1, PYTHON_VER2]
    conf.matrix = {
        "asv_dummy_test_package_1": [DUMMY1_VERSION, None],
        "asv_dummy_test_package_2": DUMMY2_VERSIONS
    }
    environments = list(environment.get_environments(conf, None))

    assert len(environments) == 2 * 2 * 2

    # Only test the first two environments, since this is so time
    # consuming
    for env in environments[:2]:
        env.create()

        output = env.run(
            ['-c', 'import asv_dummy_test_package_1 as p, sys; sys.stdout.write(p.__version__)'],
            valid_return_codes=None)
        if 'asv_dummy_test_package_1' in env._requirements:
            assert output.startswith(str(env._requirements['asv_dummy_test_package_1']))

        output = env.run(
            ['-c', 'import asv_dummy_test_package_2 as p, sys; sys.stdout.write(p.__version__)'])
        assert output.startswith(str(env._requirements['asv_dummy_test_package_2']))


@pytest.mark.skipif((not HAS_CONDA),
                    reason="Requires conda and conda-build")
def test_large_environment_matrix(tmpdir):
    # As seen in issue #169, conda can't handle using really long
    # directory names in its environment.  This creates an environment
    # with many dependencies in order to ensure it still works.

    conf = config.Config()

    conf.env_dir = str(tmpdir.join("env"))
    conf.pythons = [PYTHON_VER1]
    for i in range(25):
        conf.matrix[f'foo{i}'] = []

    environments = list(environment.get_environments(conf, None))

    for env in environments:
        # Since *actually* installing all the dependencies would make
        # this test run a long time, we only set up the environment,
        # but don't actually install dependencies into it.  This is
        # enough to trigger the bug in #169.
        env._get_requirements = lambda *a: ([], [])
        # pip / virtualenv setup still uses
        # _install_requirements
        env._install_requirements = lambda *a: None
        env.create()


@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda and conda-build")
def test_presence_checks(tmpdir, monkeypatch):
    conf = config.Config()

    if WIN:
        # Tell conda to not use hardlinks: on Windows it's not possible
        # to delete hard links to files in use, which causes problem when
        # trying to cleanup environments during this test
        monkeypatch.setenv(str('CONDA_ALWAYS_COPY'), str('True'))

    conf.env_dir = str(tmpdir.join("env"))

    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}
    environments = list(environment.get_environments(conf, None))

    for env in environments:
        env.create()
        assert env.check_presence()

        # Check env is recreated when info file is clobbered
        info_fn = os.path.join(env._path, 'asv-env-info.json')
        data = util.load_json(info_fn)
        data['python'] = '0'
        data = util.write_json(info_fn, data)
        env._is_setup = False
        env.create()
        data = util.load_json(info_fn)
        assert data['python'] == PYTHON_VER1
        env.run(['-c', 'import os'])

        # Check env is recreated if crucial things are missing
        pip_fns = [
            os.path.join(env._path, 'bin', 'pip')
        ]
        if WIN:
            pip_fns += [
                os.path.join(env._path, 'bin', 'pip.exe'),
                os.path.join(env._path, 'Scripts', 'pip'),
                os.path.join(env._path, 'Scripts', 'pip.exe')
            ]

        some_removed = False
        for pip_fn in pip_fns:
            if os.path.isfile(pip_fn):
                some_removed = True
                os.remove(pip_fn)
        assert some_removed

        env._is_setup = False
        env.create()
        assert os.path.isfile(pip_fn)
        env.run(['-c', 'import os'])


def _sorted_dict_list(lst):
    return list(sorted(lst, key=lambda x: list(sorted(x.items()))))


def test_matrix_expand_basic():
    conf = config.Config()
    conf.environment_type = 'something'
    conf.pythons = ["2.6", "2.7"]
    conf.matrix = {
        'pkg1': None,
        'pkg2': '',
        'pkg3': [''],
        'pkg4': ['1.2', '3.4'],
        'pkg5': []
    }

    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = _sorted_dict_list([
        {('python', None): '2.6', ('req', 'pkg2'): '', ('req', 'pkg3'): '',
         ('req', 'pkg4'): '1.2', ('req', 'pkg5'): ''},
        {('python', None): '2.6', ('req', 'pkg2'): '', ('req', 'pkg3'): '',
         ('req', 'pkg4'): '3.4', ('req', 'pkg5'): ''},
        {('python', None): '2.7', ('req', 'pkg2'): '', ('req', 'pkg3'): '',
         ('req', 'pkg4'): '1.2', ('req', 'pkg5'): ''},
        {('python', None): '2.7', ('req', 'pkg2'): '', ('req', 'pkg3'): '',
         ('req', 'pkg4'): '3.4', ('req', 'pkg5'): ''},
    ])
    assert combinations == expected


def test_matrix_expand_include():
    conf = config.Config()
    conf.environment_type = 'something'
    conf.pythons = ["2.6"]
    conf.matrix = {'a': '1'}
    conf.include = [
        {'python': '3.5', 'b': '2'},
        {'sys_platform': sys.platform, 'python': '2.7', 'b': '3'},
        {'sys_platform': sys.platform + 'nope', 'python': '2.7', 'b': '3'},
        {'environment_type': 'nope', 'python': '2.7', 'b': '4'},
        {'environment_type': 'something', 'python': '2.7', 'b': '5'},
    ]

    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = _sorted_dict_list([
        {('python', None): '2.6', ('req', 'a'): '1'},
        {('python', None): '3.5', ('req', 'b'): '2'},
        {('python', None): '2.7', ('req', 'b'): '3'},
        {('python', None): '2.7', ('req', 'b'): '5'}
    ])
    assert combinations == expected

    conf.include = [
        {'b': '2'}
    ]
    with pytest.raises(util.UserError):
        list(environment.iter_matrix(conf.environment_type, conf.pythons, conf))

@pytest.mark.skipif(not (HAS_PYTHON_VER2 or HAS_CONDA),
                    reason="Requires two usable Python versions")
def test_matrix_expand_include_detect_env_type():
    conf = config.Config()
    conf.environment_type = None
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}
    conf.exclude = [{}]
    conf.include = [
        {'sys_platform': sys.platform, 'python': PYTHON_VER1},
    ]

    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = _sorted_dict_list([
        {('python', None): PYTHON_VER1},
    ])
    assert combinations == expected


def test_matrix_expand_exclude():
    conf = config.Config()
    conf.environment_type = 'something'
    conf.pythons = ["2.6", "2.7"]
    conf.matrix = {
        'a': '1',
        'b': ['1', None]
    }
    conf.include = [
        {'python': '2.7', 'b': '2', 'c': None}
    ]

    # check basics
    conf.exclude = [
        {'python': '2.7', 'b': '2'},
        {'python': '2.7', 'b': None},
        {'python': '2.6', 'a': '1'},
    ]
    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = _sorted_dict_list([
        {('python', None): '2.7', ('req', 'a'): '1', ('req', 'b'): '1'},
        {('python', None): '2.7', ('req', 'b'): '2'}
    ])
    assert combinations == expected

    # check regexp
    conf.exclude = [
        {'python': '.*', 'b': None},
    ]
    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = _sorted_dict_list([
        {('python', None): '2.6', ('req', 'a'): '1', ('req', 'b'): '1'},
        {('python', None): '2.7', ('req', 'a'): '1', ('req', 'b'): '1'},
        {('python', None): '2.7', ('req', 'b'): '2'}
    ])
    assert combinations == expected

    # check environment_type as key
    conf.exclude = [
        {'environment_type': 'some.*'},
    ]
    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = [
        {('python', None): '2.7', ('req', 'b'): '2'}
    ]
    assert combinations == expected

    # check sys_platform as key
    conf.exclude = [
        {'sys_platform': sys.platform},
    ]
    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = [
        {('python', None): '2.7', ('req', 'b'): '2'}
    ]
    assert combinations == expected

    # check inverted regex
    conf.exclude = [
        {'python': '(?!2.6).*'}
    ]
    combinations = _sorted_dict_list(environment.iter_matrix(
        conf.environment_type, conf.pythons, conf))
    expected = _sorted_dict_list([
        {('python', None): '2.6', ('req', 'a'): '1', ('req', 'b'): '1'},
        {('python', None): '2.6', ('req', 'a'): '1'},
        {('python', None): '2.7', ('req', 'b'): '2'}
    ])
    assert combinations == expected


def test_iter_env_matrix_combinations():
    conf = config.Config()
    conf.environment_type = 'something'
    conf.pythons = ["2.6"]
    conf.matrix = {}
    conf.include = []

    # (matrix, expected)
    env_matrices = [
        ({'var0': ['val0', 'val1'], 'var1': ['val2', 'val3']},
         [{'var0': 'val0', 'var1': 'val2'},
          {'var0': 'val0', 'var1': 'val3'},
          {'var0': 'val1', 'var1': 'val2'},
          {'var0': 'val1', 'var1': 'val3'}]),
        ({'var0': ['val0', 'val1'], 'var1': ['val2', None]},
         [{'var0': 'val0', 'var1': 'val2'}, {'var0': 'val0'},
          {'var0': 'val1', 'var1': 'val2'}, {'var0': 'val1'}]),
        ({'var0': ['val0', 'val1']},
         [{'var0': 'val0'}, {'var0': 'val1'}]),
        ({}, [{}]),
    ]

    for matrix, expected in env_matrices:
        conf.matrix = {'env': matrix}
        expected = [{('env', key): value for key, value in item.items()}
                    for item in expected]
        for m in expected:
            m['python', None] = "2.6"
        result = _sorted_dict_list(environment.iter_matrix(conf.environment_type,
                                                           conf.pythons, conf))
        assert result == _sorted_dict_list(expected)


@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda and conda-build")
def test_conda_pip_install(tmpdir, dummy_packages):
    # test that we can install with pip into a conda environment.
    conf = config.Config()

    conf.env_dir = str(tmpdir.join("env"))

    conf.environment_type = "conda"
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {
        "pip+asv_dummy_test_package_2": DUMMY2_VERSIONS[0]
    }
    environments = list(environment.get_environments(conf, None))

    assert len(environments) == 1 * 1 * 1

    for env in environments:
        env.create()

        output = env.run(
            ['-c', 'import asv_dummy_test_package_2 as p, sys; sys.stdout.write(p.__version__)'])
        assert output.startswith(DUMMY2_VERSIONS[0])


@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda and conda-build")
def test_conda_environment_file(tmpdir, dummy_packages):
    env_file_name = str(tmpdir.join("environment.yml"))
    with open(env_file_name, "w") as temp_environment_file:
        temp_environment_file.write('name: test_conda_envs\ndependencies:'
                                    '\n  - asv_dummy_test_package_2')

    conf = config.Config()
    conf.env_dir = str(tmpdir.join("env"))
    conf.environment_type = "conda"
    conf.pythons = [PYTHON_VER1]
    conf.conda_environment_file = env_file_name
    conf.matrix = {
        "asv_dummy_test_package_1": [DUMMY1_VERSION]
    }

    environments = list(environment.get_environments(conf, None))

    assert len(environments) == 1 * 1 * 1

    for env in environments:
        env.create()

        output = env.run(
            ['-c', 'import asv_dummy_test_package_1 as p, sys; sys.stdout.write(p.__version__)'])
        assert output.startswith(str(DUMMY1_VERSION))

        output = env.run(
            ['-c', 'import asv_dummy_test_package_2 as p, sys; sys.stdout.write(p.__version__)'])
        assert output.startswith(str(DUMMY2_VERSIONS[1]))


@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda and conda-build")
def test_conda_run_executable(tmpdir):
    # test that we can install with pip into a conda environment.
    conf = config.Config()

    conf.env_dir = str(tmpdir.join("env"))

    conf.environment_type = "conda"
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}
    environments = list(environment.get_environments(conf, None))

    assert len(environments) == 1 * 1 * 1

    for env in environments:
        env.create()
        env.run_executable('conda', ['info'])


@pytest.mark.skipif(not (HAS_PYTHON_VER2 or HAS_CONDA),
                    reason="Requires two usable Python versions")
def test_environment_select():
    conf = config.Config()
    conf.environment_type = "conda"
    conf.pythons = ["2.7", "3.5"]
    conf.matrix = {
        "six": ["1.10"],
    }
    conf.include = [
        {'environment_type': 'conda', 'python': '1.9'}
    ]

    # Check default environment config
    environments = list(environment.get_environments(conf, None))
    items = sorted([(env.tool_name, env.python) for env in environments])
    assert items == [('conda', '1.9'), ('conda', '2.7'), ('conda', '3.5')]

    if HAS_VIRTUALENV:
        # Virtualenv plugin fails on initialization if not available,
        # so these tests pass only if virtualenv is present

        conf.pythons = [PYTHON_VER1]

        # Check default python specifiers
        environments = list(environment.get_environments(conf, ["conda", "virtualenv"]))
        items = sorted((env.tool_name, env.python) for env in environments)
        assert items == [('conda', '1.9'), ('conda', PYTHON_VER1), ('virtualenv', PYTHON_VER1)]

        # Check specific python specifiers
        environments = list(environment.get_environments(conf,
                                                         ["conda:3.5",
                                                          "virtualenv:" + PYTHON_VER1]))
        items = sorted((env.tool_name, env.python) for env in environments)
        assert items == [('conda', '3.5'), ('virtualenv', PYTHON_VER1)]

    # Check same specifier
    environments = list(environment.get_environments(conf, ["existing:same", ":same", "existing"]))
    items = [env.tool_name for env in environments]
    assert items == ['existing', 'existing', 'existing']

    # Check autodetect existing
    executable = os.path.relpath(os.path.abspath(sys.executable))
    environments = list(environment.get_environments(conf, ["existing",
                                                            ":same",
                                                            ":" + executable]))
    assert len(environments) == 3
    for env in environments:
        assert env.tool_name == "existing"
        assert env.python == "{0[0]}.{0[1]}".format(sys.version_info)
        assert os.path.normcase(
            os.path.abspath(env._executable)
        ) == os.path.normcase(os.path.abspath(sys.executable))

    # Select by environment name
    conf.pythons = ["2.7"]
    environments = list(environment.get_environments(conf, ["conda-py2.7-six1.10"]))
    assert len(environments) == 1
    assert environments[0].python == "2.7"
    assert environments[0].tool_name == "conda"
    assert environments[0].requirements == {'six': '1.10'}

    # Check interaction with exclude
    conf.exclude = [{'environment_type': "conda"}]
    environments = list(environment.get_environments(conf, ["conda-py2.7-six1.10"]))
    assert len(environments) == 0

    conf.exclude = [{'environment_type': 'matches nothing'}]
    environments = list(environment.get_environments(conf, ["conda-py2.7-six1.10"]))
    assert len(environments) == 1


@pytest.mark.skipif(not (HAS_PYTHON_VER2 or HAS_CONDA),
                    reason="Requires two usable Python versions")
def test_environment_select_autodetect():
    conf = config.Config()
    conf.environment_type = "conda"
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {
        "six": ["1.10"],
    }

    # Check autodetect
    environments = list(environment.get_environments(conf, [":" + PYTHON_VER1]))
    assert len(environments) == 1
    assert environments[0].python == PYTHON_VER1
    assert environments[0].tool_name in ("virtualenv", "conda")

    # Check interaction with exclude
    conf.exclude = [{'environment_type': 'matches nothing'}]
    environments = list(environment.get_environments(conf, [":" + PYTHON_VER1]))
    assert len(environments) == 1

    conf.exclude = [{'environment_type': 'virtualenv|conda'}]
    environments = list(environment.get_environments(conf, [":" + PYTHON_VER1]))
    assert len(environments) == 1

    conf.exclude = [{'environment_type': 'conda'}]
    environments = list(environment.get_environments(conf, ["conda:" + PYTHON_VER1]))
    assert len(environments) == 1

@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda")
def test_matrix_empty():
    conf = config.Config()
    conf.environment_type = ""
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}

    # Check default environment config
    environments = list(environment.get_environments(conf, None))
    items = [env.python for env in environments]
    assert items == [PYTHON_VER1]


@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda")
def test_matrix_existing():
    conf = config.Config()
    conf.environment_type = "existing"
    conf.pythons = ["same"]
    conf.matrix = {'foo': ['a', 'b'], 'bar': ['c', 'd']}

    # ExistingEnvironment should ignore the matrix
    environments = list(environment.get_environments(conf, None))
    items = [(env.tool_name, tuple(env.requirements.keys())) for env in environments]
    assert items == [('existing', ())]

    conf.exclude = {'environment_type': '.*'}
    environments = list(environment.get_environments(conf, None))
    items = [(env.tool_name, tuple(env.requirements.keys())) for env in environments]
    assert items == [('existing', ())]


# environment.yml should respect the specified order
# of channels when adding packages
@pytest.mark.skipif((not HAS_CONDA),
                    reason="Requires conda and conda-build")
@pytest.mark.parametrize("channel_list,expected_channel", [
    (["defaults", "conda-forge"], "pkgs/main"),
    (["conda-forge", "defaults"], "conda-forge"),
])
def test_conda_channel_addition(tmpdir,
                                channel_list,
                                expected_channel):
    # test that we can add conda channels to environments
    # and that we respect the specified priority order
    # of channels
    conf = config.Config()
    conf.env_dir = str(tmpdir.join("env"))
    conf.environment_type = "conda"
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}
    # these have to be valid channels
    # available for online access
    conf.conda_channels = channel_list
    environments = list(environment.get_environments(conf, None))

    # should have one environment per Python version
    assert len(environments) == 1

    # create the environments
    for env in environments:
        env.create()
        # generate JSON output from conda list
        # and parse to verify added channels
        # for current env
        # (conda info would be more direct, but
        # seems to reflect contents of condarc file,
        # which we are intentionally trying not to modify)
        conda = util.which('conda')
        print("\n**conda being used:", conda)
        out_str = str(util.check_output([conda,
                                        'list',
                                         '-p',
                                         os.path.normpath(env._path),
                                         '--json']))
        json_package_list = json.loads(out_str)
        print(json_package_list)
        for installed_package in json_package_list:
            # check only explicitly installed packages
            if installed_package['name'] not in ('python',):
                continue
            print(installed_package)
            assert installed_package['channel'] == expected_channel


@pytest.mark.skipif(not (HAS_PYPY and HAS_VIRTUALENV), reason="Requires pypy and virtualenv")
def test_pypy_virtualenv(tmpdir):
    # test that we can setup a pypy environment
    conf = config.Config()

    conf.env_dir = str(tmpdir.join("env"))

    conf.environment_type = "virtualenv"
    conf.pythons = ["pypy"]
    conf.matrix = {}
    environments = list(environment.get_environments(conf, None))

    assert len(environments) == 1

    for env in environments:
        env.create()
        output = env.run(['-c', 'import sys; print(sys.pypy_version_info)'])
        assert "(major=" in output


@pytest.mark.skipif((not HAS_CONDA), reason="Requires conda")
def test_environment_name_sanitization():
    conf = config.Config()
    conf.environment_type = "conda"
    conf.pythons = ["3.5"]
    conf.matrix = {
        "pip+git+http://github.com/space-telescope/asv.git": [],
    }

    # Check name sanitization
    environments = list(environment.get_environments(conf, []))
    assert len(environments) == 1
    assert environments[0].name == "conda-py3.5-pip+git+http___github.com_space-telescope_asv.git"


@pytest.mark.parametrize("environment_type", [
    pytest.param("conda",
                 marks=pytest.mark.skipif(not HAS_CONDA, reason="needs conda and conda-build")),
    pytest.param("virtualenv",
                 marks=pytest.mark.skipif(not (HAS_PYTHON_VER2 and HAS_VIRTUALENV),
                                          reason="needs virtualenv and python 3.7"))
])
def test_environment_environ_path(environment_type, tmpdir, monkeypatch):
    # Check that virtualenv binary dirs are in the PATH
    conf = config.Config()
    conf.env_dir = str(tmpdir.join("env"))
    conf.environment_type = environment_type
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}

    env, = environment.get_environments(conf, [])
    env.create()
    output = env.run(['-c', 'import os; print(os.environ["PATH"])'])
    paths = output.strip().split(os.pathsep)
    assert os.path.commonprefix([paths[0], conf.env_dir]) == conf.env_dir

    # Check user-site directory is not in sys.path
    output = env.run(['-c', 'import site; print(site.ENABLE_USER_SITE)'])
    usersite_in_syspath = output.strip()
    assert usersite_in_syspath == "False"

    # Check PYTHONPATH is ignored
    monkeypatch.setenv(str('PYTHONPATH'), str(tmpdir))
    output = env.run(['-c', 'import os; print(os.environ.get("PYTHONPATH", ""))'])
    assert output.strip() == ""

    monkeypatch.setenv(str('ASV_PYTHONPATH'), str("Hello python path"))
    output = env.run(['-c', 'import os; print(os.environ["PYTHONPATH"])'])
    assert output.strip() == "Hello python path"



@pytest.mark.skipif(not (HAS_PYTHON_VER2 or HAS_CONDA),
                    reason="Requires two usable Python versions")
def test_build_isolation(tmpdir):
    # build should not fail with build_cache on projects that have pyproject.toml
    tmpdir = str(tmpdir)

    # Create installable repository with pyproject.toml in it
    dvcs = generate_test_repo(tmpdir, [0], dvcs_type='git')
    fn = os.path.join(dvcs.path, 'pyproject.toml')
    with open(fn, 'w') as f:
        f.write('[build-system]\n'
                'requires = ["wheel", "setuptools"]')
    dvcs.add(fn)
    dvcs.commit("Add pyproject.toml")

    commit_hash = dvcs.get_hash(f"{util.git_default_branch()}")

    # Setup config
    conf = config.Config()
    conf.env_dir = os.path.join(tmpdir, "env")
    conf.pythons = [PYTHON_VER1]
    conf.matrix = {}
    conf.repo = os.path.abspath(dvcs.path)
    conf.build_cache_size = 8

    repo = get_repo(conf)

    env = list(environment.get_environments(conf, None))[0]
    env.create()

    # Project installation should succeed
    env.install_project(conf, repo, commit_hash)


@pytest.mark.skipif(tools.HAS_PYPY, reason="Flaky on pypy")
def test_custom_commands(tmpdir):
    # check custom install/uninstall/build commands work
    tmpdir = str(tmpdir)

    dvcs = generate_test_repo(tmpdir, [0], dvcs_type='git')

    build_py = os.path.abspath(os.path.join(tmpdir, 'build.py'))
    install_py = os.path.abspath(os.path.join(tmpdir, 'install.py'))
    uninstall_py = os.path.abspath(os.path.join(tmpdir, 'uninstall.py'))

    conf = config.Config()
    conf.env_dir = os.path.join(tmpdir, "env")
    conf.pythons = [PYTHON_VER1]
    conf.repo = os.path.abspath(dvcs.path)
    conf.matrix = {}
    conf.build_cache_size = 0

    conf.build_command = [f"python {quote(build_py)} {{build_cache_dir}}"]
    conf.install_command = [f"python {quote(install_py)} {{env_dir}} {{build_cache_dir}}"]
    conf.uninstall_command = [f"python {quote(uninstall_py)} {{env_dir}}"]

    with open(build_py, 'wb') as f:
        f.write(b"import os, sys\n"
                b"assert sys.argv[1] == os.environ['ASV_BUILD_CACHE_DIR']\n"
                b"f = open(os.path.join(os.environ['ASV_BUILD_CACHE_DIR'], 'cached'), 'wb')\n"
                b"f.write(b'data')\n"
                b"f.close()\n")

    with open(install_py, 'wb') as f:
        f.write(b"import os, sys, shutil\n"
                b"assert sys.argv[1] == os.environ['ASV_ENV_DIR']\n"
                b"assert sys.argv[2] == os.environ['ASV_BUILD_CACHE_DIR']\n"
                b"shutil.copyfile(os.path.join(os.environ['ASV_BUILD_CACHE_DIR'], 'cached'),\n"
                b"                os.path.join(os.environ['ASV_ENV_DIR'], 'installed'))\n")

    with open(uninstall_py, 'wb') as f:
        f.write(b"import os, sys\n"
                b"assert sys.argv[1] == os.environ['ASV_ENV_DIR']\n"
                b"fn = os.path.join(os.environ['ASV_ENV_DIR'], 'installed')\n"
                b"if os.path.isfile(fn): os.unlink(fn)\n")

    def get_env():
        env = list(environment.get_environments(conf, None))[0]
        env.create()
        return env

    env = get_env()
    repo = get_repo(conf)
    commit_hash = dvcs.get_branch_hashes()[0]

    cache_dir = os.path.join(env._path, 'asv-build-cache')
    cache_file = os.path.join(cache_dir, commit_hash, 'cached')
    install_file = os.path.join(env._path, 'installed')

    # Project installation should succeed with cache size 0,
    # and not leave cache files around
    env.install_project(conf, repo, commit_hash)
    assert os.path.isfile(install_file)
    assert not os.listdir(cache_dir)
    env._set_installed_commit_hash(None)

    # It should succed with nonzero cache size
    conf.build_cache_size = 1
    env = get_env()
    env.install_project(conf, repo, commit_hash)

    assert os.path.isfile(cache_file)
    assert os.path.isfile(install_file)

    # Explicitly check uninstall works
    env._uninstall_project()
    assert os.path.isfile(cache_file)
    assert not os.path.isfile(install_file)

    # Check reinstall uses cache and doesn't call build command
    conf.build_command = ['python -c "import sys; sys.exit(1)"']
    env = get_env()
    env.install_project(conf, repo, commit_hash)

    assert os.path.isfile(install_file)
    assert os.path.isfile(cache_file)

    # Bad install command should cause a failure
    conf.install_command = ['python -c "import sys; sys.exit(1)"']
    env = get_env()
    with pytest.raises(util.ProcessError):
        env.install_project(conf, repo, commit_hash)


def test_installed_commit_hash(tmpdir):
    tmpdir = str(tmpdir)

    dvcs = generate_test_repo(tmpdir, [0], dvcs_type='git')
    commit_hash = dvcs.get_branch_hashes()[0]

    conf = config.Config()
    conf.env_dir = os.path.join(tmpdir, "env")
    conf.pythons = [PYTHON_VER1]
    conf.repo = os.path.abspath(dvcs.path)
    conf.matrix = {}
    conf.build_cache_size = 0

    repo = get_repo(conf)

    def get_env():
        return list(environment.get_environments(conf, None))[0]

    env = get_env()
    env.create()

    # Check updating installed_commit_hash
    assert env.installed_commit_hash is None
    assert env._global_env_vars.get('ASV_COMMIT') is None
    env.install_project(conf, repo, commit_hash)
    assert env.installed_commit_hash == commit_hash
    assert env._global_env_vars.get('ASV_COMMIT') == commit_hash

    env = get_env()
    assert env.installed_commit_hash == commit_hash
    assert env._global_env_vars.get('ASV_COMMIT') == commit_hash

    # Configuration change results to reinstall
    env._project = "something"
    assert env.installed_commit_hash is None

    # Uninstall resets hash (but not ASV_COMMIT)
    env = get_env()
    env._uninstall_project()
    assert env.installed_commit_hash is None
    assert env._global_env_vars.get('ASV_COMMIT') is not None

    env = get_env()
    assert env.installed_commit_hash is None
    assert env._global_env_vars.get('ASV_COMMIT') is None


def test_install_success(tmpdir):
    # Check that install_project really installs the package. (gh-805)
    # This may fail if pip in install_command e.g. gets confused by an .egg-info
    # directory in its cwd to think the package is already installed.
    tmpdir = str(tmpdir)

    dvcs = generate_test_repo(tmpdir, [0], dvcs_type='git')
    commit_hash = dvcs.get_branch_hashes()[0]

    conf = config.Config()
    conf.env_dir = os.path.join(tmpdir, "env")
    conf.pythons = [PYTHON_VER1]
    conf.repo = os.path.abspath(dvcs.path)
    conf.matrix = {}
    conf.build_cache_size = 0

    repo = get_repo(conf)

    env = list(environment.get_environments(conf, None))[0]
    env.create()
    env.install_project(conf, repo, commit_hash)

    env.run(['-c', 'import asv_test_repo as t, sys; sys.exit(0 if t.dummy_value == 0 else 1)'])


def test_install_env_matrix_values(tmpdir):
    tmpdir = str(tmpdir)

    dvcs = generate_test_repo(tmpdir, [0], dvcs_type='git')
    commit_hash = dvcs.get_branch_hashes()[0]

    conf = config.Config()
    conf.env_dir = os.path.join(tmpdir, "env")
    conf.pythons = [PYTHON_VER1]
    conf.repo = os.path.abspath(dvcs.path)
    conf.matrix = {'env': {'SOME_ASV_TEST_BUILD_VALUE': '1'},
                   'env_nobuild': {'SOME_ASV_TEST_NON_BUILD_VALUE': '1'}}

    repo = get_repo(conf)

    env = list(environment.get_environments(conf, None))[0]
    env.create()
    env.install_project(conf, repo, commit_hash)

    env.run(['-c',
             'import asv_test_repo.build_time_env as t, sys; '
             'sys.exit(0 if t.env["SOME_ASV_TEST_BUILD_VALUE"] == "1" else 1)'])

    env.run(['-c',
             'import asv_test_repo.build_time_env as t, sys; '
             'sys.exit(0 if "SOME_ASV_TEST_NON_BUILD_VALUE" not in t.env else 1)'])


def test_environment_env_matrix():
    # (build_vars, non_build_vars, environ_count, build_count)
    configs = [
        ({}, {}, 1, 1),
        ({"var1": ["val1"]}, {}, 1, 1),
        ({"var1": ["val1", "val2", "val3"]}, {}, 3, 3),
        ({"var1": ["val1", "val2"], "var2": ['val3', 'val4']}, {}, 4, 4),
        ({"var1": ["val1", "val2"], "var2": ['val3', None]}, {}, 4, 4),
        ({"var1": ["val1", "val2"]}, {"var2": ['val3', None]}, 4, 2),
        ({"var1": ["val1", "val2"], "var2": ['val3', 'val4']},
         {"var3": ['val5', None]}, 8, 4),
    ]

    for build_vars, non_build_vars, environ_count, build_count in configs:
        conf = config.Config()

        conf.matrix = {
            "env": build_vars,
            "env_nobuild": non_build_vars,
        }
        environments = list(environment.get_environments(conf, None))

        assert len(environments) == environ_count
        assert len(set(e.dir_name for e in environments)) == build_count


def test__parse_matrix():
    cases = [
        ({"env": {"A": "B"}, "env_nobuild": {"C": None}, "req": {"foo": ["9"]}},
         {("env", "A"): "B", ("env_nobuild", "C"): None, ("req", "foo"): ["9"]})
    ]
    for rule, expected in cases:
        parsed = environment._parse_matrix(rule)
        assert parsed == expected


def test__parse_matrix_invalid():
    cases = [
        {"env": "1", "req": "1", "foo": "1"},
    ]
    for m in cases:
        with pytest.raises(util.UserError):
            environment._parse_matrix(m)


def test__parse_matrix_legacy():
    cases = [
        ({"foo": "1", "bar": ["2", "3"]},
         {("req", "foo"): "1", ("req", "bar"): ["2", "3"]})
    ]
    for m, expected in cases:
        parsed = environment._parse_matrix(m)
        assert parsed == expected


def test__parse_exclude_include_rule():
    cases = [
        ({"python": "2.6", "environment_type": "conda", "sys_platform": "123",
          "env": {"A": "B"}, "env_nobuild": {"C": "D"}, "req": {"foo": "9"}},
         {("python", None): "2.6", ("environment_type", None): "conda",
          ("sys_platform", None): "123", ("env", "A"): "B",
          ("env_nobuild", "C"): "D", ("req", "foo"): "9"})
    ]
    for rule, expected in cases:
        parsed = environment._parse_exclude_include_rule(rule)
        assert parsed == expected


def test__parse_exclude_include_rule_invalid():
    cases = [
        {"python": "2.6", "environment_type": "conda", "sys_platform": "123",
         "env": {"A": "B"}, "env_nobuild": {"C": "D"}, "req": {"foo": "9"}, "foo": "9"},
    ]
    for rule in cases:
        with pytest.raises(util.UserError):
            environment._parse_exclude_include_rule(rule)


def test__parse_matrix_entries():
    entries = {("python", None): "2.6", ("env", "A"): "B",
               ("env_nobuild", "C"): "D", ("req", "foo"): "9"}
    python, requirements, tagged_env_vars = environment._parse_matrix_entries(entries)
    assert python == "2.6"
    assert requirements == {"foo": "9"}
    assert tagged_env_vars == {("build", "A"): "B", ("nobuild", "C"): "D"}
