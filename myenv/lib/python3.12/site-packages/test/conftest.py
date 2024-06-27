import contextlib
import sys
import os
import shutil
import textwrap
from os.path import abspath, dirname, join, relpath

import pytest
import selenium

from asv import config, environment, repo, step_detect, util
from asv.repo import get_repo
from asv.step_detect import L1Dist

from . import tools
from .test_benchmarks import ASV_CONF_JSON, BENCHMARK_DIR
from .test_web import _rebuild_basic_html
from .tools import (DUMMY1_VERSION, DUMMY2_VERSIONS, HAS_CONDA, PYTHON_VER1, PYTHON_VER2,
                    WAIT_TIME, WIN, _build_dummy_wheels, locked_cache_dir, run_asv_with_conf)

try:
    import hglib
except ImportError:
    hglib = None

try:
    from asv import _rangemedian
    HAVE_RANGEMEDIAN = True
except ImportError:
    HAVE_RANGEMEDIAN = False

DUMMY_VALUES = (
    (6, 1),
    (6, 6),
    (6, 6),
)


def pytest_addoption(parser):
    parser.addoption("--webdriver", action="store", default="None",
                     help=("Selenium WebDriver interface to use for running the test. "
                           "Choices: None, PhantomJS, Chrome, Firefox, ChromeHeadless, "
                           "FirefoxHeadless. Alternatively, it can be arbitrary Python code "
                           "with a return statement with selenium.webdriver object, for "
                           "example 'return Chrome()'"))
    parser.addoption(
        "--runflaky", action="store_true", default=False, help="run flaky tests"
    )
    parser.addoption("--environment-type", action="store", default=None,
                     choices=("conda", "virtualenv", "mamba"),
                     help="environment_type to use in tests by default")


def generate_basic_conf(tmpdir,
                        repo_subdir='',
                        values=DUMMY_VALUES,
                        dummy_packages=True,
                        conf_version=1):
    # conf_version allows to generate different configurations with this same function
    assert conf_version in (1, 2)
    tmpdir = str(tmpdir)
    local = abspath(dirname(__file__))
    os.chdir(tmpdir)

    # Use relative paths on purpose since this is what will be in
    # actual config files

    shutil.copytree(os.path.join(local, 'benchmark'), 'benchmark')

    machine_file = join(tmpdir, 'asv-machine.json')

    shutil.copyfile(join(local, 'asv-machine.json'),
                    machine_file)

    # values not in test_dev.py copy
    repo_path = tools.generate_test_repo(tmpdir, values,
                                         subdir=repo_subdir).path

    conf_dict = {
        'env_dir': 'env',
        'benchmark_dir': 'benchmark',
        'results_dir': 'results_workflow',
        'html_dir': 'html',
        'repo': relpath(repo_path),
        'project': 'asv',
        'dvcs': 'git',
        'matrix': {
            "asv-dummy-test-package-1": [None],
            "asv-dummy-test-package-2": tools.DUMMY2_VERSIONS,
        },
    }
    if not dummy_packages:
        conf_dict['matrix'] = {}
    elif conf_version == 2:
        conf_dict['matrix'] = {
            "asv_dummy_test_package_1": [""],
            "asv_dummy_test_package_2": tools.DUMMY2_VERSIONS,
        }
    if repo_subdir:
        conf_dict['repo_subdir'] = repo_subdir

    conf = config.Config.from_json(conf_dict)

    if hasattr(sys, 'pypy_version_info'):
        conf.pythons = ["pypy{0[0]}.{0[1]}".format(sys.version_info)]

    return tmpdir, local, conf, machine_file


def pytest_sessionstart(session):
    _monkeypatch_conda_lock(session.config)

    # Unregister unwanted environment types
    env_type = session.config.getoption('environment_type')
    if env_type is not None:
        import asv.environment
        import asv.util

        for cls in asv.util.iter_subclasses(asv.environment.Environment):
            cls.matches_python_fallback = (cls.tool_name in (env_type, "existing"))


def _monkeypatch_conda_lock(config):
    import filelock

    import asv.plugins.conda
    import asv.util

    @contextlib.contextmanager
    def _conda_lock():
        conda_lock = asv.util.get_multiprocessing_lock("conda_lock")
        with conda_lock, filelock.FileLock(str(path)):
            yield

    path = config.cache.makedir('conda-lock') / 'lock'
    asv.plugins.conda._conda_lock = _conda_lock


@pytest.fixture(params=[
    "git",
    pytest.param("hg", marks=pytest.mark.skipif(hglib is None, reason="needs hglib")),
])
def two_branch_repo_case(request, tmpdir):
    r"""
    This test ensure we follow the first parent in case of merges

    The revision graph looks like this:

        @  Revision 6 (default)
        |
        | o  Revision 5 (stable)
        | |
        | o  Merge master
        |/|
        o |  Revision 4
        | |
        o |  Merge stable
        |\|
        o |  Revision 3
        | |
        | o  Revision 2
        |/
        o  Revision 1

    """
    dvcs_type = request.param
    tmpdir = str(tmpdir)
    if dvcs_type == "git":
        master = f"{util.git_default_branch()}"
    elif dvcs_type == "hg":
        master = "default"
    dvcs = tools.generate_repo_from_ops(tmpdir, dvcs_type, [
        ("commit", 1),
        ("checkout", "stable", master),
        ("commit", 2),
        ("checkout", master),
        ("commit", 3),
        ("merge", "stable"),
        ("commit", 4),
        ("checkout", "stable"),
        ("merge", master, "Merge master"),
        ("commit", 5),
        ("checkout", master),
        ("commit", 6),
    ])

    conf = config.Config()
    conf.branches = [master, "stable"]
    conf.repo = dvcs.path
    conf.project = join(tmpdir, "repo")
    r = repo.get_repo(conf)
    return dvcs, master, r, conf


@pytest.fixture
def basic_conf(tmpdir, dummy_packages):
    return generate_basic_conf(tmpdir)


@pytest.fixture
def basic_conf_2(tmpdir, dummy_packages):
    return generate_basic_conf(tmpdir, conf_version=2)


@pytest.fixture
def basic_conf_with_subdir(tmpdir, dummy_packages):
    return generate_basic_conf(tmpdir, 'some_subdir')


@pytest.fixture
def existing_env_conf(tmpdir):
    tmpdir, local, conf, machine_file = generate_basic_conf(tmpdir)
    conf.environment_type = "existing"
    conf.pythons = ["same"]
    return tmpdir, local, conf, machine_file


@pytest.fixture(scope="session")
def example_results(request):
    with locked_cache_dir(request.config, "example-results") as cache_dir:
        src = abspath(join(dirname(__file__), 'example_results'))
        dst = abspath(join(cache_dir, 'results'))

        if os.path.isdir(dst):
            return dst

        shutil.copytree(src, dst)

        src_machine = join(dirname(__file__), 'asv-machine.json')
        dst_machine = join(cache_dir, 'asv-machine.json')
        shutil.copyfile(src_machine, dst_machine)

        # Convert to current file format
        conf = config.Config.from_json({'results_dir': dst,
                                        'repo': 'none',
                                        'project': 'asv'})
        run_asv_with_conf(conf, 'update', _machine_file=dst_machine)

        return dst


@pytest.fixture(scope="session")
def browser(request, pytestconfig):
    """
    Fixture for Selenium WebDriver browser interface
    """
    driver_str = pytestconfig.getoption('webdriver')

    if driver_str == "None":
        pytest.skip("No webdriver selected for tests (use --webdriver).")

    # Evaluate the options
    def FirefoxHeadless():
        options = selenium.webdriver.FirefoxOptions()
        options.add_argument("-headless")
        return selenium.webdriver.Firefox(options=options)

    def ChromeHeadless():
        options = selenium.webdriver.ChromeOptions()
        options.add_argument('headless')
        return selenium.webdriver.Chrome(options=options)

    ns = {}
    exec("import selenium.webdriver", ns)
    exec("from selenium.webdriver import *", ns)
    ns['FirefoxHeadless'] = FirefoxHeadless
    ns['ChromeHeadless'] = ChromeHeadless

    create_driver = ns.get(driver_str, None)
    if create_driver is None:
        src = "def create_driver():\n"
        src += textwrap.indent(driver_str, "    ")
        exec(src, ns)
        create_driver = ns['create_driver']

    # Create the browser
    browser = create_driver()

    # Set timeouts
    browser.set_page_load_timeout(WAIT_TIME)
    browser.set_script_timeout(WAIT_TIME)

    # Clean up on fixture finalization
    def fin():
        browser.quit()
    request.addfinalizer(fin)

    # Set default time to wait for AJAX requests to complete
    browser.implicitly_wait(WAIT_TIME)

    return browser


@pytest.fixture
def dummy_packages(request, monkeypatch):
    """
    Build dummy wheels for required packages and set PIP_FIND_LINKS + CONDARC
    """
    to_build = [('asv_dummy_test_package_1', DUMMY1_VERSION)]
    to_build += [('asv_dummy_test_package_2', ver) for ver in DUMMY2_VERSIONS]

    tag = [PYTHON_VER1, PYTHON_VER2, to_build, HAS_CONDA]

    with locked_cache_dir(request.config, "asv-wheels", timeout=900, tag=tag) as cache_dir:
        wheel_dir = os.path.abspath(join(str(cache_dir), 'wheels'))

        monkeypatch.setenv(str('PIP_FIND_LINKS'), str('file://' + wheel_dir))

        condarc = join(wheel_dir, 'condarc')
        monkeypatch.setenv(str('CONDARC'), str(condarc))

        if os.path.isdir(wheel_dir):
            return

        tmpdir = join(str(cache_dir), "tmp")
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
        os.makedirs(tmpdir)

        try:
            os.makedirs(wheel_dir)
            _build_dummy_wheels(tmpdir, wheel_dir, to_build, build_conda=HAS_CONDA)
        except Exception:
            shutil.rmtree(wheel_dir)
            raise

        # Conda packages were installed in a local channel
        if not WIN:
            wheel_dir_str = f"file://{wheel_dir}"
        else:
            wheel_dir_str = wheel_dir

        with open(condarc, 'w') as f:
            f.write(f"channels:\n- defaults\n- {wheel_dir_str}")


@pytest.fixture(scope="session")
def basic_html(request):
    with locked_cache_dir(request.config, "asv-test_web-basic_html", timeout=900) as cache_dir:
        tmpdir = join(str(cache_dir), 'cached')
        html_dir, dvcs = _rebuild_basic_html(tmpdir)
        return html_dir, dvcs


@pytest.fixture
def benchmarks_fixture(tmpdir):
    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    shutil.copytree(BENCHMARK_DIR, 'benchmark')

    d = {}
    d.update(ASV_CONF_JSON)
    d['env_dir'] = "env"
    d['benchmark_dir'] = 'benchmark'
    d['repo'] = tools.generate_test_repo(tmpdir, [0]).path
    d['branches'] = ["master"]
    conf = config.Config.from_json(d)

    repo = get_repo(conf)
    envs = list(environment.get_environments(conf, None))
    commit_hash = repo.get_hash_from_name(repo.get_branch_name())

    return conf, repo, envs, commit_hash


@pytest.fixture(params=[
    "git",
    pytest.param("hg", marks=pytest.mark.skipif(hglib is None, reason="needs hglib")),
])
def generate_result_dir(request, tmpdir):
    tmpdir = str(tmpdir)
    dvcs_type = request.param

    def _generate_result_dir(values, commits_without_result=None):
        dvcs = tools.generate_repo_from_ops(
            tmpdir, dvcs_type, [("commit", i) for i in range(len(values))])
        commits = list(reversed(dvcs.get_branch_hashes()))
        commit_values = {}
        commits_without_result = [commits[i] for i in commits_without_result or []]
        for commit, value in zip(commits, values):
            if commit not in commits_without_result:
                commit_values[commit] = value
        conf = tools.generate_result_dir(tmpdir, dvcs, commit_values)
        repo = get_repo(conf)
        return conf, repo, commits
    return _generate_result_dir


@pytest.fixture
def show_fixture(tmpdir, example_results):
    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    conf = config.Config.from_json(
        {'results_dir': example_results,
         'repo': tools.generate_test_repo(tmpdir).path,
         'project': 'asv',
         'environment_type': "shouldn't matter what"})

    return conf


@pytest.fixture(params=[
    "python",
    pytest.param("rangemedian",
                 marks=pytest.mark.skipif(not HAVE_RANGEMEDIAN,
                                          reason="compiled asv._rangemedian required"))
])
def use_rangemedian(request):
    if request.param == "rangemedian":
        assert isinstance(step_detect.get_mu_dist([0], [1]), _rangemedian.RangeMedian)
        return True
    else:
        step_detect._rangemedian = None

        def restore():
            if HAVE_RANGEMEDIAN:
                step_detect._rangemedian = _rangemedian
        request.addfinalizer(restore)

        assert isinstance(step_detect.get_mu_dist([0], [1]), L1Dist)
        return False

def pytest_configure(config):
    config.addinivalue_line("markers",
        "flaky_pypy: Tests that are flaky on pypy.")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runflaky"):
        # --runflaky given in cli: do not skip flaky tests
        return
    skip_flaky = pytest.mark.skip(reason="need --runflaky option to run")
    for item in items:
        if "flaky" in item.keywords:
            item.add_marker(skip_flaky)
