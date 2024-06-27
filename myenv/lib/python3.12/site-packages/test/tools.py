# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This file contains utilities to generate test repositories.
"""

import datetime
import io
import os
import threading
import time
import tempfile
import textwrap
import sys
import shutil
import subprocess
import platform
import http.server
import importlib
from os.path import abspath, join, dirname, relpath, isdir
from contextlib import contextmanager
from hashlib import sha256

from filelock import FileLock

try:
    import hglib
except ImportError:
    hglib = None

import asv
from asv import util, commands, config, environment, runner
from asv.commands.preview import create_httpd
from asv.repo import get_repo
from asv.results import Results
from asv.plugins.conda import _find_conda

# Two Python versions for testing
PYTHON_VER1, PYTHON_VER2 = '3.7', platform.python_version()

# Installable library versions to use in tests
DUMMY1_VERSION = "0.14"
DUMMY2_VERSIONS = ["0.3.7", "0.3.9"]


WIN = (os.name == "nt")

try:
    util.which('pypy')
    HAS_PYPY = True
except (RuntimeError, OSError):
    HAS_PYPY = hasattr(sys, 'pypy_version_info')


def _check_conda():
    from asv.plugins.conda import _conda_lock
    conda = _find_conda()
    with _conda_lock():
        try:
            subprocess.check_call([conda, 'build', '--version'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            raise RuntimeError("conda-build is missing")


try:
    # Conda can install required Python versions on demand
    _check_conda()
    HAS_CONDA = True
except (RuntimeError, OSError):
    HAS_CONDA = False


try:
    import virtualenv  # noqa F401 checking if installed
    HAS_VIRTUALENV = True
except ImportError:
    HAS_VIRTUALENV = False


try:
    util.which(f'python{PYTHON_VER2}')
    HAS_PYTHON_VER2 = True
except (RuntimeError, OSError):
    HAS_PYTHON_VER2 = False


try:
    from selenium.common.exceptions import TimeoutException
    HAVE_WEBDRIVER = True
except ImportError:
    HAVE_WEBDRIVER = False


def _check_mamba():
    try:
        conda = _find_conda()
        importlib.import_module("libmambapy")
        subprocess.check_call(
            [conda, "build", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return True
    except (ImportError, OSError, subprocess.CalledProcessError, FileNotFoundError):
        return False


HAS_MAMBA = _check_mamba()


WAIT_TIME = 20.0


def get_default_environment_type(conf, python):
    return environment.get_environment_class(conf, python).tool_name


@contextmanager
def locked_cache_dir(config, cache_key, timeout=900, tag=None):
    base_dir = config.cache.makedir(cache_key)

    lockfile = join(str(base_dir), 'lock')
    cache_dir = join(str(base_dir), 'cache')

    lock = FileLock(lockfile)
    lock.acquire(timeout=timeout)
    try:
        # Clear cache dir contents if it was generated with different
        # asv version
        tag_fn = join(str(base_dir), 'tag.json')
        tag_content = [asv.__version__, repr(tag)]
        if os.path.isdir(cache_dir):
            try:
                if util.load_json(tag_fn) != tag_content:
                    raise ValueError()
            except (OSError, ValueError, util.UserError):
                shutil.rmtree(cache_dir)

        if not os.path.isdir(cache_dir):
            os.makedirs(cache_dir)

        yield cache_dir

        util.write_json(tag_fn, tag_content)
    finally:
        lock.release()


def run_asv(*argv, **kwargs):
    parser, subparsers = commands.make_argparser()
    args = parser.parse_args(argv)
    return args.func(args, **kwargs)


def run_asv_with_conf(conf, *argv, **kwargs):
    assert isinstance(conf, config.Config)

    parser, subparsers = commands.make_argparser()
    args = parser.parse_args(argv)

    if sys.version_info[0] >= 3:
        cls = args.func.__self__
    else:
        cls = args.func.im_self

    return cls.run_from_conf_args(conf, args, **kwargs)


# These classes are defined here, rather than using asv/plugins/git.py
# and asv/plugins/mercurial.py since here we need to perform write
# operations to the repository, and the others should be read-only for
# safety.

class Git:
    def __init__(self, path):
        self.path = abspath(path)
        self._git = util.which('git')
        self._fake_date = datetime.datetime.now(datetime.timezone.utc)

    def run_git(self, args, chdir=True, **kwargs):
        if chdir:
            cwd = self.path
        else:
            cwd = None
        kwargs['cwd'] = cwd
        return util.check_output(
            [self._git] + args, **kwargs)

    def init(self):
        self.run_git(['init'])
        self.run_git(['config', 'user.email', 'robot@asv'])
        self.run_git(['config', 'user.name', 'Robotic Swallow'])
        self.run_git(['config', 'protocol.file.allow', 'always'])

    def commit(self, message, date=None):
        if date is None:
            self._fake_date += datetime.timedelta(seconds=1)
            date = self._fake_date

        self.run_git(['commit', '--date', date.isoformat(),
                     '-m', message])

    def tag(self, number):
        self.run_git(['tag', '-a', '-m', f'Tag {number}',
                      f'tag{number}'])

    def add(self, filename):
        self.run_git(['add', relpath(filename, self.path)])

    def checkout(self, branch_name, start_commit=None):
        args = ["checkout"]
        if start_commit is not None:
            args.extend(["-b", branch_name, start_commit])
        else:
            args.append(branch_name)
        self.run_git(args)

    def merge(self, branch_name, commit_message=None):
        self.run_git(["merge", "--no-ff", "--no-commit", "-X", "theirs", branch_name])
        if commit_message is None:
            commit_message = f"Merge {branch_name}"
        self.commit(commit_message)

    def get_hash(self, name):
        return self.run_git(['rev-parse', name]).strip()

    def get_branch_hashes(self, branch=None):
        if branch is None:
            branch = util.git_default_branch()
        return [x.strip() for x in self.run_git(['rev-list', branch]).splitlines()
                if x.strip()]

    def get_commit_message(self, commit_hash):
        return self.run_git(["log", "-n", "1", "--format=%s", commit_hash]).strip()


_hg_config = """
[ui]
username = Robotic Swallow <robot@asv>
"""


class Hg:
    encoding = 'utf-8'

    def __init__(self, path):
        self._fake_date = datetime.datetime.now(datetime.timezone.utc)
        if isinstance(path, bytes):
            path = path.decode('utf-8')
        self.path = abspath(path)
        self._repo = None

    def __del__(self):
        if self._repo is not None:
            self._repo.close()
            self._repo = None

    def init(self):
        hglib.init(self.path)
        with io.open(join(self.path, '.hg', 'hgrc'), 'w', encoding="utf-8") as fd:
            fd.write(_hg_config)
        self._repo = hglib.open(self.path.encode(sys.getfilesystemencoding()),
                                encoding=self.encoding)

    def commit(self, message, date=None):
        if date is None:
            self._fake_date += datetime.timedelta(seconds=1)
            date = self._fake_date
        date = f"{util.datetime_to_timestamp(date)} 0"

        self._repo.commit(message.encode(self.encoding),
                          date=date.encode(self.encoding))

    def tag(self, number):
        self._fake_date += datetime.timedelta(seconds=1)
        date = f"{util.datetime_to_timestamp(self._fake_date)} 0"

        self._repo.tag(
            [f'tag{number}'.encode(self.encoding)],
            message=f"Tag {number}".encode(self.encoding),
            date=date.encode(self.encoding))

    def add(self, filename):
        self._repo.add([filename.encode(sys.getfilesystemencoding())])

    def checkout(self, branch_name, start_commit=None):
        if start_commit is not None:
            self._repo.update(start_commit.encode(self.encoding))
            self._repo.branch(branch_name.encode(self.encoding))
        else:
            self._repo.update(branch_name.encode(self.encoding))

    def merge(self, branch_name, commit_message=None):
        self._repo.merge(branch_name.encode(self.encoding),
                         tool=b"internal:other")
        if commit_message is None:
            commit_message = f"Merge {branch_name}"
        self.commit(commit_message)

    def get_hash(self, name):
        log = self._repo.log(name.encode(self.encoding), limit=1)
        if log:
            return log[0][1].decode(self.encoding)
        return None

    def get_branch_hashes(self, branch=None):
        if branch is None:
            branch = "default"
        log = self._repo.log(f'sort(ancestors({branch}), -rev)'.encode(self.encoding))
        return [entry[1].decode(self.encoding) for entry in log]

    def get_commit_message(self, commit_hash):
        return self._repo.log(commit_hash.encode(self.encoding))[0].desc.decode(self.encoding)


def copy_template(src, dst, dvcs, values):
    for root, dirs, files in os.walk(src):
        for dir in dirs:
            src_path = join(root, dir)
            dst_path = join(dst, relpath(src_path, src))
            if not isdir(dst_path):
                os.makedirs(dst_path)

        for file in files:
            src_path = join(root, file)
            dst_path = join(dst, relpath(src_path, src))

            try:
                with io.open(src_path, 'r', encoding='utf-8') as fd:
                    content = fd.read()
            except UnicodeDecodeError:
                # File is some sort of binary file...  just copy it
                # directly with no template substitution
                with io.open(src_path, 'rb') as fd:
                    content = fd.read()
                with io.open(dst_path, 'wb') as fd:
                    fd.write(content)
            else:
                content = content.format(**values)
                with io.open(dst_path, 'w', encoding='utf-8') as fd:
                    fd.write(content)

            dvcs.add(dst_path)


def generate_test_repo(tmpdir, values=[0], dvcs_type='git',
                       extra_branches=(), subdir=''):
    """
    Generate a test repository

    Parameters
    ----------
    tmpdir
        Repository directory
    values : list
        List of values to substitute in the template
    dvcs_type : {'git', 'hg'}
        What dvcs to use
    extra_branches : list of (start_commit, branch_name, values)
        Additional branches to generate in the repository.
        For branch start commits, use relative references, e.g.,
        the format 'main~10' or 'default~10' works both for Hg
        and Git.
    subdir
        A relative subdirectory inside the repository to copy the
        test project into.

    Returns
    -------
    dvcs : Git or Hg

    """
    if dvcs_type == 'git':
        dvcs_cls = Git
    elif dvcs_type == 'hg':
        dvcs_cls = Hg
    else:
        raise ValueError(f"Unknown dvcs type {dvcs_type}")

    template_path = join(dirname(__file__), 'test_repo_template')

    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    dvcs_path = tempfile.mkdtemp(prefix='test_repo', dir=tmpdir)
    dvcs = dvcs_cls(dvcs_path)
    dvcs.init()

    project_path = os.path.join(dvcs_path, subdir)
    if not os.path.exists(project_path):
        os.makedirs(project_path)

    for i, value in enumerate(values):
        mapping = {
            'version': i,
            'dummy_value': value
        }

        copy_template(template_path, project_path, dvcs, mapping)

        dvcs.commit(f"Revision {i}")
        dvcs.tag(i)

    if extra_branches:
        for start_commit, branch_name, values in extra_branches:
            dvcs.checkout(branch_name, start_commit)
            for i, value in enumerate(values):
                mapping = {
                    'version': f"{i}",
                    'dummy_value': value
                }
                copy_template(template_path, project_path, dvcs, mapping)
                dvcs.commit(f"Revision {branch_name}.{i}")

    return dvcs


def generate_repo_from_ops(tmpdir, dvcs_type, operations):
    if dvcs_type == 'git':
        dvcs_cls = Git
    elif dvcs_type == 'hg':
        dvcs_cls = Hg
    else:
        raise ValueError(f"Unknown dvcs type {dvcs_type}")

    template_path = join(dirname(__file__), 'test_repo_template')

    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    dvcs_path = tempfile.mkdtemp(prefix='test_repo', dir=tmpdir)
    dvcs = dvcs_cls(dvcs_path)
    dvcs.init()

    version = 0
    for op in operations:
        if op[0] == "commit":
            copy_template(template_path, dvcs_path, dvcs, {
                "version": version,
                "dummy_value": op[1],
            })
            version += 1
            dvcs.commit(f"Revision {version}", *op[2:])
        elif op[0] == "checkout":
            dvcs.checkout(*op[1:])
        elif op[0] == "merge":
            dvcs.merge(*op[1:])
        else:
            raise ValueError(f"Unknown dvcs operation {op}")

    return dvcs


def generate_result_dir(tmpdir, dvcs, values, branches=None, updated=None):
    result_dir = join(tmpdir, "results")
    os.makedirs(result_dir)
    html_dir = join(tmpdir, "html")
    machine_dir = join(result_dir, "tarzan")
    os.makedirs(machine_dir)

    if branches is None:
        branches = [None]

    conf = config.Config.from_json({
        'results_dir': result_dir,
        'html_dir': html_dir,
        'repo': dvcs.path,
        'project': 'asv',
        'branches': branches or [None],
    })
    repo = get_repo(conf)

    util.write_json(join(machine_dir, "machine.json"), {
        'machine': 'tarzan',
        'version': 1,
    })

    if updated is None:
        updated = datetime.datetime(
            1970, 1, 1,
            tzinfo = datetime.timezone.utc
        )

    benchmark_version = sha256(os.urandom(16)).hexdigest()

    params = []
    param_names = None
    for commit, value in values.items():
        if isinstance(value, dict):
            params = value["params"]
            value = value["result"]
        else:
            value = [value]
        result = Results({"machine": "tarzan"}, {}, commit,
                         repo.get_date_from_name(commit), "2.7", None, {})
        value = runner.BenchmarkResult(
            result=value,
            samples=[None] * len(value),
            number=[None] * len(value),
            errcode=0,
            stderr='',
            profile=None)
        result.add_result({"name": "time_func", "version": benchmark_version, "params": params},
                          value, started_at=updated, duration=1.0)
        result.save(result_dir)

    if params:
        param_names = [f"param{k}" for k in range(len(params))]

    util.write_json(join(result_dir, "benchmarks.json"), {
        "time_func": {
            "name": "time_func",
            "params": params or [],
            "param_names": param_names or [],
            "version": benchmark_version,
        }
    }, api_version=2)
    return conf


@contextmanager
def preview(base_path):
    """
    Context manager for ASV preview web server. Gives the base URL to use.

    Parameters
    ----------
    base_path : str
        Path to serve files from

    """

    class Handler(http.server.SimpleHTTPRequestHandler):
        def translate_path(self, path):
            # Don't serve from cwd, but from a different directory
            path = http.server.SimpleHTTPRequestHandler.translate_path(self, path)
            path = os.path.join(base_path, os.path.relpath(path, os.getcwd()))
            return util.long_path(path)

    httpd, base_url = create_httpd(Handler)

    def run():
        try:
            httpd.serve_forever()
        except Exception:
            import traceback
            traceback.print_exc()
        finally:
            httpd.server_close()

    thread = threading.Thread(target=run)
    thread.daemon = True
    thread.start()
    try:
        yield base_url
    finally:
        # Stop must be run in a separate thread, because
        # httpd.shutdown blocks until serve_forever returns.  We don't
        # want to block here --- it appears in some environments
        # problems shutting down the server may arise.
        stopper = threading.Thread(target=httpd.shutdown)
        stopper.daemon = True
        stopper.start()
        stopper.join(5.0)


def get_with_retry(browser, url):
    for j in range(2):
        try:
            return browser.get(url)
        except TimeoutException:
            time.sleep(2)

    return browser.get(url)


def _build_dummy_wheels(tmpdir, wheel_dir, to_build, build_conda=False):
    # Build fake wheels for testing

    for name, version in to_build:
        build_dir = join(tmpdir, name + '-' + version)
        os.makedirs(build_dir)

        with open(join(build_dir, 'setup.py'), 'w') as f:
            f.write("from setuptools import setup; "
                    "setup(name='{name}', version='{version}', packages=['{name}'])"
                    "".format(name=name, version=version))
        os.makedirs(join(build_dir, name))
        with open(join(build_dir, name, '__init__.py'), 'w') as f:
            f.write(f"__version__ = '{version}'")

        subprocess.check_call([sys.executable, '-mpip', 'wheel',
                               '--build-option=--universal',
                               '-w', wheel_dir,
                               '.'],
                              cwd=build_dir)

        if build_conda:
            _build_dummy_conda_pkg(name, version, build_dir, wheel_dir)


def _build_dummy_conda_pkg(name, version, build_dir, dst):
    # Build fake conda packages for testing
    from asv.plugins.conda import _conda_lock

    build_dir = os.path.abspath(build_dir)

    with open(join(build_dir, 'meta.yaml'), 'w') as f:
        f.write(textwrap.dedent("""\
        package:
          name: "{name}"
          version: "{version}"
        source:
          path: {build_dir}
        build:
          number: 0
          script: "python -m pip install . --no-deps --ignore-installed "
        requirements:
          host:
            - pip
            - python
          run:
            - python
        about:
          license: BSD
          summary: Dummy test package
        """.format(name=name,
                   version=version,
                   build_dir=util.shlex_quote(build_dir))))

    conda = _find_conda()

    for pyver in [PYTHON_VER1, PYTHON_VER2]:
        with _conda_lock():
            subprocess.check_call([conda, 'build',
                                   '--output-folder=' + dst,
                                   '--no-anaconda-upload',
                                   '--python=' + pyver,
                                   '.'],
                                  cwd=build_dir)
