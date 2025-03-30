# Licensed under a 3-clause BSD style license - see LICENSE.rst
import sys
import re
import os

from packaging.version import Version

from .. import environment, util
from ..console import log

WIN = (os.name == "nt")


class Virtualenv(environment.Environment):
    """
    Manage an environment using virtualenv.
    """
    tool_name = "virtualenv"

    def __init__(self, conf, python, requirements, tagged_env_vars):
        """
        Parameters
        ----------
        conf : Config instance

        python : str
            Version of Python.  Must be of the form "MAJOR.MINOR".

        executable : str
            Path to Python executable.

        requirements : dict
            Dictionary mapping a PyPI package name to a version
            identifier string.
        """
        executable = Virtualenv._find_python(python)
        if executable is None:
            raise environment.EnvironmentUnavailable(
                f"No executable found for python {python}")

        self._executable = executable
        self._python = python
        self._requirements = requirements
        super(Virtualenv, self).__init__(conf,
                                         python,
                                         requirements,
                                         tagged_env_vars)

        try:
            import virtualenv  # noqa F401 unused, but required to test whether virtualenv is installed or not
        except ImportError:
            raise environment.EnvironmentUnavailable(
                "virtualenv package not installed")

    @staticmethod
    def _find_python(python):
        """Find Python executable for the given Python version"""
        is_pypy = python.startswith("pypy")

        # Parse python specifier
        if is_pypy:
            executable = python
            if python == 'pypy':
                python_version = '2'
            else:
                python_version = python[4:]
        else:
            python_version = python
            executable = f"python{python_version}"

        # Find Python executable on path
        try:
            return util.which(executable)
        except OSError:
            pass

        # Maybe the current one is correct?
        current_is_pypy = hasattr(sys, 'pypy_version_info')
        current_versions = [f'{sys.version_info[0]}',
                            f'{sys.version_info[0]}.{sys.version_info[1]}']

        if is_pypy == current_is_pypy and python_version in current_versions:
            return sys.executable

        return None

    @property
    def name(self):
        """
        Get a name to uniquely identify this environment.
        """
        python = self._python
        if self._python.startswith('pypy'):
            # get_env_name adds py-prefix
            python = python[2:]
        return environment.get_env_name(self.tool_name,
                                        python,
                                        self._requirements,
                                        self._tagged_env_vars)

    @classmethod
    def matches(self, python):
        if not (re.match(r'^[0-9].*$', python) or re.match(r'^pypy[0-9.]*$', python)):
            # The python name should be a version number, or pypy+number
            return False

        try:
            import virtualenv
        except ImportError:
            return False
        else:
            if Version(virtualenv.__version__) == Version('1.11.0'):
                log.warning(
                    "asv is not compatible with virtualenv 1.11 due to a bug in "
                    "setuptools.")
            if Version(virtualenv.__version__) < Version('1.10'):
                log.warning(
                    "If using virtualenv, it much be at least version 1.10")

        executable = Virtualenv._find_python(python)
        return executable is not None

    def _setup(self):
        """
        Setup the environment on disk using virtualenv.
        Then, all of the requirements are installed into
        it using `pip install`.
        """
        env = dict(os.environ)
        env.update(self.build_env_vars)

        log.info(f"Creating virtualenv for {self.name}")
        util.check_call([
            sys.executable,
            "-mvirtualenv",
            "--wheel=bundle",
            "--setuptools=bundle",
            "-p",
            self._executable,
            self._path], env=env)

        log.info(f"Installing requirements for {self.name}")
        self._install_requirements()

    def _install_requirements(self):
        pip_args = ['install', '-v', 'wheel', 'pip>=8']

        env = dict(os.environ)
        env.update(self.build_env_vars)

        self._run_pip(pip_args, env=env)
        pip_args = []

        for key, val in {**self._requirements,
                         **self._base_requirements}.items():
            if key.startswith("pip+"):
                pip_args.append(f"{key[4:]} {val}")
            else:
                pip_args.append(f"{key} {val}")

        for declaration in pip_args:
            parsed_declaration = util.ParsedPipDeclaration(declaration)
            pip_call = util.construct_pip_call(self._run_pip, parsed_declaration)
            pip_call()

    def _run_pip(self, args, **kwargs):
        # Run pip via python -m pip, so that it works on Windows when
        # upgrading pip itself, and avoids shebang length limit on Linux
        return self.run_executable('python', ['-mpip'] + list(args), **kwargs)

    def run(self, args, **kwargs):
        joined_args = ' '.join(args)
        log.debug(f"Running '{joined_args}' in {self.name}")
        return self.run_executable('python', args, **kwargs)
