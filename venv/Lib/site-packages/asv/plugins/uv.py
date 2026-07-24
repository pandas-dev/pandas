# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import re

from .. import environment, util
from ..console import log

WIN = os.name == "nt"


class Uv(environment.Environment):
    """
    Manage an environment using 'uv venv'.
    """

    tool_name = "uv"

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

        self._python = python
        self._requirements = requirements
        super().__init__(conf, python, requirements, tagged_env_vars)

        try:
            self._uv_path = util.which("uv")
        except OSError:
            self._uv_path = None

        if not self._uv_path:
            raise environment.EnvironmentUnavailable("uv command not found")

    @property
    def name(self):
        """
        Get a name to uniquely identify this environment.
        """
        python = self._python
        if self._python.startswith('pypy'):
            # get_env_name adds py-prefix
            python = python[2:]
        return environment.get_env_name(
            self.tool_name, python, self._requirements, self._tagged_env_vars
        )

    @classmethod
    def matches(self, python):
        if not (re.match(r'^[0-9].*$', python) or re.match(r'^pypy[0-9.]*$', python)):
            # The python name should be a version number, or pypy+number
            return False

        if not self._uv_path:
            log.warning(
                "asv requires the 'uv' command to be available when using the 'uv' environment_type."
            )
            return False

        return True

    def _setup(self):
        """
        Setup the environment on disk using 'uv venv'.
        Then, all of the requirements are installed into
        it using `pip install`.
        """
        env = dict(os.environ)
        env.update(self.build_env_vars)

        # Adjust the environment variables to use the virtualenv
        self._venv_env_vars = {"VIRTUAL_ENV": self._path}
        if "PATH" in env:
            self._venv_env_vars["PATH"] = f"{self._path}/bin:{env['PATH']}"
        else:
            self._venv_env_vars["PATH"] = f"{self._path}/bin"
        env.update(self._venv_env_vars)

        log.info(f"Creating virtualenv for {self.name}")

        util.check_call(
            [
                'uv',
                'venv',
                f'--python={self._python}',
                '--no-project',
                '--seed',
                self._path,
            ],
            env=env,
        )

        log.info(f"Installing requirements for {self.name}")
        self._install_requirements()

    def _install_requirements(self):
        pip_args = ['install', '-v', 'wheel', 'pip>=8']

        env = dict(os.environ)
        env.update(self.build_env_vars)

        self._run_pip(pip_args, env=env)
        pip_args = []

        for key, val in {**self._requirements, **self._base_requirements}.items():
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
        return self.run_executable('python', ['-m', 'pip'] + list(args), **kwargs)

    def run(self, args, **kwargs):
        joined_args = ' '.join(args)
        log.debug(f"Running '{joined_args}' in {self.name}")
        return self.run_executable('python', args, **kwargs)
