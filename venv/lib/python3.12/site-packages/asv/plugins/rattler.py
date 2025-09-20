import asyncio
import os
import sys
from pathlib import Path

from yaml import load

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

from rattler import ChannelPriority, MatchSpec, VirtualPackage, install, solve

from .. import environment, util
from ..console import log


class Rattler(environment.Environment):
    """
    Manage an environment using py-rattler.

    Dependencies are installed using py-rattler. The benchmarked
    project is installed using the build command specified.
    """

    tool_name = "rattler"

    def __init__(self, conf, python, requirements, tagged_env_vars):
        """
        Parameters
        ----------
        conf : Config instance

        python : str
            Version of Python.  Must be of the form "MAJOR.MINOR".

        requirements : dict
            Dictionary mapping a PyPI package name to a version
            identifier string.
        """
        self._python = python
        self._requirements = requirements
        self._channels = conf.conda_channels
        self._environment_file = None

        if conf.conda_environment_file == "IGNORE":
            log.debug(
                "Skipping environment file due to conda_environment_file set to IGNORE"
            )
            self._environment_file = None
        elif not conf.conda_environment_file:
            if (Path("environment.yml")).exists():
                log.debug("Using environment.yml")
                self._environment_file = "environment.yml"
        else:
            if (Path(conf.conda_environment_file)).exists():
                log.debug(f"Using {conf.conda_environment_file}")
                self._environment_file = conf.conda_environment_file
            else:
                log.debug(
                    f"Environment file {conf.conda_environment_file} not found, ignoring"
                )

        super().__init__(conf, python, requirements, tagged_env_vars)
        # Rattler configuration things
        self._pkg_cache = f"{self._env_dir}/pkgs"

        self._channel_priority = ChannelPriority.Strict
        condarc_path = Path(os.getenv("CONDARC", ""))
        if condarc_path.is_file() and os.getenv("ASV_USE_CONDARC"):
            log.debug(f"Loading environment configuration from {condarc_path}")
            with condarc_path.open() as f:
                condarc_data = load(f, Loader=Loader) or {}

            if "channels" in condarc_data:
                self._channels = condarc_data["channels"] + self._channels

            if "channel_priority" in condarc_data:
                priority_str = condarc_data["channel_priority"]
                priority_map = {
                    "strict": ChannelPriority.Strict,
                    "flexible": ChannelPriority.Flexible,
                    "disabled": ChannelPriority.Disabled,
                }
                if priority_str in priority_map:
                    self._channel_priority = priority_map[priority_str]
                    log.debug(f"Set channel priority to {priority_str}")
                else:
                    log.warning(
                        f"Unknown channel_priority '{priority_str}' in .condarc"
                    )

    def _setup(self):
        asyncio.run(self._async_setup())

    async def _async_setup(self):
        log.info(f"Creating environment for {self.name}")

        _args, pip_args = self._get_requirements()
        _pkgs = ["python", "wheel", "pip"]  # baseline, overwritten by env file
        env = dict(os.environ)
        env.update(self.build_env_vars)
        if self._environment_file:
            # For named environments
            env_file_name = self._environment_file
            env_data = load(Path(env_file_name).open(), Loader=Loader)
            _pkgs = [x for x in env_data.get("dependencies", []) if isinstance(x, str)]
            self._channels += [
                x for x in env_data.get("channels", []) if isinstance(x, str)
            ]
            self._channels = list(dict.fromkeys(self._channels).keys())
            # Handle possible pip keys
            pip_maybe = [
                x for x in env_data.get("dependencies", []) if isinstance(x, dict)
            ]
            if len(pip_maybe) == 1:
                try:
                    pip_args += pip_maybe[0]["pip"]
                except KeyError:
                    raise KeyError("Only pip is supported as a secondary key")
        _pkgs += _args
        _pkgs = [util.replace_cpython_version(pkg, self._python) for pkg in _pkgs]
        specs = [MatchSpec(p) for p in _pkgs]
        if hasattr(VirtualPackage, "current"):
            virtual_packages = VirtualPackage.current()
        else:
            virtual_packages = VirtualPackage.detect()

        # Expand the 'defaults' meta-channel as rattler requires concrete
        # channel URLs.
        expanded_channels = []
        for channel in self._channels:
            if channel == "defaults":
                log.debug("Expanding 'defaults' meta-channel")
                expanded_channels.extend(
                    [
                        "https://repo.anaconda.com/pkgs/main",
                        "https://repo.anaconda.com/pkgs/r",
                    ]
                )
                if sys.platform.startswith("win"):
                    expanded_channels.append("https://repo.anaconda.com/pkgs/msys2")
            else:
                expanded_channels.append(channel)
        final_channels = list(dict.fromkeys(expanded_channels).keys())

        solved_records = await solve(
            # Channels to use for solving
            channels=final_channels,
            # The specs to solve for
            specs=specs,
            # Virtual packages define the specifications of the environment
            virtual_packages=virtual_packages,
            channel_priority=self._channel_priority,
        )
        await install(records=solved_records, target_prefix=self._path)
        if pip_args:
            for declaration in pip_args:
                parsed_declaration = util.ParsedPipDeclaration(declaration)
                pip_call = util.construct_pip_call(self._run_pip, parsed_declaration)
                pip_call()

    def _get_requirements(self):
        _args = []
        pip_args = []

        for key, val in {**self._requirements, **self._base_requirements}.items():
            if key.startswith("pip+"):
                pip_args.append(f"{key[4:]} {val}")
            else:
                if val:
                    _args.append(f"{key}={val}")
                else:
                    _args.append(key)

        return _args, pip_args

    def run_executable(self, executable, args, **kwargs):
        return super().run_executable(executable, args, **kwargs)

    def run(self, args, **kwargs):
        log.debug(f"Running '{' '.join(args)}' in {self.name}")
        return self.run_executable("python", args, **kwargs)

    def _run_pip(self, args, **kwargs):
        # Run pip via python -m pip, so that it works on Windows when
        # upgrading pip itself, and avoids shebang length limit on Linux
        return self.run_executable("python", ["-m", "pip"] + list(args), **kwargs)
