from __future__ import annotations

import logging
from contextlib import contextmanager
from subprocess import Popen
from typing import TYPE_CHECKING

from virtualenv.seed.embed.base_embed import BaseEmbed
from virtualenv.seed.wheels import Version, get_wheel, pip_wheel_env_run
from virtualenv.util.subprocess import LogCmd

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

    from virtualenv.config.cli.parser import VirtualEnvOptions
    from virtualenv.create.creator import Creator

LOGGER = logging.getLogger(__name__)


class PipInvoke(BaseEmbed):
    def __init__(self, options: VirtualEnvOptions) -> None:
        super().__init__(options)

    def run(self, creator: Creator) -> None:
        if not self.enabled:
            return
        for_py_version = creator.interpreter.version_release_str
        with self.get_pip_install_cmd(creator.exe, for_py_version) as cmd:
            env = pip_wheel_env_run(self.extra_search_dir, self.app_data, self.env)
            self._execute(cmd, env)

    @staticmethod
    def _execute(cmd: list[str], env: dict[str, str]) -> Popen[bytes]:
        LOGGER.debug("pip seed by running: %s", LogCmd(cmd, env))
        process = Popen(cmd, env=env)
        process.communicate()
        if process.returncode != 0:
            msg = f"failed seed with code {process.returncode}"
            raise RuntimeError(msg)
        return process

    @contextmanager
    def get_pip_install_cmd(self, exe: Path, for_py_version: str) -> Generator[list[str], None, None]:
        cmd = [
            str(exe),
            "-m",
            "pip",
            "-q",
            "install",
            "--only-binary",
            ":all:",
            "--disable-pip-version-check",
            "--ignore-installed",
        ]
        if not self.download:
            cmd.append("--no-index")
        folders = set()
        for dist, version in self.distribution_to_versions().items():
            wheel = get_wheel(
                distribution=dist,
                version=version,
                for_py_version=for_py_version,
                search_dirs=self.extra_search_dir,
                download=False,
                app_data=self.app_data,
                do_periodic_update=self.periodic_update,
                env=self.env,
            )
            if wheel is None:
                msg = f"could not get wheel for distribution {dist}"
                raise RuntimeError(msg)
            folders.add(str(wheel.path.parent))
            cmd.append(Version.as_pip_req(dist, wheel.version))
        for folder in sorted(folders):
            cmd.extend(["--find-links", str(folder)])
        yield cmd


__all__ = [
    "PipInvoke",
]
