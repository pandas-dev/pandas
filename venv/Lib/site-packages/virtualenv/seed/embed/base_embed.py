from __future__ import annotations

import logging
from abc import ABC
from argparse import SUPPRESS
from pathlib import Path
from typing import TYPE_CHECKING

from virtualenv.seed.seeder import Seeder
from virtualenv.seed.wheels import Version
from virtualenv.seed.wheels.embed import MIN, OLDEST_SUPPORTED

if TYPE_CHECKING:
    from argparse import ArgumentParser

    from python_discovery import PythonInfo

    from virtualenv.app_data.base import AppData
    from virtualenv.config.cli.parser import VirtualEnvOptions

LOGGER = logging.getLogger(__name__)
PERIODIC_UPDATE_ON_BY_DEFAULT = True


class BaseEmbed(Seeder, ABC):
    def __init__(self, options: VirtualEnvOptions) -> None:
        super().__init__(options, enabled=options.no_seed is False)

        self.download = options.download
        self.extra_search_dir = [i.resolve() for i in options.extra_search_dir if i.exists()]

        self.pip_version = options.pip
        self.setuptools_version = options.setuptools

        # virtualenv no longer bundles wheel; the parsed default stays None (unused) so the
        # warning below fires only when you pass --wheel or --no-wheel
        self.wheel_version = options.wheel or "none"

        self.no_pip = options.no_pip
        self.no_setuptools = options.no_setuptools
        self.app_data = options.app_data
        self.periodic_update = not options.no_periodic_update

        if options.wheel is not None or options.no_wheel:
            LOGGER.warning(
                "DEPRECATION: the --wheel and --no-wheel options do nothing; virtualenv no longer bundles wheel. "
                "They will be removed in a release after 2026-12. Stop passing them.",
            )
        self.no_wheel = True

        if not self.distribution_to_versions():
            self.enabled = False

    @classmethod
    def distributions(cls) -> dict[str, str]:
        return {
            "pip": Version.bundle,
            "setuptools": Version.bundle,
            "wheel": Version.bundle,
        }

    def distribution_to_versions(self) -> dict[str, str]:
        return {
            distribution: getattr(self, f"{distribution}_version")
            for distribution in self.distributions()
            if getattr(self, f"no_{distribution}", None) is False and getattr(self, f"{distribution}_version") != "none"
        }

    @classmethod
    def cannot_seed(cls, interpreter: PythonInfo) -> str | None:
        """Explain why the bundled wheels cannot seed the target Python version.

        The embedded pip/setuptools stopped shipping for Pythons below :data:`OLDEST_SUPPORTED`, so seeding one would
        install an incompatible wheel.

        :param interpreter: the interpreter to be seeded

        :returns: ``None`` when the bundled wheels still support the target, otherwise a message naming the target and
            the remedies

        """
        if interpreter.version_info[:2] >= OLDEST_SUPPORTED:
            return None
        target = f"{interpreter.version_info.major}.{interpreter.version_info.minor}"
        return (
            f"the bundled seeder no longer ships pip/setuptools for Python {target}; the oldest supported target is "
            f"Python {MIN} - pass --no-seed for an empty environment, use a seeder that provides Python {target} "
            f"wheels, or install an older virtualenv release"
        )

    @classmethod
    def add_parser_arguments(cls, parser: ArgumentParser, interpreter: PythonInfo, app_data: AppData) -> None:  # ruff:ignore[unused-class-method-argument]
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            "--no-download",
            "--never-download",
            dest="download",
            action="store_false",
            help=f"pass to disable download of the latest {'/'.join(cls.distributions())} from PyPI",
            default=True,
        )
        group.add_argument(
            "--download",
            dest="download",
            action="store_true",
            help=f"pass to enable download of the latest {'/'.join(cls.distributions())} from PyPI",
            default=False,
        )
        parser.add_argument(
            "--extra-search-dir",
            metavar="d",
            type=Path,
            nargs="+",
            help="a path containing wheels to extend the internal wheel list (can be set 1+ times)",
            default=[],
        )
        for distribution, default in cls.distributions().items():
            help_ = f"version of {distribution} to install as seed: embed, bundle, none or exact version"
            if interpreter.version_info[:2] >= (3, 12) and distribution == "setuptools":
                default = "none"  # ruff:ignore[redefined-loop-name]
            if distribution == "wheel":
                default = None  # ruff:ignore[redefined-loop-name]
                help_ = SUPPRESS
            parser.add_argument(
                f"--{distribution}",
                dest=distribution,
                metavar="version",
                help=help_,
                default=default,
            )
        for distribution in cls.distributions():
            help_ = f"do not install {distribution}"
            if distribution == "wheel":
                help_ = SUPPRESS
            parser.add_argument(
                f"--no-{distribution}",
                dest=f"no_{distribution}",
                action="store_true",
                help=help_,
                default=False,
            )
        parser.add_argument(
            "--no-periodic-update",
            dest="no_periodic_update",
            action="store_true",
            help="disable the periodic (once every 14 days) update of the embedded wheels",
            default=not PERIODIC_UPDATE_ON_BY_DEFAULT,
        )

    def __repr__(self) -> str:
        result = self.__class__.__name__
        result += "("
        if self.extra_search_dir:
            result += f"extra_search_dir={', '.join(str(i) for i in self.extra_search_dir)},"
        result += f"download={self.download},"
        for distribution in self.distributions():
            if getattr(self, f"no_{distribution}", None):
                continue
            version = getattr(self, f"{distribution}_version", None)
            if version == "none":
                continue
            ver = f"={version or 'latest'}"
            result += f" {distribution}{ver},"
        return result[:-1] + ")"


__all__ = [
    "BaseEmbed",
]
