from __future__ import annotations

import errno
import logging
import os
import sys
from timeit import default_timer
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from virtualenv.config.cli.parser import VirtualEnvOptions
    from virtualenv.run.session import Session

LOGGER = logging.getLogger(__name__)


def run(
    args: list[str] | None = None, options: VirtualEnvOptions | None = None, env: MutableMapping[str, str] | None = None
) -> None:
    env = os.environ if env is None else env
    start = default_timer()
    from virtualenv.run import cli_run  # ruff:ignore[import-outside-top-level]
    from virtualenv.util.error import ProcessCallFailedError  # ruff:ignore[import-outside-top-level]

    if args is None:
        args = sys.argv[1:]
    try:
        session = cli_run(args, options, env=env)
        LOGGER.warning(LogSession(session, start))
    except ProcessCallFailedError as exception:
        print(f"subprocess call failed for {exception.cmd} with code {exception.code}")  # ruff:ignore[print]
        print(exception.out, file=sys.stdout, end="")  # ruff:ignore[print]
        print(exception.err, file=sys.stderr, end="")  # ruff:ignore[print]
        raise SystemExit(exception.code)  # ruff:ignore[raise-without-from-inside-except]
    except OSError as exception:
        if exception.errno == errno.EMFILE:
            print(  # ruff:ignore[print]
                "OSError: [Errno 24] Too many open files. You may need to increase your OS open files limit.\n"
                "  On macOS/Linux, try 'ulimit -n 2048'.\n"
                "  For Windows, this is not a common issue, but you can try to close some applications.",
                file=sys.stderr,
            )
        raise


class LogSession:
    def __init__(self, session: Session, start: float) -> None:
        self.session = session
        self.start = start

    def __str__(self) -> str:
        spec = self.session.creator.interpreter.spec
        elapsed = (default_timer() - self.start) * 1000
        lines = [
            f"created virtual environment {spec} in {elapsed:.0f}ms",
            f"  creator {self.session.creator!s}",
        ]
        if self.session.seeder.enabled:
            lines.append(f"  seeder {self.session.seeder!s}")
            path = self.session.creator.purelib.iterdir()
            packages = sorted("==".join(i.stem.split("-")) for i in path if i.suffix == ".dist-info")
            lines.append(f"    added seed packages: {', '.join(packages)}")

        if self.session.activators:
            lines.append(f"  activators {','.join(i.__class__.__name__ for i in self.session.activators)}")
        return "\n".join(lines)


def run_with_catch(args: list[str] | None = None, env: MutableMapping[str, str] | None = None) -> None:
    from virtualenv.config.cli.parser import VirtualEnvOptions  # ruff:ignore[import-outside-top-level]

    env = os.environ if env is None else env
    options = VirtualEnvOptions()
    try:
        run(args, options, env)
    except (KeyboardInterrupt, SystemExit, Exception) as exception:  # ruff:ignore[blind-except]
        try:
            _exit_for_exception(options, exception)
        finally:
            for handler in LOGGER.handlers:  # force flush of log messages before the trace is printed
                handler.flush()


def _exit_for_exception(options: VirtualEnvOptions, exception: BaseException) -> None:
    if getattr(options, "with_traceback", False):
        raise exception
    if not (isinstance(exception, SystemExit) and exception.code == 0):
        LOGGER.error("%s: %s", type(exception).__name__, exception)
    sys.exit(exception.code if isinstance(exception, SystemExit) else 1)


if __name__ == "__main__":  # pragma: no cov
    run_with_catch()  # pragma: no cov
