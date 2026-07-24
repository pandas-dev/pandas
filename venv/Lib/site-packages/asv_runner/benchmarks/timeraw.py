import os
import re
import subprocess
import sys
import textwrap
from hashlib import sha256

from ._base import _get_first_attr, code_fingerprint
from .time import TimeBenchmark


def _normalize_timeraw_env(env):
    """
    Normalize a timeraw ``env`` mapping once (load time).

    #### Parameters
    **env** (`dict` or `None`)
    : Extra environment variables for the timed subprocess, or ``None``.

    #### Returns
    **out** (`dict` or `None`)
    : Mapping of ``str`` keys to ``str`` values, or ``None`` if unused/empty.

    #### Raises
    **TypeError**
    : If ``env`` is not a ``dict``, or any value is ``None`` (cannot be a real
    env entry; refusing to coerce to the string ``"None"``).
    """
    if env is None:
        return None
    if not isinstance(env, dict):
        raise TypeError(
            f"timeraw benchmark attribute 'env' must be a dict, got {type(env)!r}"
        )
    if not env:
        return None
    out = {}
    for key, value in env.items():
        if value is None:
            raise TypeError(
                f"timeraw env values must not be None (key {key!r}); "
                "omit the key instead"
            )
        out[str(key)] = str(value)
    return out


def _env_fingerprint(env):
    """
    Stable, collision-resistant string for ``env`` in benchmark version identity.

    Uses ``repr(sorted(items))`` so values containing newlines or ``=`` cannot
    alias another mapping (e.g. ``{"A": "1\\nB=2"}`` vs ``{"A": "1", "B": "2"}``).
    """
    if not env:
        return ""
    # Sorted (key, value) pairs; repr is unambiguous for str pairs on 3.7+.
    return repr(sorted(env.items()))


class _SeparateProcessTimer:
    """
    Timer that runs a statement in a separate Python process via ``timeit``.

    **env** must already be normalized (or ``None``); not re-validated here.
    """

    subprocess_tmpl = textwrap.dedent(
        '''
        from __future__ import print_function
        from timeit import timeit, default_timer as timer
        print(repr(timeit(stmt="""{stmt}""", setup="""{setup}""",
                    number={number}, timer=timer)))
    '''
    ).strip()

    def __init__(self, func, env=None):
        self.func = func
        self.env = env

    def _child_environ(self):
        """Explicit child environ: parent mapping plus optional overrides."""
        child = dict(os.environ)
        if self.env:
            child.update(self.env)
        return child

    def timeit(self, number):
        stmt = self.func()
        if isinstance(stmt, tuple):
            stmt, setup = stmt
        else:
            setup = ""
        stmt = textwrap.dedent(stmt)
        setup = textwrap.dedent(setup)
        stmt = stmt.replace(r'"""', r"\"\"\"")
        setup = setup.replace(r'"""', r"\"\"\"")

        code = self.subprocess_tmpl.format(stmt=stmt, setup=setup, number=number)

        evaler = textwrap.dedent(
            """
            import sys
            code = sys.stdin.read()
            exec(code)
            """
        )

        proc = subprocess.Popen(
            [sys.executable, "-c", evaler],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=self._child_environ(),
        )
        stdout, stderr = proc.communicate(input=code.encode("utf-8"))
        # Forward timed-process stderr to this process so asv can record it
        # (e.g. timeraw_count writes markers for number*repeat checks).
        if stderr:
            sys.stderr.write(stderr.decode("utf-8", errors="replace"))
            sys.stderr.flush()
        if proc.returncode != 0:
            raise RuntimeError(f"Subprocess failed: {stderr.decode()}")

        return float(stdout.decode("utf-8").strip())


class TimerawBenchmark(TimeBenchmark):
    """
    Timing benchmark executed once per sample in a separate process.

    Set ``env = {"KEY": "value"}`` on the function or class (or use
    ``@benchmark(env={...})``) to inject variables into the **timed** child
    (airspeed-velocity/asv#1471). Default ``version`` includes a fingerprint of
    that mapping so env-only changes invalidate results.
    """

    name_regex = re.compile("^(Timeraw[A-Z_].+)|(timeraw_.+)$")

    def __init__(self, name, func, attr_sources):
        TimeBenchmark.__init__(self, name, func, attr_sources)
        explicit_version = _get_first_attr(attr_sources, "version", None)
        if explicit_version is None:
            env_fp = _env_fingerprint(
                _normalize_timeraw_env(_get_first_attr(attr_sources, "env", None))
            )
            if env_fp:
                payload = self.code + "\n# timeraw_env\n" + env_fp
                self.version = sha256(payload.encode("utf-8")).hexdigest()
                token = code_fingerprint(payload)
                self.version_alts = (token,) if token != self.version else ()

    def _load_vars(self):
        TimeBenchmark._load_vars(self)
        self.number = int(_get_first_attr(self._attr_sources, "number", 1))
        self._timeraw_env = _normalize_timeraw_env(
            _get_first_attr(self._attr_sources, "env", None)
        )
        del self.timer

    def _get_timer(self, *param):
        if param:

            def func():
                return self.func(*param)

        else:
            func = self.func
        return _SeparateProcessTimer(func, env=self._timeraw_env)

    def do_profile(self, filename=None):
        raise ValueError("Raw timing benchmarks cannot be profiled")


export_as_benchmark = [TimerawBenchmark]
