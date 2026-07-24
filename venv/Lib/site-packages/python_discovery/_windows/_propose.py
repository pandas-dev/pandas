from __future__ import annotations

from typing import TYPE_CHECKING

from python_discovery._py_info import PythonInfo
from python_discovery._py_spec import PythonSpec

from ._pep514 import discover_pythons

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping

    from python_discovery._cache import PyInfoCache

_IMPLEMENTATION_BY_ORG: dict[str, str] = {
    "ContinuumAnalytics": "CPython",
    "PythonCore": "CPython",
}


class Pep514PythonInfo(PythonInfo):
    """A Python information acquired from PEP-514."""


def propose_interpreters(
    spec: PythonSpec,
    cache: PyInfoCache | None,
    env: Mapping[str, str],
) -> Generator[PythonInfo, None, None]:
    existing = list(discover_pythons())
    existing.sort(
        key=lambda i: (
            *tuple(-1 if j is None else j for j in i[1:4]),
            1 if i[0] == "PythonCore" else 0,
        ),
        reverse=True,
    )

    for name, major, minor, arch, threaded, exe, _ in existing:
        implementation = _IMPLEMENTATION_BY_ORG.get(name, name)

        skip_pre_filter = implementation.lower() != "cpython"
        registry_spec = PythonSpec("", implementation, major, minor, None, arch, exe, free_threaded=threaded)
        if skip_pre_filter or registry_spec.satisfies(spec):
            interpreter = Pep514PythonInfo.from_exe(exe, cache, env=env, raise_on_error=False)
            if interpreter is not None and interpreter.satisfies(spec, impl_must_match=True):
                yield interpreter


__all__ = [
    "Pep514PythonInfo",
    "propose_interpreters",
]
