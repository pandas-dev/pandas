# This file must have the same content for mypyc/build_setup.py and lib-rt/build_setup.py,
# it exists to work around absence of support for per-file compile flags in setuptools.
# The version in mypyc/ is the source of truth, and should be copied to lib-rt if modified.

import os
import platform
import sys

try:
    # Import setuptools so that it monkey-patch overrides distutils
    import setuptools  # noqa: F401
except ImportError:
    pass

if sys.version_info >= (3, 12):
    # From setuptools' monkeypatch
    from distutils import ccompiler  # type: ignore[import-not-found]
else:
    from distutils import ccompiler

EXTRA_FLAGS_PER_COMPILER_TYPE_PER_PATH_COMPONENT = {
    "msvc": {
        "base64/arch/sse42": ["/arch:SSE4.2"],
        "base64/arch/avx2": ["/arch:AVX2"],
        "base64/arch/avx": ["/arch:AVX"],
    }
}

ccompiler.CCompiler.__spawn = ccompiler.CCompiler.spawn  # type: ignore[attr-defined]
X86_64 = platform.machine() in ("x86_64", "AMD64", "amd64")
PYODIDE = "PYODIDE" in os.environ
NO_EXTRA_FLAGS = "MYPYC_NO_EXTRA_FLAGS" in os.environ


def spawn(self, cmd, **kwargs) -> None:  # type: ignore[no-untyped-def]
    new_cmd = list(cmd)
    if PYODIDE:
        for argument in reversed(new_cmd):
            if not str(argument).endswith(".c"):
                continue
            if "base64/arch/" in str(argument):
                new_cmd.extend(["-msimd128"])
    elif not NO_EXTRA_FLAGS:
        compiler_type: str = self.compiler_type
        extra_options = EXTRA_FLAGS_PER_COMPILER_TYPE_PER_PATH_COMPONENT.get(compiler_type, None)
        if X86_64 and extra_options is not None:
            # filenames are closer to the end of command line
            for argument in reversed(new_cmd):
                # Check if the matching argument contains a source filename.
                if not str(argument).endswith(".c"):
                    continue

                for path in extra_options.keys():
                    if path in str(argument):
                        if compiler_type == "bcpp":
                            compiler = new_cmd.pop()
                            # Borland accepts a source file name at the end,
                            # insert the options before it
                            new_cmd.extend(extra_options[path])
                            new_cmd.append(compiler)
                        else:
                            new_cmd.extend(extra_options[path])

                        # path component is found, no need to search any further
                        break
    self.__spawn(new_cmd, **kwargs)


ccompiler.CCompiler.spawn = spawn  # type: ignore[method-assign]
