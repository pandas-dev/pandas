"""Well-known (and otherwise) constants used by ``jupyterlite-pyodide-kernel``"""

### Pyodide-specific values
#: the key for PyPI-compatible API responses pointing to wheels
PIPLITE_URLS = "pipliteUrls"
DISABLE_PYPI_FALLBACK = "disablePyPIFallback"
#: the schema for piplite-compatible wheel index
PIPLITE_INDEX_SCHEMA = "piplite.v0.schema.json"
#: the schema for the Pyodide kernel settings
KERNEL_SETTINGS_SCHEMA = "kernel.v0.schema.json"
#: where we put wheels, for now
PYPI_WHEELS = "pypi"
#: the plugin id for the Pyodide kernel labextension
PYODIDE_KERNEL_PLUGIN_ID = "@jupyterlite/pyodide-kernel-extension:kernel"
#: the npm name of the Pyodide kernel
PYODIDE_KERNEL_NPM_NAME = PYODIDE_KERNEL_PLUGIN_ID.split(":")[0]
#: the package.json key for piplite
PKG_JSON_PIPLITE = "piplite"
#: the package.json/piplite key for wheels
PKG_JSON_WHEELDIR = "wheelDir"

#: the jupyter-lite.json config key for the Pyodide base URL
PYODIDE_URL = "pyodideUrl"

#: directory name and filenames for the Pyodide distribution
PYODIDE = "pyodide"
#: the ES module entry point, required by Pyodide 314.0.0 and later
PYODIDE_MJS = "pyodide.mjs"
PYODIDE_LOCK_STEM = "pyodide-lock"
PYODIDE_LOCK = f"{PYODIDE_LOCK_STEM}.json"
PYODIDE_URL_ENV_VAR = "JUPYTERLITE_PYODIDE_URL"

#: probably only compatible with this version of pyodide
PYODIDE_VERSION = "314.0.1"

#: probably only compatible with this version of python in browser
PYODIDE_PYTHON_VERSION = "3.14"

#: the only kind of noarch wheel piplite understands
NOARCH_WHL = "py3-none-any.whl"

#: [pyodide <0.26] Emscripten platform tag (emscripten_*_wasm32)
EMSCRIPTEN_ABI_WHL = "emscripten_*_wasm32.whl"

#: variable alias for EMSCRIPTEN_ABI_WHL
WASM_WHL = EMSCRIPTEN_ABI_WHL

#: [pyodide <0.30] [micropip <0.11.1] Pyodide platform tag (pyodide_*_wasm32)
PYODIDE_ABI_WHL = "pyodide_*_wasm32.whl"

#: [pyodide >=0.30] [micropip >=0.11.1] PEP 783 platform tag (pyemscripten_*_wasm32)
#: See: https://peps.python.org/pep-0783/
#:      https://pyodide.org/en/stable/development/abi.html
PYEMSCRIPTEN_ABI_WHL = "pyemscripten_*_wasm32.whl"

#: all wheel filename patterns recognised by the piplite addon
ALL_WHL = [NOARCH_WHL, WASM_WHL, PYODIDE_ABI_WHL, PYEMSCRIPTEN_ABI_WHL]

#: a best-effort wheel filename pattern
RE_WHEEL_DIST_NAME = r"(?P<name>[a-zA-Z\d][a-z\d_\-\.]*[^\-])-[\d\.]+.*\.whl"

#: the default fallback URL prefix for pyodide packages
PYODIDE_CDN_URL = f"https://cdn.jsdelivr.net/pyodide/v{PYODIDE_VERSION}/full"

#: the default fallback URL for a lockfile
PYODIDE_LOCK_DEFAULT_URL = f"{PYODIDE_CDN_URL}/{PYODIDE_LOCK}"

#: the path to ``pyodide-lock``-downloaded wheels
PYODIDE_UV_WHEELS = "_uv_wheels"

#: configuration key for the loadPyodide options
LOAD_PYODIDE_OPTIONS = "loadPyodideOptions"

#: configuration key for the lockfile URL
OPTION_LOCK_FILE_URL = "lockFileURL"

#: configuration key for preloaded packages
OPTION_PACKAGES = "packages"

#: the python project file
PYPROJECT_TOML = "pyproject.toml"

#: key in a ``pyproject.toml`` for named dependency groups
PEP_735_DEP_GROUPS = "dependency-groups"

#: self-referential group in ``dependency-groups`` member
PEP_735_INC_GROUP = "include-group"
