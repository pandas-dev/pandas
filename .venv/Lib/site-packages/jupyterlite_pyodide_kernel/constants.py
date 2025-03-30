"""Well-known (and otherwise) constants used by ``jupyterlite-pyodide-kernel``"""

### Pyodide-specific values
#: the key for PyPI-compatible API responses pointing to wheels
PIPLITE_URLS = "pipliteUrls"
DISABLE_PYPI_FALLBACK = "disablePyPIFallback"
#: the schema for piplite-compatible wheel index
PIPLITE_INDEX_SCHEMA = "piplite.v0.schema.json"
#: the schema for piplite-compatible wheel index
KERNEL_SETTINGS_SCHEMA = "kernel.v0.schema.json"
#: where we put wheels, for now
PYPI_WHEELS = "pypi"
#: the plugin id for the Pydodide kernel labextension
PYODIDE_KERNEL_PLUGIN_ID = "@jupyterlite/pyodide-kernel-extension:kernel"
#: the npm name of the Pyodide kernel
PYODIDE_KERNEL_NPM_NAME = PYODIDE_KERNEL_PLUGIN_ID.split(":")[0]
#: the package.json key for piplite
PKG_JSON_PIPLITE = "piplite"
#: the package.json/piplite key for wheels
PKG_JSON_WHEELDIR = "wheelDir"

#: where we put wheels, for now
PYODIDE_URL = "pyodideUrl"

#: where we put pyodide, for now
PYODIDE = "pyodide"
PYODIDE_JS = "pyodide.js"
PYODIDE_LOCK = "pyodide-lock.json"
PYODIDE_URL_ENV_VAR = "JUPYTERLITE_PYODIDE_URL"

#: probably only compatible with this version of pyodide
PYODIDE_VERSION = "0.27.1"

#: the only kind of noarch wheel piplite understands
NOARCH_WHL = "py3-none-any.whl"

#: the only kind of binary wheel piplite previously understood
EMSCRIPTEN_ABI_WHL = "emscripten_*_wasm32.whl"

#: legacy variable alias
WASM_WHL = EMSCRIPTEN_ABI_WHL

#: the Pyodide ABI wheel is the same as the Emscripten
#: ABI wheel, but with a different platform tag, i.e.,
#  YYYY_buildnumber.
PYODIDE_ABI_WHL = "pyodide_*_wasm32.whl"

ALL_WHL = [NOARCH_WHL, WASM_WHL, PYODIDE_ABI_WHL]
