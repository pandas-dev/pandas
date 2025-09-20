"""Well-known (and otherwise) constants used by JupyterLite"""

import shutil
from pathlib import Path

#: a default permission for directories
MOD_DIRECTORY = 0o755

#: a default permission for files
MOD_FILE = 0o644

#: a locale for reproducible file sorting
C_LOCALE = "C"

#: the encoding for pretty much every file written and read by jupyterlite
UTF8 = dict(encoding="utf-8")

#: default arguments for normalized JSON
JSON_FMT = dict(sort_keys=True, indent=2)

# the root of this project
ROOT = Path(__file__).parent

#: all of the archives
ALL_APP_ARCHIVES = sorted(ROOT.glob("jupyterlite-*.tgz"))

#: the extension point for addons, including core
ADDON_ENTRYPOINT = "jupyterlite.addon.v0"

### other parties' well-known paths
#: a predictably-serveable HTML file
INDEX_HTML = "index.html"

#: settings overrides. used JupyterLab build system, usually goes in
#: $PREFIX/share/jupyter/lab/
OVERRIDES_JSON = "overrides.json"

#: the canonical location within an env (or archive) for labextensions
SHARE_LABEXTENSIONS = "share/jupyter/labextensions"

#: the canonical location of labextension metadata
PACKAGE_JSON = "package.json"

#: the generally-used listing of pip requirements
REQUIREMENTS_TXT = "requirements.txt"

#: output equivalent to `sha256sum *` for providing a local bill-of-data
SHA256SUMS = "SHA256SUMS"

#: a script DOM ID on most jupyter pages
JUPYTER_CONFIG_DATA = "jupyter-config-data"

#: configuration key for prebuilt extensions
FEDERATED_EXTENSIONS = "federated_extensions"

#: configuration key for disabled extensions
DISABLED_EXTENSIONS = "disabledExtensions"

#: configuration key for extension settings overrides
SETTINGS_OVERRIDES = "settingsOverrides"

#: configuration key for file types
SETTINGS_FILE_TYPES = "fileTypes"

#: the top-level key for lite plugin settings
LITE_PLUGIN_SETTINGS = "litePluginSettings"

### jupyterlite "well-known" paths

#: our schema
JUPYTERLITE_SCHEMA = "jupyterlite.schema.v0.json"

#: our configuration file
JUPYTERLITE_JSON = "jupyter-lite.json"

#: our configuration file
JUPYTERLITE_IPYNB = "jupyter-lite.ipynb"
JUPYTERLITE_METADATA = "jupyter-lite"
JUPYTER_LITE_CONFIG = "jupyter_lite_config.json"

#: Needs a better canonical location
DEFAULT_OUTPUT_DIR = "_output"

#: commonly-used filenames for response fixtures, e.g. settings
ALL_JSON = "all.json"
ALL_FEDERATED_JSON = "all_federated.json"

### Environment Variables

#: a canonical environment variable for triggering reproducible builds
SOURCE_DATE_EPOCH = "SOURCE_DATE_EPOCH"

#: this is arrived at by inspection
NPM_SOURCE_DATE_EPOCH = 499162500

#: known zip extensions
EXTENSION_ZIP = (".whl", ".zip", ".conda")

#: known compressed tar extensions
EXTENSION_TAR = (".tgz", ".tar.bz2", ".tar.gz")

### URLs

#: the Jupyter API route for Contents API
API_CONTENTS = "api/contents"
API_TRANSLATIONS = "api/translations"
LAB_EXTENSIONS = "extensions"

#: our doit task-based plugin system
HOOKS = [
    "status",
    "init",
    "build",
    "check",
    "serve",
    "archive",
]

#: the name of the previous hook
HOOK_PARENTS = dict(
    build="init",
    check="build",
    serve="build",
    archive="build",
)

#: the lifecycle stages inside a hook
PHASES = ["pre_", "", "post_"]


#: extensions to be considered sourcemaps
SOURCEMAPS = [".js.map", ".mjs.map"]
SOURCEMAP_IGNORE_PATTERNS = shutil.ignore_patterns(*[f"*{p}" for p in SOURCEMAPS])

#: enough file types to serve all our demo files
DEFAULT_FILE_TYPES = dict(
    text=dict(
        css=[[".css"], ["text/css"]],
        csv=[[".csv"], ["text/csv"]],
        fasta=[[".fasta"], ["text/plain"]],
        html=[[".html"], ["text/html"]],
        ical=[[".ical", ".ics", ".ifb", ".icalendar"], ["text/calendar"]],
        js=[[".js", ".mjs"], ["application/javascript"]],
        manifest=[[".manifest"], ["text/cache-manifest"]],
        md=[[".md", ".markdown"], ["text/markdown"]],
        plain=[[".txt"], ["text/plain"]],
        py=[[".py"], ["text/x-python", "application/x-python-code"]],
        svg=[[".svg"], ["image/svg+xml"]],
        toml=[[".toml"], ["application/toml"]],
        vue=[[".vue"], ["text/plain"]],
        xml=[[".xml"], ["application/xml"]],
        yaml=[[".yaml", ".yml"], ["application/x-yaml"]],
    ),
    json=dict(
        geojson=[[".geojson"], ["application/geo+json"]],
        ipynb=[[".ipynb"], ["application/x-ipynb+json"]],
        jsmap=[[".map"], ["application/json"]],
        json=[[".json"], ["application/json"]],
    ),
    base64=dict(
        gzip=[[".tgz", ".gz", ".gzip"], ["application/gzip"]],
        ico=[[".ico"], ["image/x-icon"]],
        jpeg=[[".jpeg", ".jpg"], ["image/jpeg"]],
        pdf=[[".pdf"], ["application/pdf"]],
        png=[[".png"], ["image/png"]],
        wasm=[[".wasm"], ["application/wasm"]],
        wheel=[[".whl"], ["octet/stream", "application/x-wheel+zip"]],
    ),
)
