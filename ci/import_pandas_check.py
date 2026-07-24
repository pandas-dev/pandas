import sys

import pandas  # noqa: F401

blocklist = {
    "bs4",
    "gcsfs",
    "html5lib",
    "http",
    "ipython",
    "jinja2",
    "hypothesis",
    "lxml",
    "matplotlib",
    "openpyxl",
    "py",
    "pytest",
    "s3fs",
    "scipy",
    "tables",
    "urllib.request",
    "xlrd",
    "xlsxwriter",
}

# GH#28227 for some of these check for top-level modules, while others are
#  more specific (e.g. urllib.request)
import_mods = {m.split(".")[0] for m in sys.modules} | set(sys.modules)
mods = blocklist & import_mods
if mods:
    raise Exception(f"pandas should not import {', '.join(mods)}")

# GH#41432 importing this module registers the pyarrow extension types for
# Period and Interval, so pyarrow recognizes them when reading pandas data
# directly rather than only after a parquet or feather round-trip.
if (
    "pyarrow" in sys.modules
    and "pandas.core.arrays.arrow.extension_types" not in sys.modules
):
    raise Exception("pandas should register its pyarrow extension types on import")
