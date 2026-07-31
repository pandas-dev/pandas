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
