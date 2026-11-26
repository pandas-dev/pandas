import sys
import warnings

# pytz < 2023.3.post1 calls datetime.utcfromtimestamp at import, which Python 3.12
# deprecates. Remove once the pytz minimum version is at least 2023.3.post1
# (pinned lower in https://github.com/pandas-dev/pandas/pull/65133).
warnings.filterwarnings(
    "ignore",
    message=r"datetime\.datetime\.utcfromtimestamp\(\) is deprecated",
    category=DeprecationWarning,
    module=r"pytz\.tzinfo",
)

import pandas  # noqa: F401

blocklist = {
    "bs4",
    "gcsfs",
    "html5lib",
    "http",
    "ipython",
    "jinja2",
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
