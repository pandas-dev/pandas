import codecs
import json
import locale
import os
import platform
import struct
import sys
from typing import List, Optional, Tuple, Union

from pandas.compat._optional import VERSIONS, _get_version, import_optional_dependency


def _get_commit_hash() -> Optional[str]:
    """
    Use vendored versioneer code to get git hash, which handles
    git worktree correctly.
    """
    from pandas._version import get_versions

    versions = get_versions()
    return versions["full-revisionid"]


def get_sys_info() -> List[Tuple[str, Optional[Union[str, int]]]]:
    """
    Returns system information as a list
    """
    blob: List[Tuple[str, Optional[Union[str, int]]]] = []

    # get full commit hash
    commit = _get_commit_hash()

    blob.append(("commit", commit))

    try:
        (sysname, nodename, release, version, machine, processor) = platform.uname()
        blob.extend(
            [
                ("python", ".".join(map(str, sys.version_info))),
                ("python-bits", struct.calcsize("P") * 8),
                ("OS", f"{sysname}"),
                ("OS-release", f"{release}"),
                # FIXME: dont leave commented-out
                # ("Version", f"{version}"),
                ("machine", f"{machine}"),
                ("processor", f"{processor}"),
                ("byteorder", f"{sys.byteorder}"),
                ("LC_ALL", f"{os.environ.get('LC_ALL', 'None')}"),
                ("LANG", f"{os.environ.get('LANG', 'None')}"),
                ("LOCALE", ".".join(map(str, locale.getlocale()))),
            ]
        )
    except (KeyError, ValueError):
        pass

    return blob


def show_versions(as_json=False):
    sys_info = get_sys_info()
    deps = [
        "pandas",
        # required
        "numpy",
        "pytz",
        "dateutil",
        # install / build,
        "pip",
        "setuptools",
        "Cython",
        # test
        "pytest",
        "hypothesis",
        # docs
        "sphinx",
        # Other, need a min version
        "blosc",
        "feather",
        "xlsxwriter",
        "lxml.etree",
        "html5lib",
        "pymysql",
        "psycopg2",
        "jinja2",
        # Other, not imported.
        "IPython",
        "pandas_datareader",
    ]

    deps.extend(list(VERSIONS))
    deps_blob = []

    for modname in deps:
        mod = import_optional_dependency(
            modname, raise_on_missing=False, on_version="ignore"
        )
        ver: Optional[str]
        if mod:
            ver = _get_version(mod)
        else:
            ver = None
        deps_blob.append((modname, ver))

    if as_json:
        j = dict(system=dict(sys_info), dependencies=dict(deps_blob))

        if as_json is True:
            print(j)
        else:
            with codecs.open(as_json, "wb", encoding="utf8") as f:
                json.dump(j, f, indent=2)

    else:
        maxlen = max(len(x) for x in deps)
        print("\nINSTALLED VERSIONS")
        print("------------------")
        for k, stat in sys_info:
            print(f"{k:<{maxlen}}: {stat}")
        print("")
        for k, stat in deps_blob:
            print(f"{k:<{maxlen}}: {stat}")


def main() -> int:
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option(
        "-j",
        "--json",
        metavar="FILE",
        nargs=1,
        help="Save output as JSON into file, pass in '-' to output to stdout",
    )

    (options, args) = parser.parse_args()

    if options.json == "-":
        options.json = True

    show_versions(as_json=options.json)

    return 0


if __name__ == "__main__":
    sys.exit(main())
