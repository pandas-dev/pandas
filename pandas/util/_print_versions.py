import codecs
import locale
import os
import platform
import struct
import subprocess
import sys

from pandas.compat._optional import VERSIONS, _get_version, import_optional_dependency


def get_sys_info():
    "Returns system information as a dict"

    blob = []

    # get full commit hash
    commit = None
    if os.path.isdir(".git") and os.path.isdir("pandas"):
        try:
            pipe = subprocess.Popen(
                'git log --format="%H" -n 1'.split(" "),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            so, serr = pipe.communicate()
        except (OSError, ValueError):
            pass
        else:
            if pipe.returncode == 0:
                commit = so
                try:
                    commit = so.decode("utf-8")
                except ValueError:
                    pass
                commit = commit.strip().strip('"')

    blob.append(("commit", commit))

    try:
        (sysname, nodename, release, version, machine, processor) = platform.uname()
        blob.extend(
            [
                ("python", ".".join(map(str, sys.version_info))),
                ("python-bits", struct.calcsize("P") * 8),
                ("OS", "{sysname}".format(sysname=sysname)),
                ("OS-release", "{release}".format(release=release)),
                # ("Version", "{version}".format(version=version)),
                ("machine", "{machine}".format(machine=machine)),
                ("processor", "{processor}".format(processor=processor)),
                ("byteorder", "{byteorder}".format(byteorder=sys.byteorder)),
                ("LC_ALL", "{lc}".format(lc=os.environ.get("LC_ALL", "None"))),
                ("LANG", "{lang}".format(lang=os.environ.get("LANG", "None"))),
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
        if mod:
            ver = _get_version(mod)
        else:
            ver = None
        deps_blob.append((modname, ver))

    if as_json:
        try:
            import json
        except ImportError:
            import simplejson as json

        j = dict(system=dict(sys_info), dependencies=dict(deps_blob))

        if as_json is True:
            print(j)
        else:
            with codecs.open(as_json, "wb", encoding="utf8") as f:
                json.dump(j, f, indent=2)

    else:
        maxlen = max(len(x) for x in deps)
        tpl = "{{k:<{maxlen}}}: {{stat}}".format(maxlen=maxlen)
        print("\nINSTALLED VERSIONS")
        print("------------------")
        for k, stat in sys_info:
            print(tpl.format(k=k, stat=stat))
        print("")
        for k, stat in deps_blob:
            print(tpl.format(k=k, stat=stat))


def main():
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option(
        "-j",
        "--json",
        metavar="FILE",
        nargs=1,
        help="Save output as JSON into file, pass in " "'-' to output to stdout",
    )

    (options, args) = parser.parse_args()

    if options.json == "-":
        options.json = True

    show_versions(as_json=options.json)

    return 0


if __name__ == "__main__":
    sys.exit(main())
