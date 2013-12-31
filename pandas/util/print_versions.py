import os
import platform
import sys
import struct

def get_sys_info():
    "Returns system information as a dict"

    # list of tuples over dict because OrderedDict not in 2.6, least
    # resistance.
    blob = []
    try:
        sysname, nodename, release, version, machine, processor = platform.uname()
        blob = [
            ("python", "%d.%d.%d.%s.%s" % sys.version_info[:]),
            ("python-bits", struct.calcsize("P") * 8),
            ("OS","%s" % (sysname)),
            ("OS-release", "%s" % (release)),
            # ("Version", "%s" % (version)),
            # ("Machine", "%s" % (machine)),
            ("processor", "%s" % (processor)),
            ("byteorder", "%s" % sys.byteorder),
            ("LC_ALL", "%s" % os.environ.get('LC_ALL', "None")),
            ("LANG", "%s" % os.environ.get('LANG', "None")),

        ]
    except:
        pass

    return blob


def show_versions(as_json=False):
    import imp
    sys_info = get_sys_info()

    deps = [
        # (MODULE_NAME, f(mod) -> mod version)
        ("pandas", lambda mod: mod.__version__),
        ("Cython", lambda mod: mod.__version__),
        ("numpy", lambda mod: mod.version.version),
        ("scipy", lambda mod: mod.version.version),
        ("statsmodels", lambda mod: mod.__version__),
        ("patsy", lambda mod: mod.__version__),
        ("scikits.timeseries", lambda mod: mod.__version__),
        ("dateutil", lambda mod: mod.__version__),
        ("pytz", lambda mod: mod.VERSION),
        ("bottleneck", lambda mod: mod.__version__),
        ("tables", lambda mod: mod.__version__),
        ("numexpr", lambda mod: mod.__version__),
        ("matplotlib", lambda mod: mod.__version__),
        ("openpyxl", lambda mod: mod.__version__),
        ("xlrd", lambda mod: mod.__VERSION__),
        ("xlwt", lambda mod: mod.__VERSION__),
        ("xlsxwriter", lambda mod: mod.__version__),
        ("sqlalchemy", lambda mod: mod.__version__),
        ("lxml", lambda mod: mod.etree.__version__),
        ("bs4", lambda mod: mod.__version__),
        ("html5lib", lambda mod: mod.__version__),
        ("bq", lambda mod: mod._VersionNumber()),
        ("apiclient", lambda mod: mod.__version__),
    ]

    deps_blob = list()
    for (modname, ver_f) in deps:
        try:
            mod = imp.load_module(modname, *imp.find_module(modname))
            ver = ver_f(mod)
            deps_blob.append((modname, ver))
        except:
            deps_blob.append((modname, None))

    if (as_json):
        # 2.6-safe
        try:
            import json
        except:
            import simplejson as json

        print(json.dumps(dict(system=dict(sys_info), dependencies=dict(deps_blob)), indent=2))

    else:

        print("\nINSTALLED VERSIONS")
        print("------------------")

        for k, stat in sys_info:
            print("%s: %s" % (k, stat))

        print("")
        for k, stat in deps_blob:
            print("%s: %s" % (k, stat))


if __name__ == "__main__":
    show_versions(as_json=False)
