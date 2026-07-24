"""Utility functions for printing version information."""

import contextlib
import importlib
import locale
import os
import platform
import struct
import subprocess
import sys


def get_sys_info():
    """Returns system information as a dict"""

    blob = []

    # get full commit hash
    commit = None
    if os.path.isdir(".git") and os.path.isdir("xarray"):
        try:
            pipe = subprocess.Popen(
                ("git", "log", '--format="%H"', "-n", "1"),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            so, _ = pipe.communicate()
        except Exception:
            pass
        else:
            if pipe.returncode == 0:
                commit = so
                with contextlib.suppress(ValueError):
                    commit = so.decode("utf-8")
                commit = commit.strip().strip('"')

    blob.append(("commit", commit))

    try:
        (sysname, _nodename, release, _version, machine, processor) = platform.uname()
        blob.extend(
            [
                ("python", sys.version),
                ("python-bits", struct.calcsize("P") * 8),
                ("OS", f"{sysname}"),
                ("OS-release", f"{release}"),
                # ("Version", f"{version}"),
                ("machine", f"{machine}"),
                ("processor", f"{processor}"),
                ("byteorder", f"{sys.byteorder}"),
                ("LC_ALL", f"{os.environ.get('LC_ALL', 'None')}"),
                ("LANG", f"{os.environ.get('LANG', 'None')}"),
                ("LOCALE", f"{locale.getlocale()}"),
            ]
        )
    except Exception:
        pass

    return blob


def netcdf_and_hdf5_versions():
    libhdf5_version = None
    libnetcdf_version = None
    try:
        import netCDF4

        libhdf5_version = netCDF4.__hdf5libversion__
        libnetcdf_version = netCDF4.__netcdf4libversion__
    except ImportError:
        try:
            import h5py

            libhdf5_version = h5py.version.hdf5_version
        except ImportError:
            pass
    return [("libhdf5", libhdf5_version), ("libnetcdf", libnetcdf_version)]


def show_versions(file=sys.stdout):
    """print the versions of xarray and its dependencies

    Parameters
    ----------
    file : file-like, optional
        print to the given file-like object. Defaults to sys.stdout.
    """
    sys_info = get_sys_info()

    try:
        sys_info.extend(netcdf_and_hdf5_versions())
    except Exception as e:
        print(f"Error collecting netcdf / hdf5 version: {e}")

    deps = [
        # (MODULE_NAME, f(mod) -> mod version)
        ("xarray", lambda mod: mod.__version__),
        ("pandas", lambda mod: mod.__version__),
        ("numpy", lambda mod: mod.__version__),
        ("scipy", lambda mod: mod.__version__),
        # xarray optionals
        ("netCDF4", lambda mod: mod.__version__),
        ("pydap", lambda mod: mod.__version__),
        ("h5netcdf", lambda mod: mod.__version__),
        ("h5py", lambda mod: mod.__version__),
        ("zarr", lambda mod: mod.__version__),
        ("cftime", lambda mod: mod.__version__),
        ("nc_time_axis", lambda mod: mod.__version__),
        ("iris", lambda mod: mod.__version__),
        ("bottleneck", lambda mod: mod.__version__),
        ("dask", lambda mod: mod.__version__),
        ("distributed", lambda mod: mod.__version__),
        ("matplotlib", lambda mod: mod.__version__),
        ("cartopy", lambda mod: mod.__version__),
        ("seaborn", lambda mod: mod.__version__),
        ("numbagg", lambda mod: mod.__version__),
        ("fsspec", lambda mod: mod.__version__),
        ("cupy", lambda mod: mod.__version__),
        ("pint", lambda mod: mod.__version__),
        ("sparse", lambda mod: mod.__version__),
        ("flox", lambda mod: mod.__version__),
        ("numpy_groupies", lambda mod: mod.__version__),
        # xarray setup/test
        ("setuptools", lambda mod: mod.__version__),
        ("pip", lambda mod: mod.__version__),
        ("conda", lambda mod: mod.__version__),
        ("pytest", lambda mod: mod.__version__),
        ("mypy", lambda mod: importlib.metadata.version(mod.__name__)),
        # Misc.
        ("IPython", lambda mod: mod.__version__),
        ("sphinx", lambda mod: mod.__version__),
    ]

    deps_blob = []
    for modname, ver_f in deps:
        try:
            if modname in sys.modules:
                mod = sys.modules[modname]
            else:
                mod = importlib.import_module(modname)
        except Exception:
            deps_blob.append((modname, None))
        else:
            try:
                ver = ver_f(mod)
                deps_blob.append((modname, ver))
            except Exception:
                deps_blob.append((modname, "installed"))

    print("\nINSTALLED VERSIONS", file=file)
    print("------------------", file=file)

    for k, stat in sys_info:
        print(f"{k}: {stat}", file=file)

    print("", file=file)
    for k, stat in deps_blob:
        print(f"{k}: {stat}", file=file)


if __name__ == "__main__":
    show_versions()
