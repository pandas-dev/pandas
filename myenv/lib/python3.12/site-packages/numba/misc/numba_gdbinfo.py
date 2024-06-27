"""Module for displaying information about Numba's gdb set up"""
from collections import namedtuple
import os
import re
import subprocess
from textwrap import dedent
from numba import config

# Container for the output of the gdb info data collection
_fields = ('binary_loc, extension_loc, py_ver, np_ver, supported')
_gdb_info = namedtuple('_gdb_info', _fields)


class _GDBTestWrapper():
    """Wraps the gdb binary and has methods for checking what the gdb binary
    has support for (Python and NumPy)."""

    def __init__(self,):
        gdb_binary = config.GDB_BINARY
        if gdb_binary is None:
            msg = ("No valid binary could be found for gdb named: "
                   f"{config.GDB_BINARY}")
            raise ValueError(msg)
        self._gdb_binary = gdb_binary

    def _run_cmd(self, cmd=()):
        gdb_call = [self.gdb_binary, '-q',]
        for x in cmd:
            gdb_call.append('-ex')
            gdb_call.append(x)
        gdb_call.extend(['-ex', 'q'])
        return subprocess.run(gdb_call, capture_output=True, timeout=10,
                              text=True)

    @property
    def gdb_binary(self):
        return self._gdb_binary

    @classmethod
    def success(cls, status):
        return status.returncode == 0

    def check_launch(self):
        """Checks that gdb will launch ok"""
        return self._run_cmd()

    def check_python(self):
        cmd = ("python from __future__ import print_function; "
               "import sys; print(sys.version_info[:2])")
        return self._run_cmd((cmd,))

    def check_numpy(self):
        cmd = ("python from __future__ import print_function; "
               "import types; import numpy; "
               "print(isinstance(numpy, types.ModuleType))")
        return self._run_cmd((cmd,))

    def check_numpy_version(self):
        cmd = ("python from __future__ import print_function; "
               "import types; import numpy;"
               "print(numpy.__version__)")
        return self._run_cmd((cmd,))


def collect_gdbinfo():
    """Prints information to stdout about the gdb setup that Numba has found"""

    # State flags:
    gdb_state = None
    gdb_has_python = False
    gdb_has_numpy = False
    gdb_python_version = 'No Python support'
    gdb_python_numpy_version = "No NumPy support"

    # There are so many ways for gdb to not be working as expected. Surround
    # the "is it working" tests with try/except and if there's an exception
    # store it for processing later.
    try:
        # Check gdb exists
        gdb_wrapper = _GDBTestWrapper()

        # Check gdb works
        status = gdb_wrapper.check_launch()
        if not gdb_wrapper.success(status):
            msg = (f"gdb at '{gdb_wrapper.gdb_binary}' does not appear to work."
                   f"\nstdout: {status.stdout}\nstderr: {status.stderr}")
            raise ValueError(msg)
        gdb_state = gdb_wrapper.gdb_binary
    except Exception as e:
        gdb_state = f"Testing gdb binary failed. Reported Error: {e}"
    else:
        # Got this far, so gdb works, start checking what it supports
        status = gdb_wrapper.check_python()
        if gdb_wrapper.success(status):
            version_match = re.match(r'\((\d+),\s+(\d+)\)',
                                     status.stdout.strip())
            if version_match is not None:
                pymajor, pyminor = version_match.groups()
                gdb_python_version = f"{pymajor}.{pyminor}"
                gdb_has_python = True

                status = gdb_wrapper.check_numpy()
                if gdb_wrapper.success(status):
                    if "Traceback" not in status.stderr.strip():
                        if status.stdout.strip() == 'True':
                            gdb_has_numpy = True
                            gdb_python_numpy_version = "Unknown"
                            # NumPy is present find the version
                            status = gdb_wrapper.check_numpy_version()
                            if gdb_wrapper.success(status):
                                if "Traceback" not in status.stderr.strip():
                                    gdb_python_numpy_version = \
                                        status.stdout.strip()

    # Work out what level of print-extension support is present in this gdb
    if gdb_has_python:
        if gdb_has_numpy:
            print_ext_supported = "Full (Python and NumPy supported)"
        else:
            print_ext_supported = "Partial (Python only, no NumPy support)"
    else:
        print_ext_supported = "None"

    # Work out print ext location
    print_ext_file = "gdb_print_extension.py"
    print_ext_path = os.path.join(os.path.dirname(__file__), print_ext_file)

    # return!
    return _gdb_info(gdb_state, print_ext_path, gdb_python_version,
                     gdb_python_numpy_version, print_ext_supported)


def display_gdbinfo(sep_pos=45):
    """Displays the information collected by collect_gdbinfo.
    """
    gdb_info = collect_gdbinfo()
    print('-' * 80)
    fmt = f'%-{sep_pos}s : %-s'
    # Display the information
    print(fmt % ("Binary location", gdb_info.binary_loc))
    print(fmt % ("Print extension location", gdb_info.extension_loc))
    print(fmt % ("Python version", gdb_info.py_ver))
    print(fmt % ("NumPy version", gdb_info.np_ver))
    print(fmt % ("Numba printing extension support", gdb_info.supported))

    print("")
    print("To load the Numba gdb printing extension, execute the following "
          "from the gdb prompt:")
    print(f"\nsource {gdb_info.extension_loc}\n")
    print('-' * 80)
    warn = """
    =============================================================
    IMPORTANT: Before sharing you should remove any information
    in the above that you wish to keep private e.g. paths.
    =============================================================
    """
    print(dedent(warn))


if __name__ == '__main__':
    display_gdbinfo()
