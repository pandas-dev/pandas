# Integrated from the cffi-buildtool project by Rose Davidson
# (https://github.com/inklesspen/cffi-buildtool), under the following
# license:
#
# MIT License
#
# Copyright (c) 2024, Rose Davidson
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice (including the
# next paragraph) shall be included in all copies or substantial portions
# of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""Implementation of the ``cffi-gen-src`` command-line tool.

This module is private; the command line is the only supported
interface. Two subcommands:

``exec-python``
    Execute a Python script that constructs a :class:`cffi.FFI`
    (the same kind of script that the CFFI docs' "Main mode of usage"
    describes) and emit the generated C source.

``read-sources``
    Create the :class:`cffi.FFI` from a separate ``cdef`` file and C
    source prelude, then emit the generated C source.
"""

import argparse
import io
import os
import sys

from .api import FFI


def _execfile(pysrc, filename, globs):
    compiled = compile(source=pysrc, filename=filename, mode='exec')
    exec(compiled, globs, globs)


def find_ffi_in_python_script(pysrc, filename, ffivar):
    """Execute ``pysrc`` and return the :class:`FFI` object it defines.

    The script is executed with ``__name__`` set to ``"cffi.gen_src"``,
    so a trailing ``if __name__ == "__main__": ffibuilder.compile()``
    block in the script is skipped.

    ``ffivar`` is the name bound by the script to the :class:`FFI`
    object, or to a callable that returns one.

    Raises :class:`NameError` if the name is not bound by the script,
    or :class:`TypeError` if the name does not resolve to an
    :class:`FFI` instance.
    """
    filename = os.path.abspath(filename)
    globs = {'__name__': 'cffi.gen_src', '__file__': filename}
    old_path = sys.path[:]
    sys.path.insert(0, os.path.dirname(filename))
    try:
        _execfile(pysrc, filename, globs)
        if ffivar not in globs:
            raise NameError(
                "Expected to find the FFI object with the name %r, "
                "but it was not found." % (ffivar,)
            )
        ffi = globs[ffivar]
        if not isinstance(ffi, FFI) and callable(ffi):
            # Maybe it's a callable that returns a FFI
            ffi = ffi()
        if not isinstance(ffi, FFI):
            raise TypeError(
                "Found an object with the name %r but it was not an "
                "instance of cffi.api.FFI" % (ffivar,)
            )
        return ffi
    finally:
        sys.path[:] = old_path


def make_ffi_from_sources(modulename, cdef, csrc):
    """Create an :class:`FFI` from ``cdef`` text and a C source prelude."""
    ffibuilder = FFI()
    ffibuilder.cdef(cdef)
    ffibuilder.set_source(modulename, csrc)
    return ffibuilder


def generate_c_source(ffi):
    """Return the C source that :meth:`FFI.emit_c_code` would write."""
    output = io.StringIO()
    ffi.emit_c_code(output)
    return output.getvalue()


def write_c_source(output, generated):
    if output == '-':
        sys.stdout.write(generated)
        return
    with open(output, 'w', encoding='utf-8') as f:
        f.write(generated)


def same_input_file(file1, file2):
    if file1 is file2:
        return True
    try:
        return os.path.samefile(file1.name, file2.name)
    except OSError:
        return False


def exec_python(*, output, pyfile, ffi_var):
    with pyfile:
        ffi = find_ffi_in_python_script(pyfile.read(), pyfile.name, ffi_var)
    generated = generate_c_source(ffi)
    write_c_source(output, generated)


def read_sources(*, output, module_name, cdef_input, csrc_input):
    with csrc_input, cdef_input:
        csrc = csrc_input.read()
        cdef = cdef_input.read()
    ffi = make_ffi_from_sources(module_name, cdef, csrc)
    generated = generate_c_source(ffi)
    write_c_source(output, generated)


def _prog():
    # The same parser serves both documented invocations; make --help
    # and usage messages show the one that was actually used.
    argv0 = os.path.basename(sys.argv[0]) if sys.argv else ''
    if argv0.startswith('cffi-gen-src'):
        return 'cffi-gen-src'
    return 'python -m cffi.gen_src'


parser = argparse.ArgumentParser(
    prog=_prog(),
    description='Generate CFFI C source for extension modules.',
)
subparsers = parser.add_subparsers(dest='mode')

exec_python_parser = subparsers.add_parser(
    'exec-python',
    help='Execute a Python script that defines an FFI object',
)
exec_python_parser.add_argument(
    '--ffi-var',
    default='ffibuilder',
    help="Name of the FFI object in the Python script; defaults to 'ffibuilder'.",
)
exec_python_parser.add_argument(
    'pyfile',
    type=argparse.FileType('r', encoding='utf-8'),
    help='Path to the Python script',
)
exec_python_parser.add_argument(
    'output',
    help='Output path for the C source',
)

read_sources_parser = subparsers.add_parser(
    'read-sources',
    help='Read cdef and C source prelude files that define an FFI object',
)
read_sources_parser.add_argument(
    'module_name',
    help='Full name of the generated module, including packages',
)
read_sources_parser.add_argument(
    'cdef',
    type=argparse.FileType('r', encoding='utf-8'),
    help='File containing C definitions',
)
read_sources_parser.add_argument(
    'csrc',
    type=argparse.FileType('r', encoding='utf-8'),
    help='File containing C source prelude',
)
read_sources_parser.add_argument(
    'output',
    help='Output path for the C source',
)


def run(args=None):
    args = parser.parse_args(args=args)
    if args.mode == 'exec-python':
        exec_python(output=args.output, pyfile=args.pyfile, ffi_var=args.ffi_var)
    elif args.mode == 'read-sources':
        if same_input_file(args.cdef, args.csrc):
            parser.error('cdef and csrc are the same file and should not be')
        read_sources(
            output=args.output,
            module_name=args.module_name,
            cdef_input=args.cdef,
            csrc_input=args.csrc,
        )
    else:
        parser.error('a subcommand is required: exec-python or read-sources')
    parser.exit(0)
