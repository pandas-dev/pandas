from __future__ import division, absolute_import, print_function

from distutils.unixccompiler import UnixCCompiler
from numpy.distutils.exec_command import find_executable

class IntelCCompiler(UnixCCompiler):
    """ A modified Intel compiler compatible with an gcc built Python."""
    compiler_type = 'intel'
    cc_exe = 'icc'
    cc_args = 'fPIC'

    def __init__ (self, verbose=0, dry_run=0, force=0):
        UnixCCompiler.__init__ (self, verbose, dry_run, force)
        self.cc_exe = 'icc -fPIC'
        compiler = self.cc_exe
        self.set_executables(compiler=compiler,
                             compiler_so=compiler,
                             compiler_cxx=compiler,
                             linker_exe=compiler,
                             linker_so=compiler + ' -shared')

class IntelItaniumCCompiler(IntelCCompiler):
    compiler_type = 'intele'

    # On Itanium, the Intel Compiler used to be called ecc, let's search for
    # it (now it's also icc, so ecc is last in the search).
    for cc_exe in map(find_executable, ['icc', 'ecc']):
        if cc_exe:
            break

class IntelEM64TCCompiler(UnixCCompiler):
    """ A modified Intel x86_64 compiler compatible with a 64bit gcc built Python.
    """
    compiler_type = 'intelem'
    cc_exe = 'icc -m64 -fPIC'
    cc_args = "-fPIC"
    def __init__ (self, verbose=0, dry_run=0, force=0):
        UnixCCompiler.__init__ (self, verbose, dry_run, force)
        self.cc_exe = 'icc -m64 -fPIC'
        compiler = self.cc_exe
        self.set_executables(compiler=compiler,
                             compiler_so=compiler,
                             compiler_cxx=compiler,
                             linker_exe=compiler,
                             linker_so=compiler + ' -shared')
