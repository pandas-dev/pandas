"""
unixccompiler - can handle very long argument lists for ar.

"""
from __future__ import division, absolute_import, print_function

import os

from distutils.errors import DistutilsExecError, CompileError
from distutils.unixccompiler import *
from numpy.distutils.ccompiler import replace_method
from numpy.distutils.compat import get_exception

if sys.version_info[0] < 3:
    from . import log
else:
    from numpy.distutils import log

# Note that UnixCCompiler._compile appeared in Python 2.3
def UnixCCompiler__compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
    """Compile a single source files with a Unix-style compiler."""
    # HP ad-hoc fix, see ticket 1383
    ccomp = self.compiler_so
    if ccomp[0] == 'aCC':
        # remove flags that will trigger ANSI-C mode for aCC
        if '-Ae' in ccomp:
            ccomp.remove('-Ae')
        if '-Aa' in ccomp:
            ccomp.remove('-Aa')
        # add flags for (almost) sane C++ handling
        ccomp += ['-AA']
        self.compiler_so = ccomp

    display = '%s: %s' % (os.path.basename(self.compiler_so[0]), src)
    try:
        self.spawn(self.compiler_so + cc_args + [src, '-o', obj] +
                   extra_postargs, display = display)
    except DistutilsExecError:
        msg = str(get_exception())
        raise CompileError(msg)

replace_method(UnixCCompiler, '_compile', UnixCCompiler__compile)


def UnixCCompiler_create_static_lib(self, objects, output_libname,
                                    output_dir=None, debug=0, target_lang=None):
    """
    Build a static library in a separate sub-process.

    Parameters
    ----------
    objects : list or tuple of str
        List of paths to object files used to build the static library.
    output_libname : str
        The library name as an absolute or relative (if `output_dir` is used)
        path.
    output_dir : str, optional
        The path to the output directory. Default is None, in which case
        the ``output_dir`` attribute of the UnixCCompiler instance.
    debug : bool, optional
        This parameter is not used.
    target_lang : str, optional
        This parameter is not used.

    Returns
    -------
    None

    """
    objects, output_dir = self._fix_object_args(objects, output_dir)

    output_filename = \
                    self.library_filename(output_libname, output_dir=output_dir)

    if self._need_link(objects, output_filename):
        try:
            # previous .a may be screwed up; best to remove it first
            # and recreate.
            # Also, ar on OS X doesn't handle updating universal archives
            os.unlink(output_filename)
        except (IOError, OSError):
            pass
        self.mkpath(os.path.dirname(output_filename))
        tmp_objects = objects + self.objects
        while tmp_objects:
            objects = tmp_objects[:50]
            tmp_objects = tmp_objects[50:]
            display = '%s: adding %d object files to %s' % (
                           os.path.basename(self.archiver[0]),
                           len(objects), output_filename)
            self.spawn(self.archiver + [output_filename] + objects,
                       display = display)

        # Not many Unices required ranlib anymore -- SunOS 4.x is, I
        # think the only major Unix that does.  Maybe we need some
        # platform intelligence here to skip ranlib if it's not
        # needed -- or maybe Python's configure script took care of
        # it for us, hence the check for leading colon.
        if self.ranlib:
            display = '%s:@ %s' % (os.path.basename(self.ranlib[0]),
                                   output_filename)
            try:
                self.spawn(self.ranlib + [output_filename],
                           display = display)
            except DistutilsExecError:
                msg = str(get_exception())
                raise LibError(msg)
    else:
        log.debug("skipping %s (up-to-date)", output_filename)
    return

replace_method(UnixCCompiler, 'create_static_lib',
               UnixCCompiler_create_static_lib)
