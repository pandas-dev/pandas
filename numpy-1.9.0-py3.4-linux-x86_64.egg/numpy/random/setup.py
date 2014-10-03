from __future__ import division, print_function

from os.path import join, split, dirname
import os
import sys
from distutils.dep_util import newer
from distutils.msvccompiler import get_build_version as get_msvc_build_version

def needs_mingw_ftime_workaround():
    # We need the mingw workaround for _ftime if the msvc runtime version is
    # 7.1 or above and we build with mingw ...
    # ... but we can't easily detect compiler version outside distutils command
    # context, so we will need to detect in randomkit whether we build with gcc
    msver = get_msvc_build_version()
    if msver and msver >= 8:
        return True

    return False

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration, get_mathlibs
    config = Configuration('random', parent_package, top_path)

    def generate_libraries(ext, build_dir):
        config_cmd = config.get_config_cmd()
        libs = get_mathlibs()
        tc = testcode_wincrypt()
        if config_cmd.try_run(tc):
            libs.append('Advapi32')
        ext.libraries.extend(libs)
        return None

    # enable unix large file support on 32 bit systems
    # (64 bit off_t, lseek -> lseek64 etc.)
    defs = [('_FILE_OFFSET_BITS', '64'),
            ('_LARGEFILE_SOURCE', '1'),
            ('_LARGEFILE64_SOURCE', '1')]
    if needs_mingw_ftime_workaround():
        defs.append(("NPY_NEEDS_MINGW_TIME_WORKAROUND", None))

    libs = []
    # Configure mtrand
    config.add_extension('mtrand',
                         sources=[join('mtrand', x) for x in
                                  ['mtrand.c', 'randomkit.c', 'initarray.c',
                                   'distributions.c']]+[generate_libraries],
                         libraries=libs,
                         depends=[join('mtrand', '*.h'),
                                  join('mtrand', '*.pyx'),
                                  join('mtrand', '*.pxi'),],
                         define_macros=defs,
                         )

    config.add_data_files(('.', join('mtrand', 'randomkit.h')))
    config.add_data_dir('tests')

    return config

def testcode_wincrypt():
    return """\
/* check to see if _WIN32 is defined */
int main(int argc, char *argv[])
{
#ifdef _WIN32
    return 0;
#else
    return 1;
#endif
}
"""

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
