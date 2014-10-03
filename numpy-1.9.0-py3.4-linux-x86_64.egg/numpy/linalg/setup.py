from __future__ import division, print_function

import os
import sys

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('linalg', parent_package, top_path)

    config.add_data_dir('tests')

    # Configure lapack_lite

    src_dir = 'lapack_lite'
    lapack_lite_src = [
        os.path.join(src_dir, 'python_xerbla.c'),
        os.path.join(src_dir, 'zlapack_lite.c'),
        os.path.join(src_dir, 'dlapack_lite.c'),
        os.path.join(src_dir, 'blas_lite.c'),
        os.path.join(src_dir, 'dlamch.c'),
        os.path.join(src_dir, 'f2c_lite.c'),
        os.path.join(src_dir, 'f2c.h'),
    ]

    lapack_info = get_info('lapack_opt', 0) # and {}
    def get_lapack_lite_sources(ext, build_dir):
        if not lapack_info:
            print("### Warning:  Using unoptimized lapack ###")
            return ext.depends[:-1]
        else:
            if sys.platform=='win32':
                print("### Warning:  python_xerbla.c is disabled ###")
                return ext.depends[:1]
            return ext.depends[:2]

    config.add_extension('lapack_lite',
                         sources = [get_lapack_lite_sources],
                         depends = ['lapack_litemodule.c'] + lapack_lite_src,
                         extra_info = lapack_info
                         )

    # umath_linalg module

    config.add_extension('_umath_linalg',
                         sources = [get_lapack_lite_sources],
                         depends =  ['umath_linalg.c.src'] + lapack_lite_src,
                         extra_info = lapack_info,
                         libraries = ['npymath'],
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
