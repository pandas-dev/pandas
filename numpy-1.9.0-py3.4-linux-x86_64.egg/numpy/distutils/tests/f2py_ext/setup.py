#!/usr/bin/env python
from __future__ import division, print_function

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('f2py_ext', parent_package, top_path)
    config.add_extension('fib2', ['src/fib2.pyf', 'src/fib1.f'])
    config.add_data_dir('tests')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
