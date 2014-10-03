#!/usr/bin/env python
from __future__ import division, print_function

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('testnumpydistutils', parent_package, top_path)
    config.add_subpackage('pyrex_ext')
    config.add_subpackage('f2py_ext')
    #config.add_subpackage('f2py_f90_ext')
    config.add_subpackage('swig_ext')
    config.add_subpackage('gen_ext')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
