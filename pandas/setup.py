#!/usr/bin/env python

import numpy

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pandas', parent_package, top_path)
    config.add_subpackage('core')
    config.add_subpackage('io')
    config.add_subpackage('rpy')
    config.add_subpackage('sandbox')
    config.add_subpackage('stats')
    config.add_subpackage('util')
    config.add_data_dir('tests')

    config.add_extension('_tseries',
                         sources=['src/tseries.c'],
                         include_dirs=[numpy.get_include()])
    config.add_extension('_sparse',
                         sources=['src/sparse.c'],
                         include_dirs=[numpy.get_include()])
    return config

if __name__ == '__main__':
    print('This is the wrong setup.py file to run')

