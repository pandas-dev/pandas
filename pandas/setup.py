#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pandas', parent_package, top_path)
    config.add_subpackage('core')
    config.add_subpackage('stats')
    config.add_subpackage('io')
    config.add_subpackage('util')
    config.add_subpackage('lib')
    return config

if __name__ == '__main__':
    print('This is the wrong setup.py file to run')

