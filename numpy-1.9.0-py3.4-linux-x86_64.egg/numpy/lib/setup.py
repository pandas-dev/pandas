from __future__ import division, print_function

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('lib', parent_package, top_path)

    config.add_include_dirs(join('..', 'core', 'include'))

    config.add_extension('_compiled_base',
                         sources=[join('src', '_compiled_base.c')]
                         )

    config.add_data_dir('benchmarks')
    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
