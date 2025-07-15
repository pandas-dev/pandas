from setuptools import distutils
from source_module import cc


setup = distutils.core.setup


def run_setup():
    setup(ext_modules=[cc.distutils_extension()])


if __name__ == '__main__':
    run_setup()
