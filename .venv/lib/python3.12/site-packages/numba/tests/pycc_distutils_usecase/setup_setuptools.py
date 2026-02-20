from setuptools import setup

from source_module import cc


def run_setup():
    setup(ext_modules=[cc.distutils_extension()])


if __name__ == '__main__':
    run_setup()
