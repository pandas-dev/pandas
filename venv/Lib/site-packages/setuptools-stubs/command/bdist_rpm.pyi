from setuptools.dist import Distribution

from .._distutils.command import bdist_rpm as orig

class bdist_rpm(orig.bdist_rpm):
    distribution: Distribution  # override distutils.dist.Distribution with setuptools.dist.Distribution
    def run(self) -> None: ...
