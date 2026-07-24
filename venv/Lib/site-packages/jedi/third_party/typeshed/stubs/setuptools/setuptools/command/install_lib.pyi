from _typeshed import StrPath, Unused

from setuptools.dist import Distribution

from .._distutils.command import install_lib as orig

class install_lib(orig.install_lib):
    distribution: Distribution  # override distutils.dist.Distribution with setuptools.dist.Distribution
    def run(self) -> None: ...
    def get_exclusions(self): ...
    def copy_tree(
        self,
        infile: StrPath,
        outfile: str,
        preserve_mode: bool = True,
        preserve_times: bool = True,
        preserve_symlinks: bool = False,
        level: Unused = 1,
    ): ...
    def get_outputs(self): ...
