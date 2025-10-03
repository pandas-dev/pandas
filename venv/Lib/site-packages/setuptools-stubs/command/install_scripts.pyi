from setuptools.dist import Distribution

from .._distutils.command import install_scripts as orig

class install_scripts(orig.install_scripts):
    distribution: Distribution  # override distutils.dist.Distribution with setuptools.dist.Distribution
    no_ep: bool
    def initialize_options(self) -> None: ...
    outfiles: list[str]
    def run(self) -> None: ...
    def write_script(self, script_name, contents, mode: str = "t", *ignored) -> None: ...
