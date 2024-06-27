from .._distutils.command import bdist_rpm as orig

class bdist_rpm(orig.bdist_rpm):
    def run(self) -> None: ...
