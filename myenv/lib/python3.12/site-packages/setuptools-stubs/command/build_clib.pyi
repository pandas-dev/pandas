from .._distutils.command import build_clib as orig

class build_clib(orig.build_clib):
    def build_libraries(self, libraries) -> None: ...
