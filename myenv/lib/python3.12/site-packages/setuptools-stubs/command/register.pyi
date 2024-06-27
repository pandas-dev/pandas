from .._distutils.command import register as orig

class register(orig.register):
    def run(self) -> None: ...
