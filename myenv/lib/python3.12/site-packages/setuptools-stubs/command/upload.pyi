from .._distutils.command import upload as orig

class upload(orig.upload):
    def run(self) -> None: ...
