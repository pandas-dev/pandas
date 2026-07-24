from gdb import Progspace
from gdb.missing_files import MissingFileHandler

class MissingObjfileHandler(MissingFileHandler):
    def __call__(self, buildid: str, filename: str) -> bool | str | None: ...

def register_handler(locus: Progspace | None, handler: MissingObjfileHandler, replace: bool = False) -> None: ...
