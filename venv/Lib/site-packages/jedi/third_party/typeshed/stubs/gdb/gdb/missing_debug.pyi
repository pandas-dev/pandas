import gdb
from gdb import Progspace
from gdb.missing_files import MissingFileHandler

class MissingDebugHandler(MissingFileHandler):
    def __call__(self, objfile: gdb.Objfile) -> bool | str | None: ...

def register_handler(locus: Progspace | None, handler: MissingDebugHandler, replace: bool = False) -> None: ...
