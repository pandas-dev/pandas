import sys
from html import escape as escape
from io import FileIO, TextIOWrapper
from queue import Empty as Empty, Queue as Queue
from typing import IO

if sys.version_info >= (3, 13):
    # NOTE: These are the only attributes we realistically care about
    class cgi_FieldStorage:
        filename: str
        file: IO[bytes]
        def make_file(self) -> TextIOWrapper | FileIO: ...

    def parse_header(line: str) -> tuple[str, dict[str, str]]: ...

else:
    from cgi import FieldStorage as _cgi_FieldStorage, parse_header as parse_header

    class cgi_FieldStorage(_cgi_FieldStorage):
        # NOTE: The only kinds of objects of this type the user is exposed to
        #       will contain a file with a filename. We're technically lying
        #       if people create their own instances, but that shouldn't happen
        filename: str
        file: IO[bytes]
