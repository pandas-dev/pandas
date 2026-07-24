from typing import Final

from .html5parser import HTMLParser as HTMLParser, parse as parse, parseFragment as parseFragment
from .serializer import serialize as serialize
from .treebuilders import getTreeBuilder as getTreeBuilder
from .treewalkers import getTreeWalker as getTreeWalker

__all__ = ["HTMLParser", "parse", "parseFragment", "getTreeBuilder", "getTreeWalker", "serialize"]

__version__: Final[str]
