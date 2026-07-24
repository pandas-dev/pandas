import argparse
from logging import Logger

from ..formatting.base import BaseFormatter
from .finder import LoadedPlugin

LOG: Logger

def make(reporters: dict[str, LoadedPlugin], options: argparse.Namespace) -> BaseFormatter: ...
