import argparse
import configparser
from collections.abc import Sequence
from logging import Logger

from .manager import OptionManager

LOG: Logger

def aggregate_options(
    manager: OptionManager, cfg: configparser.RawConfigParser, cfg_dir: str, argv: Sequence[str] | None
) -> argparse.Namespace: ...
