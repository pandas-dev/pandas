import configparser
from logging import Logger
from typing import Any

from .manager import OptionManager

LOG: Logger

def load_config(config: str | None, extra: list[str], *, isolated: bool = False) -> tuple[configparser.RawConfigParser, str]: ...
def parse_config(option_manager: OptionManager, cfg: configparser.RawConfigParser, cfg_dir: str) -> dict[str, Any]: ...
