from typing import Final

from gunicorn.glogging import Logger as GLogger

SD_LISTEN_FDS_START: Final[int]

def listen_fds(unset_environment: bool = True) -> int: ...
def sd_notify(state: str, logger: GLogger, unset_environment: bool = False) -> None: ...
