from gunicorn.arbiter import Arbiter
from gunicorn.ctl.handlers import CommandHandlers

class ControlSocketServer:
    arbiter: Arbiter
    socket_path: str
    socket_mode: int
    handlers: CommandHandlers
    def __init__(self, arbiter: Arbiter, socket_path: str, socket_mode: int = 0o600) -> None: ...
    def start(self) -> None: ...
    def stop(self) -> None: ...
