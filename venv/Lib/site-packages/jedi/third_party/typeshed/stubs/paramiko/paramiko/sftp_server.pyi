from _typeshed import FileDescriptorOrPath
from logging import Logger
from typing import Any

from paramiko.channel import Channel
from paramiko.server import ServerInterface, SubsystemHandler
from paramiko.sftp import BaseSFTP
from paramiko.sftp_attr import SFTPAttributes
from paramiko.sftp_handle import SFTPHandle
from paramiko.sftp_si import SFTPServerInterface
from paramiko.transport import Transport

class SFTPServer(BaseSFTP, SubsystemHandler):
    logger: Logger
    ultra_debug: bool
    next_handle: int
    file_table: dict[bytes, SFTPHandle]
    folder_table: dict[bytes, SFTPHandle]
    server: SFTPServerInterface
    sock: Channel | None
    def __init__(
        self,
        channel: Channel,
        name: str,
        server: ServerInterface,
        sftp_si: type[SFTPServerInterface] = ...,
        *args: Any,
        **kwargs: Any,
    ) -> None: ...
    def start_subsystem(self, name: str, transport: Transport, channel: Channel) -> None: ...
    def finish_subsystem(self) -> None: ...
    @staticmethod
    def convert_errno(e: int) -> int: ...
    @staticmethod
    def set_file_attr(filename: FileDescriptorOrPath, attr: SFTPAttributes) -> None: ...
