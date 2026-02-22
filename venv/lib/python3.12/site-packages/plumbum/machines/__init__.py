from __future__ import annotations

from plumbum.machines.local import LocalCommand, LocalMachine, local
from plumbum.machines.remote import BaseRemoteMachine, RemoteCommand
from plumbum.machines.ssh_machine import PuttyMachine, SshMachine

__all__ = (
    "BaseRemoteMachine",
    "LocalCommand",
    "LocalMachine",
    "PuttyMachine",
    "RemoteCommand",
    "SshMachine",
    "local",
)
