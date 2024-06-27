from plumbum.machines.local import LocalCommand, LocalMachine, local
from plumbum.machines.remote import BaseRemoteMachine, RemoteCommand
from plumbum.machines.ssh_machine import PuttyMachine, SshMachine

__all__ = (
    "LocalCommand",
    "LocalMachine",
    "local",
    "BaseRemoteMachine",
    "RemoteCommand",
    "PuttyMachine",
    "SshMachine",
)
