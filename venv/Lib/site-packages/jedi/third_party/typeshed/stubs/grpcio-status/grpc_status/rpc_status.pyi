import grpc

from . import _async as aio

# Returns a google.rpc.status.Status message corresponding to a given grpc.Call.
def from_call(call: grpc.Call): ...

# Convert a google.rpc.status.Status message to grpc.Status.
def to_status(status) -> grpc.Status: ...

__all__ = ["from_call", "to_status", "aio"]
