"""Exceptions raised by the memorydb service."""

from typing import List

from moto.core.exceptions import JsonRESTError


class MemoryDBClientError(JsonRESTError):
    code = 400


class ClusterAlreadyExistsFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("ClusterAlreadyExistsFault", msg)


class InvalidSubnetError(MemoryDBClientError):
    def __init__(self, subnet_identifier: List[str]):
        super().__init__("InvalidSubnetError", f"Subnet {subnet_identifier} not found.")


class SubnetGroupAlreadyExistsFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("SubnetGroupAlreadyExistsFault", msg)


class ClusterNotFoundFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("ClusterNotFoundFault", msg)


class SnapshotAlreadyExistsFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("SnapshotAlreadyExistsFault", msg)


class SnapshotNotFoundFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("SnapshotNotFoundFault", msg)


class SubnetGroupNotFoundFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("SubnetGroupNotFoundFault", msg)


class TagNotFoundFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("TagNotFoundFault", msg)


class InvalidParameterValueException(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("InvalidParameterValueException", msg)


class SubnetGroupInUseFault(MemoryDBClientError):
    def __init__(self, msg: str):
        super().__init__("SubnetGroupInUseFault", msg)
