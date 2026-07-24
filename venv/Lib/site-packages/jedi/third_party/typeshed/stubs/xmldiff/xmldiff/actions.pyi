from _typeshed import Incomplete
from typing import NamedTuple

class DeleteNode(NamedTuple):
    node: str

class InsertNode(NamedTuple):
    target: Incomplete
    tag: str
    position: int

class RenameNode(NamedTuple):
    node: str
    tag: Incomplete

class MoveNode(NamedTuple):
    node: str
    target: Incomplete
    position: int

class UpdateTextIn(NamedTuple):
    node: str
    text: Incomplete

class UpdateTextAfter(NamedTuple):
    node: str
    text: Incomplete

class UpdateAttrib(NamedTuple):
    node: str
    name: str
    value: Incomplete

class DeleteAttrib(NamedTuple):
    node: str
    name: str

class InsertAttrib(NamedTuple):
    node: str
    name: str
    value: Incomplete

class RenameAttrib(NamedTuple):
    node: str
    oldname: str
    newname: str

class InsertComment(NamedTuple):
    target: Incomplete
    position: Incomplete
    text: Incomplete

class InsertNamespace(NamedTuple):
    prefix: str
    uri: str

class DeleteNamespace(NamedTuple):
    prefix: str
