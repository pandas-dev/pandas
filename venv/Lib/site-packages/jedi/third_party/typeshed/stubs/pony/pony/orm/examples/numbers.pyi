from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Numbers(Entity):
    __slots__ = ()
    id: Incomplete
    int8: Incomplete
    int16: Incomplete
    int24: Incomplete
    int32: Incomplete
    int64: Incomplete
    uint8: Incomplete
    uint16: Incomplete
    uint24: Incomplete
    uint32: Incomplete

def populate_database() -> None: ...
def test_data() -> None: ...
