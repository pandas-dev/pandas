from _typeshed import Incomplete
from collections import defaultdict
from collections.abc import Iterable
from datetime import date, datetime
from decimal import Decimal

from pony.orm.core import Database, Entity

class Bag:
    database: Database
    session_cache: Incomplete
    entity_configs: dict[Entity, tuple[Incomplete, bool]]
    objects: defaultdict[type[Entity], set[Entity]]
    vars: dict[Incomplete, Incomplete]
    dicts: defaultdict[Incomplete, dict[Incomplete, Incomplete]]
    def __init__(bag, database: Database) -> None: ...
    def config(
        bag,
        entity: Entity,
        only=None,
        exclude=None,
        with_collections: bool = True,
        with_lazy: bool = False,
        related_objects: bool = True,
    ) -> tuple[Incomplete, bool]: ...
    def put(bag, x: Entity | Iterable[Entity]) -> None: ...
    def to_dict(bag): ...
    def to_json(bag) -> str: ...

def to_dict(objects): ...
def to_json(objects) -> str: ...
def json_converter(x: datetime | date | Decimal) -> str: ...
