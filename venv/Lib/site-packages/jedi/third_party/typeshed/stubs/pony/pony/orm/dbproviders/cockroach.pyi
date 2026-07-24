from _typeshed import Incomplete
from typing import ClassVar

from pony.orm import dbapiprovider
from pony.orm.dbproviders.postgres import (
    PGArrayConverter,
    PGBlobConverter,
    PGColumn,
    PGIntConverter,
    PGProvider,
    PGSchema,
    PGSQLBuilder,
    PGTimedeltaConverter,
    PGTranslator,
)

NoneType: type[None]

class CRColumn(PGColumn):
    auto_template: ClassVar[str]

class CRSchema(PGSchema):
    column_class: ClassVar[type[CRColumn]]

class CRTranslator(PGTranslator): ...
class CRSQLBuilder(PGSQLBuilder): ...

class CRIntConverter(PGIntConverter):
    signed_types: Incomplete
    unsigned_types: Incomplete

class CRBlobConverter(PGBlobConverter):
    def sql_type(converter): ...

class CRTimedeltaConverter(PGTimedeltaConverter): ...

class PGUuidConverter(dbapiprovider.UuidConverter):
    def py2sql(converter, val): ...

class CRArrayConverter(PGArrayConverter):
    array_types: Incomplete

class CRProvider(PGProvider):
    dbschema_cls: ClassVar[type[CRSchema]]
    translator_cls: ClassVar[type[CRTranslator]]
    sqlbuilder_cls: ClassVar[type[CRSQLBuilder]]
    array_converter_cls: ClassVar[type[CRArrayConverter]]

provider_cls = CRProvider
