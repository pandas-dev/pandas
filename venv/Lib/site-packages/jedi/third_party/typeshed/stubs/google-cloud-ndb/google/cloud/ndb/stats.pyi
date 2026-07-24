from _typeshed import Incomplete

from google.cloud.ndb import model

class BaseStatistic(model.Model):
    STORED_KIND_NAME: str
    bytes: Incomplete
    count: Incomplete
    timestamp: Incomplete

class BaseKindStatistic(BaseStatistic):
    STORED_KIND_NAME: str
    kind_name: Incomplete
    entity_bytes: Incomplete

class GlobalStat(BaseStatistic):
    STORED_KIND_NAME: str
    entity_bytes: Incomplete
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete
    composite_index_bytes: Incomplete
    composite_index_count: Incomplete

class NamespaceStat(BaseStatistic):
    STORED_KIND_NAME: str
    subject_namespace: Incomplete
    entity_bytes: Incomplete
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete
    composite_index_bytes: Incomplete
    composite_index_count: Incomplete

class KindStat(BaseKindStatistic):
    STORED_KIND_NAME: str
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete
    composite_index_bytes: Incomplete
    composite_index_count: Incomplete

class KindRootEntityStat(BaseKindStatistic):
    STORED_KIND_NAME: str

class KindNonRootEntityStat(BaseKindStatistic):
    STORED_KIND_NAME: str

class PropertyTypeStat(BaseStatistic):
    STORED_KIND_NAME: str
    property_type: Incomplete
    entity_bytes: Incomplete
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete

class KindPropertyTypeStat(BaseKindStatistic):
    STORED_KIND_NAME: str
    property_type: Incomplete
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete

class KindPropertyNameStat(BaseKindStatistic):
    STORED_KIND_NAME: str
    property_name: Incomplete
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete

class KindPropertyNamePropertyTypeStat(BaseKindStatistic):
    STORED_KIND_NAME: str
    property_type: Incomplete
    property_name: Incomplete
    builtin_index_bytes: Incomplete
    builtin_index_count: Incomplete

class KindCompositeIndexStat(BaseStatistic):
    STORED_KIND_NAME: str
    index_id: Incomplete
    kind_name: Incomplete

class NamespaceGlobalStat(GlobalStat):
    STORED_KIND_NAME: str

class NamespaceKindStat(KindStat):
    STORED_KIND_NAME: str

class NamespaceKindRootEntityStat(KindRootEntityStat):
    STORED_KIND_NAME: str

class NamespaceKindNonRootEntityStat(KindNonRootEntityStat):
    STORED_KIND_NAME: str

class NamespacePropertyTypeStat(PropertyTypeStat):
    STORED_KIND_NAME: str

class NamespaceKindPropertyTypeStat(KindPropertyTypeStat):
    STORED_KIND_NAME: str

class NamespaceKindPropertyNameStat(KindPropertyNameStat):
    STORED_KIND_NAME: str

class NamespaceKindPropertyNamePropertyTypeStat(KindPropertyNamePropertyTypeStat):
    STORED_KIND_NAME: str

class NamespaceKindCompositeIndexStat(KindCompositeIndexStat):
    STORED_KIND_NAME: str
