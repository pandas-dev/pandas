from typing import Any, Dict, Iterable, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ResourceNotFound


class TimestreamTable(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        table_name: str,
        db_name: str,
        retention_properties: Dict[str, int],
        magnetic_store_write_properties: Dict[str, Any],
        schema: Dict[str, Any],
    ):
        self.region_name = region_name
        self.name = table_name
        self.db_name = db_name
        self.retention_properties = retention_properties or {
            "MemoryStoreRetentionPeriodInHours": 6,
            "MagneticStoreRetentionPeriodInDays": 73000,
        }
        self.magnetic_store_write_properties = magnetic_store_write_properties or {}
        self.schema = schema or {
            "CompositePartitionKey": [
                {
                    "Type": "MEASURE",
                    "Name": "",
                    "EnforcementInRecord": "",
                }
            ]
        }
        self.records: List[Dict[str, Any]] = []
        self.arn = f"arn:{get_partition(self.region_name)}:timestream:{self.region_name}:{account_id}:database/{self.db_name}/table/{self.name}"

    def update(
        self,
        retention_properties: Dict[str, int],
        magnetic_store_write_properties: Dict[str, Any],
        schema: Dict[str, Any],
    ) -> None:
        self.retention_properties = retention_properties
        if magnetic_store_write_properties is not None:
            self.magnetic_store_write_properties = magnetic_store_write_properties
        if schema is not None:
            self.schema = schema

    def write_records(self, records: List[Dict[str, Any]]) -> None:
        self.records.extend(records)

    def description(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "TableName": self.name,
            "DatabaseName": self.db_name,
            "TableStatus": "ACTIVE",
            "RetentionProperties": self.retention_properties,
            "MagneticStoreWriteProperties": self.magnetic_store_write_properties,
            "Schema": self.schema,
        }


class TimestreamDatabase(BaseModel):
    def __init__(
        self, account_id: str, region_name: str, database_name: str, kms_key_id: str
    ):
        self.account_id = account_id
        self.region_name = region_name
        self.name = database_name
        self.kms_key_id = (
            kms_key_id
            or f"arn:{get_partition(region_name)}:kms:{region_name}:{account_id}:key/default_key"
        )
        self.arn = f"arn:{get_partition(self.region_name)}:timestream:{self.region_name}:{account_id}:database/{self.name}"
        self.created_on = unix_time()
        self.updated_on = unix_time()
        self.tables: Dict[str, TimestreamTable] = dict()

    def update(self, kms_key_id: str) -> None:
        self.kms_key_id = kms_key_id

    def create_table(
        self,
        table_name: str,
        retention_properties: Dict[str, int],
        magnetic_store_write_properties: Dict[str, Any],
        schema: Dict[str, Any],
    ) -> TimestreamTable:
        table = TimestreamTable(
            account_id=self.account_id,
            region_name=self.region_name,
            table_name=table_name,
            db_name=self.name,
            retention_properties=retention_properties,
            magnetic_store_write_properties=magnetic_store_write_properties,
            schema=schema,
        )
        self.tables[table_name] = table
        return table

    def update_table(
        self,
        table_name: str,
        retention_properties: Dict[str, int],
        magnetic_store_write_properties: Dict[str, Any],
        schema: Dict[str, Any],
    ) -> TimestreamTable:
        table = self.tables[table_name]
        table.update(
            retention_properties=retention_properties,
            magnetic_store_write_properties=magnetic_store_write_properties,
            schema=schema,
        )
        return table

    def delete_table(self, table_name: str) -> None:
        self.tables.pop(table_name, None)

    def describe_table(self, table_name: str) -> TimestreamTable:
        if table_name not in self.tables:
            raise ResourceNotFound(f"The table {table_name} does not exist.")
        return self.tables[table_name]

    def list_tables(self) -> Iterable[TimestreamTable]:
        return self.tables.values()

    def description(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "DatabaseName": self.name,
            "TableCount": len(self.tables.keys()),
            "KmsKeyId": self.kms_key_id,
            "CreationTime": self.created_on,
            "LastUpdatedTime": self.updated_on,
        }


class TimestreamWriteBackend(BaseBackend):
    """
    When using the decorators, you can use the following internal API to verify records have arrived:

    .. sourcecode:: python

        from moto.core import DEFAULT_ACCOUNT_ID as ACCOUNT_ID
        from moto.timestreamwrite.models import timestreamwrite_backends

        backend = timestreamwrite_backends[ACCOUNT_ID]["us-east-1"]
        records = backend.databases["mydatabase"].tables["mytable"].records

    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.databases: Dict[str, TimestreamDatabase] = dict()
        self.tagging_service = TaggingService()

    def create_database(
        self, database_name: str, kms_key_id: str, tags: List[Dict[str, str]]
    ) -> TimestreamDatabase:
        database = TimestreamDatabase(
            self.account_id, self.region_name, database_name, kms_key_id
        )
        self.databases[database_name] = database
        self.tagging_service.tag_resource(database.arn, tags)
        return database

    def delete_database(self, database_name: str) -> None:
        del self.databases[database_name]

    def describe_database(self, database_name: str) -> TimestreamDatabase:
        if database_name not in self.databases:
            raise ResourceNotFound(f"The database {database_name} does not exist.")
        return self.databases[database_name]

    def list_databases(self) -> Iterable[TimestreamDatabase]:
        return self.databases.values()

    def update_database(
        self, database_name: str, kms_key_id: str
    ) -> TimestreamDatabase:
        database = self.databases[database_name]
        database.update(kms_key_id=kms_key_id)
        return database

    def create_table(
        self,
        database_name: str,
        table_name: str,
        retention_properties: Dict[str, int],
        tags: List[Dict[str, str]],
        magnetic_store_write_properties: Dict[str, Any],
        schema: Dict[str, Any],
    ) -> TimestreamTable:
        database = self.describe_database(database_name)
        table = database.create_table(
            table_name,
            retention_properties,
            magnetic_store_write_properties,
            schema,
        )
        self.tagging_service.tag_resource(table.arn, tags)
        return table

    def delete_table(self, database_name: str, table_name: str) -> None:
        database = self.describe_database(database_name)
        database.delete_table(table_name)

    def describe_table(self, database_name: str, table_name: str) -> TimestreamTable:
        database = self.describe_database(database_name)
        return database.describe_table(table_name)

    def list_tables(self, database_name: str) -> Iterable[TimestreamTable]:
        database = self.describe_database(database_name)
        return database.list_tables()

    def update_table(
        self,
        database_name: str,
        table_name: str,
        retention_properties: Dict[str, int],
        magnetic_store_write_properties: Dict[str, Any],
        schema: Dict[str, Any],
    ) -> TimestreamTable:
        database = self.describe_database(database_name)
        return database.update_table(
            table_name,
            retention_properties,
            magnetic_store_write_properties,
            schema,
        )

    def write_records(
        self, database_name: str, table_name: str, records: List[Dict[str, Any]]
    ) -> None:
        database = self.describe_database(database_name)
        table = database.describe_table(table_name)
        table.write_records(records)

    def describe_endpoints(self) -> Dict[str, List[Dict[str, Any]]]:
        # https://docs.aws.amazon.com/timestream/latest/developerguide/Using-API.endpoint-discovery.how-it-works.html
        # Usually, the address look like this:
        # ingest-cell1.timestream.us-east-1.amazonaws.com
        # Where 'cell1' can be any number, 'cell2', 'cell3', etc - whichever endpoint happens to be available for that particular account
        # We don't implement a cellular architecture in Moto though, so let's keep it simple
        return {
            "Endpoints": [
                {
                    "Address": f"ingest.timestream.{self.region_name}.amazonaws.com",
                    "CachePeriodInMinutes": 1440,
                }
            ]
        }

    def list_tags_for_resource(
        self, resource_arn: str
    ) -> Dict[str, List[Dict[str, str]]]:
        return self.tagging_service.list_tags_for_resource(resource_arn)

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagging_service.tag_resource(resource_arn, tags)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagging_service.untag_resource_using_names(resource_arn, tag_keys)


# Boto does not return any regions at the time of writing (20/10/2021)
# Hardcoding the known regions for now
# Thanks, Jeff
timestreamwrite_backends = BackendDict(
    TimestreamWriteBackend,
    "timestream-write",
    additional_regions=[
        "us-east-1",
        "us-east-2",
        "us-west-2",
        "eu-central-1",
        "eu-west-1",
        "ap-southeast-2",
        "ap-northeast-1",
    ],
)
