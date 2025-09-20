import json

from moto.core.responses import BaseResponse
from moto.timestreamquery.models import TimestreamQueryBackend, timestreamquery_backends

from .models import TimestreamWriteBackend, timestreamwrite_backends


class TimestreamWriteResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="timestream-write")

    @property
    def timestreamquery_backend(self) -> TimestreamQueryBackend:
        return timestreamquery_backends[self.current_account][self.region]

    @property
    def timestreamwrite_backend(self) -> TimestreamWriteBackend:
        """Return backend instance specific for this region."""
        return timestreamwrite_backends[self.current_account][self.region]

    def create_database(self) -> str:
        database_name = self._get_param("DatabaseName")
        kms_key_id = self._get_param("KmsKeyId")
        tags = self._get_param("Tags")
        database = self.timestreamwrite_backend.create_database(
            database_name=database_name, kms_key_id=kms_key_id, tags=tags
        )
        return json.dumps(dict(Database=database.description()))

    def delete_database(self) -> str:
        database_name = self._get_param("DatabaseName")
        self.timestreamwrite_backend.delete_database(database_name=database_name)
        return "{}"

    def describe_database(self) -> str:
        database_name = self._get_param("DatabaseName")
        database = self.timestreamwrite_backend.describe_database(
            database_name=database_name
        )
        return json.dumps(dict(Database=database.description()))

    def update_database(self) -> str:
        database_name = self._get_param("DatabaseName")
        kms_key_id = self._get_param("KmsKeyId")
        database = self.timestreamwrite_backend.update_database(
            database_name, kms_key_id
        )
        return json.dumps(dict(Database=database.description()))

    def list_databases(self) -> str:
        all_dbs = self.timestreamwrite_backend.list_databases()
        return json.dumps(dict(Databases=[db.description() for db in all_dbs]))

    def create_table(self) -> str:
        database_name = self._get_param("DatabaseName")
        table_name = self._get_param("TableName")
        retention_properties = self._get_param("RetentionProperties")
        tags = self._get_param("Tags")
        magnetic_store_write_properties = self._get_param(
            "MagneticStoreWriteProperties"
        )
        schema = self._get_param("Schema")
        table = self.timestreamwrite_backend.create_table(
            database_name,
            table_name,
            retention_properties,
            tags,
            magnetic_store_write_properties,
            schema=schema,
        )
        return json.dumps(dict(Table=table.description()))

    def delete_table(self) -> str:
        database_name = self._get_param("DatabaseName")
        table_name = self._get_param("TableName")
        self.timestreamwrite_backend.delete_table(database_name, table_name)
        return "{}"

    def describe_table(self) -> str:
        database_name = self._get_param("DatabaseName")
        table_name = self._get_param("TableName")
        table = self.timestreamwrite_backend.describe_table(database_name, table_name)
        return json.dumps(dict(Table=table.description()))

    def list_tables(self) -> str:
        database_name = self._get_param("DatabaseName")
        tables = self.timestreamwrite_backend.list_tables(database_name)
        return json.dumps(dict(Tables=[t.description() for t in tables]))

    def update_table(self) -> str:
        database_name = self._get_param("DatabaseName")
        table_name = self._get_param("TableName")
        retention_properties = self._get_param("RetentionProperties")
        magnetic_store_write_properties = self._get_param(
            "MagneticStoreWriteProperties"
        )
        schema = self._get_param("Schema")
        table = self.timestreamwrite_backend.update_table(
            database_name,
            table_name,
            retention_properties,
            magnetic_store_write_properties,
            schema=schema,
        )
        return json.dumps(dict(Table=table.description()))

    def write_records(self) -> str:
        database_name = self._get_param("DatabaseName")
        table_name = self._get_param("TableName")
        records = self._get_param("Records")
        self.timestreamwrite_backend.write_records(database_name, table_name, records)
        resp = {
            "RecordsIngested": {
                "Total": len(records),
                "MemoryStore": len(records),
                "MagneticStore": 0,
            }
        }
        return json.dumps(resp)

    def describe_endpoints(self) -> str:
        resp = self.timestreamwrite_backend.describe_endpoints()
        return json.dumps(resp)

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_param("ResourceARN")
        tags = self.timestreamwrite_backend.list_tags_for_resource(resource_arn)
        return json.dumps(tags)

    def tag_resource(self) -> str:
        resource_arn = self._get_param("ResourceARN")
        tags = self._get_param("Tags")
        self.timestreamwrite_backend.tag_resource(resource_arn, tags)
        return "{}"

    def untag_resource(self) -> str:
        resource_arn = self._get_param("ResourceARN")
        tag_keys = self._get_param("TagKeys")
        self.timestreamwrite_backend.untag_resource(resource_arn, tag_keys)
        return "{}"

    # AWS uses the same path/headers for TimestreamWrite and TimestreamQuery
    # The only difference is the host, but we don't have access to that in ServerMode
    # Keep the `Query`-responses here
    # So we don't have to jump through all kinds of hoops to explain the difference between Query and Write to MotoServer

    def create_scheduled_query(self) -> str:
        name = self._get_param("Name")
        query_string = self._get_param("QueryString")
        schedule_configuration = self._get_param("ScheduleConfiguration")
        notification_configuration = self._get_param("NotificationConfiguration")
        target_configuration = self._get_param("TargetConfiguration")
        scheduled_query_execution_role_arn = self._get_param(
            "ScheduledQueryExecutionRoleArn"
        )
        tags = self._get_param("Tags")
        kms_key_id = self._get_param("KmsKeyId")
        error_report_configuration = self._get_param("ErrorReportConfiguration")
        scheduled_query = self.timestreamquery_backend.create_scheduled_query(
            name=name,
            query_string=query_string,
            schedule_configuration=schedule_configuration,
            notification_configuration=notification_configuration,
            target_configuration=target_configuration,
            scheduled_query_execution_role_arn=scheduled_query_execution_role_arn,
            tags=tags,
            kms_key_id=kms_key_id,
            error_report_configuration=error_report_configuration,
        )
        return json.dumps(dict(Arn=scheduled_query.arn))

    def delete_scheduled_query(self) -> str:
        scheduled_query_arn = self._get_param("ScheduledQueryArn")
        self.timestreamquery_backend.delete_scheduled_query(
            scheduled_query_arn=scheduled_query_arn,
        )
        return "{}"

    def update_scheduled_query(self) -> str:
        scheduled_query_arn = self._get_param("ScheduledQueryArn")
        state = self._get_param("State")
        self.timestreamquery_backend.update_scheduled_query(
            scheduled_query_arn=scheduled_query_arn,
            state=state,
        )
        return "{}"

    def query(self) -> str:
        query_string = self._get_param("QueryString")
        result = self.timestreamquery_backend.query(query_string=query_string)
        return json.dumps(result)

    def describe_scheduled_query(self) -> str:
        scheduled_query_arn = self._get_param("ScheduledQueryArn")
        scheduled_query = self.timestreamquery_backend.describe_scheduled_query(
            scheduled_query_arn=scheduled_query_arn,
        )
        return json.dumps(dict(ScheduledQuery=scheduled_query.description()))
