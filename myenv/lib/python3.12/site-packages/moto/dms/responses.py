import json

from moto.core.responses import BaseResponse

from .models import DatabaseMigrationServiceBackend, dms_backends


class DatabaseMigrationServiceResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="dms")

    @property
    def dms_backend(self) -> DatabaseMigrationServiceBackend:
        return dms_backends[self.current_account][self.region]

    def create_replication_task(self) -> str:
        replication_task_identifier = self._get_param("ReplicationTaskIdentifier")
        source_endpoint_arn = self._get_param("SourceEndpointArn")
        target_endpoint_arn = self._get_param("TargetEndpointArn")
        replication_instance_arn = self._get_param("ReplicationInstanceArn")
        migration_type = self._get_param("MigrationType")
        table_mappings = self._get_param("TableMappings")
        replication_task_settings = self._get_param("ReplicationTaskSettings")
        replication_task = self.dms_backend.create_replication_task(
            replication_task_identifier=replication_task_identifier,
            source_endpoint_arn=source_endpoint_arn,
            target_endpoint_arn=target_endpoint_arn,
            replication_instance_arn=replication_instance_arn,
            migration_type=migration_type,
            table_mappings=table_mappings,
            replication_task_settings=replication_task_settings,
        )

        return json.dumps({"ReplicationTask": replication_task.to_dict()})

    def start_replication_task(self) -> str:
        replication_task_arn = self._get_param("ReplicationTaskArn")
        replication_task = self.dms_backend.start_replication_task(
            replication_task_arn=replication_task_arn
        )

        return json.dumps({"ReplicationTask": replication_task.to_dict()})

    def stop_replication_task(self) -> str:
        replication_task_arn = self._get_param("ReplicationTaskArn")
        replication_task = self.dms_backend.stop_replication_task(
            replication_task_arn=replication_task_arn
        )

        return json.dumps({"ReplicationTask": replication_task.to_dict()})

    def delete_replication_task(self) -> str:
        replication_task_arn = self._get_param("ReplicationTaskArn")
        replication_task = self.dms_backend.delete_replication_task(
            replication_task_arn=replication_task_arn
        )

        return json.dumps({"ReplicationTask": replication_task.to_dict()})

    def describe_replication_tasks(self) -> str:
        filters = self._get_list_prefix("Filters.member")
        max_records = self._get_int_param("MaxRecords")
        replication_tasks = self.dms_backend.describe_replication_tasks(
            filters=filters, max_records=max_records
        )

        return json.dumps(
            dict(ReplicationTasks=[t.to_dict() for t in replication_tasks])
        )
