from typing import Any, Dict, Iterable, List


def match_task_arn(task: Dict[str, Any], arns: List[str]) -> bool:
    return task["ReplicationTaskArn"] in arns


def match_task_id(task: Dict[str, Any], ids: List[str]) -> bool:
    return task["ReplicationTaskIdentifier"] in ids


def match_task_migration_type(task: Dict[str, Any], migration_types: List[str]) -> bool:
    return task["MigrationType"] in migration_types


def match_task_endpoint_arn(task: Dict[str, Any], endpoint_arns: List[str]) -> bool:
    return (
        task["SourceEndpointArn"] in endpoint_arns
        or task["TargetEndpointArn"] in endpoint_arns
    )


def match_task_replication_instance_arn(
    task: Dict[str, Any], replication_instance_arns: List[str]
) -> bool:
    return task["ReplicationInstanceArn"] in replication_instance_arns


task_filter_functions = {
    "replication-task-arn": match_task_arn,
    "replication-task-id": match_task_id,
    "migration-type": match_task_migration_type,
    "endpoint-arn": match_task_endpoint_arn,
    "replication-instance-arn": match_task_replication_instance_arn,
}


def filter_tasks(tasks: Iterable[Any], filters: List[Dict[str, Any]]) -> Any:
    matching_tasks = tasks

    for f in filters:
        filter_function = task_filter_functions.get(f["Name"])

        if not filter_function:
            continue

        matching_tasks = filter(
            lambda task: filter_function(task, f["Values"]), matching_tasks
        )

    return matching_tasks
