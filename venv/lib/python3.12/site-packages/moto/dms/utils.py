from __future__ import annotations

import random
import string
from typing import TYPE_CHECKING, Any, Dict, List

if TYPE_CHECKING:
    from moto.dms.models import FakeReplicationTask


def random_id(uppercase: bool = True, length: int = 13) -> str:
    ascii_set = string.ascii_uppercase if uppercase else string.ascii_lowercase
    chars = list(range(10)) + list(ascii_set)
    resource_id = random.choice(ascii_set) + "".join(
        str(random.choice(chars)) for _ in range(length - 1)
    )
    return resource_id


def match_task_arn(task: FakeReplicationTask, arns: List[str]) -> bool:
    return task.arn in arns


def match_task_id(task: FakeReplicationTask, ids: List[str]) -> bool:
    return task.id in ids


def match_task_migration_type(
    task: FakeReplicationTask, migration_types: List[str]
) -> bool:
    return task.migration_type in migration_types


def match_task_endpoint_arn(
    task: FakeReplicationTask, endpoint_arns: List[str]
) -> bool:
    return (
        task.source_endpoint_arn in endpoint_arns
        or task.target_endpoint_arn in endpoint_arns
    )


def match_task_replication_instance_arn(
    task: FakeReplicationTask, replication_instance_arns: List[str]
) -> bool:
    return task.replication_instance_arn in replication_instance_arns


task_filter_functions = {
    "replication-task-arn": match_task_arn,
    "replication-task-id": match_task_id,
    "migration-type": match_task_migration_type,
    "endpoint-arn": match_task_endpoint_arn,
    "replication-instance-arn": match_task_replication_instance_arn,
}


def filter_tasks(
    tasks: List[FakeReplicationTask], filters: List[Dict[str, Any]]
) -> Any:
    matching_tasks = tasks

    for f in filters:
        filter_function = task_filter_functions.get(f["Name"])

        if not filter_function:
            continue

        matching_tasks = list(
            filter(lambda task: filter_function(task, f["Values"]), matching_tasks)
        )

    return matching_tasks
