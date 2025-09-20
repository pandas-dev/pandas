from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from moto.utilities.utils import get_partition

from .exceptions import ValidationError


def make_arn_for_compute_env(account_id: str, name: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:batch:{region_name}:{account_id}:compute-environment/{name}"


def make_arn_for_job_queue(account_id: str, name: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:batch:{region_name}:{account_id}:job-queue/{name}"


def make_arn_for_job(account_id: str, job_id: str, region_name: str) -> str:
    return f"arn:{get_partition(region_name)}:batch:{region_name}:{account_id}:job/{job_id}"


def make_arn_for_task_def(
    account_id: str, name: str, revision: int, region_name: str
) -> str:
    return f"arn:{get_partition(region_name)}:batch:{region_name}:{account_id}:job-definition/{name}:{revision}"


def lowercase_first_key(some_dict: Dict[str, Any]) -> Dict[str, Any]:
    new_dict: Dict[str, Any] = {}
    for key, value in some_dict.items():
        new_key = key[0].lower() + key[1:]
        try:
            if isinstance(value, dict):
                new_dict[new_key] = lowercase_first_key(value)
            elif all([isinstance(v, dict) for v in value]):
                new_dict[new_key] = [lowercase_first_key(v) for v in value]
            else:
                new_dict[new_key] = value
        except TypeError:
            new_dict[new_key] = value

    return new_dict


def validate_job_status(target_job_status: str, valid_job_statuses: List[str]) -> None:
    if target_job_status not in valid_job_statuses:
        raise ValidationError(
            (
                "1 validation error detected: Value at 'current_status' failed "
                "to satisfy constraint: Member must satisfy enum value set: {valid_statues}"
            ).format(valid_statues=valid_job_statuses)
        )


class JobStatus(str, Enum):
    SUBMITTED = "SUBMITTED"
    PENDING = "PENDING"
    RUNNABLE = "RUNNABLE"
    STARTING = "STARTING"
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"

    @classmethod
    def job_statuses(self) -> List[str]:
        return sorted([item.value for item in JobStatus])

    @classmethod
    def is_job_already_started(self, current_status: str) -> bool:
        validate_job_status(current_status, JobStatus.job_statuses())
        return current_status not in [
            JobStatus.SUBMITTED,
            JobStatus.PENDING,
            JobStatus.RUNNABLE,
            JobStatus.STARTING,
        ]

    @classmethod
    def is_job_before_starting(self, current_status: str) -> bool:
        validate_job_status(current_status, JobStatus.job_statuses())
        return current_status in [
            JobStatus.SUBMITTED,
            JobStatus.PENDING,
            JobStatus.RUNNABLE,
        ]

    @classmethod
    def status_transitions(self) -> List[Tuple[Optional[str], str]]:
        return [
            (JobStatus.SUBMITTED.value, JobStatus.PENDING.value),
            (JobStatus.PENDING.value, JobStatus.RUNNABLE.value),
            (JobStatus.RUNNABLE.value, JobStatus.STARTING),
            (JobStatus.STARTING.value, JobStatus.RUNNING.value),
        ]
