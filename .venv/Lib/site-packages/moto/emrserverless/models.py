"""EMRServerlessBackend class with methods for supported APIs."""

import inspect
import re
from datetime import datetime
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_without_milliseconds
from moto.emrcontainers.utils import paginated_list
from moto.utilities.utils import get_partition

from .exceptions import (
    AccessDeniedException,
    ResourceNotFoundException,
    ValidationException,
)
from .utils import (
    default_auto_start_configuration,
    default_auto_stop_configuration,
    random_appplication_id,
    random_job_id,
)

APPLICATION_ARN_TEMPLATE = "arn:{partition}:emr-serverless:{region}:{account_id}:/applications/{application_id}"
JOB_RUN_ARN_TEMPLATE = "arn:{partition}:emr-serverless:{region}:{account_id}:/applications/{application_id}/jobruns/{job_run_id}"

# Defaults used for creating an EMR Serverless application
APPLICATION_STATUS = "STARTED"
JOB_STATUS = "SUCCESS"


class FakeApplication(BaseModel):
    def __init__(
        self,
        name: str,
        release_label: str,
        application_type: str,
        client_token: str,
        account_id: str,
        region_name: str,
        initial_capacity: str,
        maximum_capacity: str,
        tags: Dict[str, str],
        auto_start_configuration: str,
        auto_stop_configuration: str,
        network_configuration: str,
    ):
        # Provided parameters
        self.name = name
        self.release_label = release_label
        self.application_type = application_type.capitalize()
        self.client_token = client_token
        self.initial_capacity = initial_capacity
        self.maximum_capacity = maximum_capacity
        self.auto_start_configuration = (
            auto_start_configuration or default_auto_start_configuration()
        )
        self.auto_stop_configuration = (
            auto_stop_configuration or default_auto_stop_configuration()
        )
        self.network_configuration = network_configuration
        self.tags: Dict[str, str] = tags or {}

        # Service-generated-parameters
        self.id = random_appplication_id()
        self.arn = APPLICATION_ARN_TEMPLATE.format(
            partition="aws",
            region=region_name,
            account_id=account_id,
            application_id=self.id,
        )
        self.state = APPLICATION_STATUS
        self.state_details = ""
        self.created_at = iso_8601_datetime_without_milliseconds(
            datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
        )
        self.updated_at = self.created_at

    def __iter__(self) -> Iterator[Tuple[str, Any]]:
        yield "applicationId", self.id
        yield "name", self.name
        yield "arn", self.arn
        yield (
            "autoStartConfig",
            self.auto_start_configuration,
        )
        yield (
            "autoStopConfig",
            self.auto_stop_configuration,
        )

    def to_dict(self) -> Dict[str, Any]:
        """
        Dictionary representation of an EMR Serverless Application.
        When used in `list-applications`, capacity, auto-start/stop configs, and tags are not returned. https://docs.aws.amazon.com/emr-serverless/latest/APIReference/API_ListApplications.html
        When used in `get-application`, more details are returned. https://docs.aws.amazon.com/emr-serverless/latest/APIReference/API_GetApplication.html#API_GetApplication_ResponseSyntax
        """
        caller_methods = inspect.stack()[1].function
        caller_methods_type = caller_methods.split("_")[0]

        if caller_methods_type in ["get", "update"]:
            response = {
                "applicationId": self.id,
                "name": self.name,
                "arn": self.arn,
                "releaseLabel": self.release_label,
                "type": self.application_type,
                "state": self.state,
                "stateDetails": self.state_details,
                "createdAt": self.created_at,
                "updatedAt": self.updated_at,
                "autoStartConfiguration": self.auto_start_configuration,
                "autoStopConfiguration": self.auto_stop_configuration,
                "tags": self.tags,
            }
        else:
            response = {
                "id": self.id,
                "name": self.name,
                "arn": self.arn,
                "releaseLabel": self.release_label,
                "type": self.application_type,
                "state": self.state,
                "stateDetails": self.state_details,
                "createdAt": self.created_at,
                "updatedAt": self.updated_at,
            }

        if self.network_configuration:
            response.update({"networkConfiguration": self.network_configuration})
        if self.initial_capacity:
            response.update({"initialCapacity": self.initial_capacity})
        if self.maximum_capacity:
            response.update({"maximumCapacity": self.maximum_capacity})

        return response


class FakeJobRun(BaseModel):
    def __init__(
        self,
        application_id: str,
        client_token: str,
        execution_role_arn: str,
        account_id: str,
        region_name: str,
        release_label: str,
        application_type: str,
        job_driver: Optional[Dict[str, Dict[str, Union[str, List[str]]]]],
        configuration_overrides: Optional[Dict[str, Union[List[Any], Dict[str, Any]]]],
        tags: Optional[Dict[str, str]],
        network_configuration: Optional[Dict[str, List[str]]],
        execution_timeout_minutes: Optional[int],
        name: Optional[str],
    ):
        self.name = name
        self.application_id = application_id
        self.client_token = client_token
        self.execution_role_arn = execution_role_arn
        self.job_driver = job_driver
        self.configuration_overrides = configuration_overrides
        self.network_configuration = network_configuration
        self.execution_timeout_minutes = execution_timeout_minutes or 720

        # Service-generated-parameters
        self.id = random_job_id()

        self.arn = JOB_RUN_ARN_TEMPLATE.format(
            partition="aws",
            account_id=account_id,
            application_id=self.application_id,
            region=region_name,
            job_run_id=self.id,
        )

        self.release_label = release_label
        self.application_type = application_type

        self.state = JOB_STATUS
        self.state_details: Optional[str] = None

        self.created_by: Optional[str] = None

        self.created_at: str = iso_8601_datetime_without_milliseconds(
            datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
        )
        self.updated_at: str = self.created_at

        self.total_execution_duration_seconds: int = 0
        self.billed_resource_utilization: Dict[str, float] = {
            "vCPUHour": 0.0,
            "memoryGBHour": 0.0,
            "storageGBHour": 0.0,
        }

        self.tags = tags

    def to_dict(self, caller_methods_type: str) -> Dict[str, Any]:
        # The response structure is different for get/update and list
        if caller_methods_type in ["get", "update"]:
            response = {
                "applicationId": self.application_id,
                "jobRunId": self.id,
                "name": self.name,
                "arn": self.arn,
                "createdBy": self.created_by,
                "createdAt": self.created_at,
                "updatedAt": self.updated_at,
                "executionRole": self.execution_role_arn,
                "state": self.state,
                "stateDetails": self.state_details,
                "releaseLabel": self.release_label,
                "configurationOverrides": self.configuration_overrides,
                "jobDriver": self.job_driver,
                "tags": self.tags,
                "networkConfiguration": self.network_configuration,
                "totalExecutionDurationSeconds": self.total_execution_duration_seconds,
                "executionTimeoutMinutes": self.execution_timeout_minutes,
                "billedResourceUtilization": self.billed_resource_utilization,
            }
        else:
            response = {
                "applicationId": self.application_id,
                "id": self.id,
                "name": self.name,
                "arn": self.arn,
                "createdBy": self.created_by,
                "createdAt": self.created_at,
                "updatedAt": self.updated_at,
                "executionRole": self.execution_role_arn,
                "state": self.state,
                "stateDetails": self.state_details,
                "releaseLabel": self.release_label,
                "type": self.application_type,
            }

        return response


class EMRServerlessBackend(BaseBackend):
    """Implementation of EMRServerless APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.region_name = region_name
        self.partition = get_partition(region_name)
        self.applications: Dict[str, FakeApplication] = dict()
        self.job_runs: Dict[str, List[FakeJobRun]] = (
            dict()
        )  # {application_id: [job_run1, job_run2]}

    def create_application(
        self,
        name: str,
        release_label: str,
        application_type: str,
        client_token: str,
        initial_capacity: str,
        maximum_capacity: str,
        tags: Dict[str, str],
        auto_start_configuration: str,
        auto_stop_configuration: str,
        network_configuration: str,
    ) -> FakeApplication:
        if application_type not in ["HIVE", "SPARK"]:
            raise ValidationException(f"Unsupported engine {application_type}")

        if not re.match(r"emr-[0-9]{1}\.[0-9]{1,2}\.0(" "|-[0-9]{8})", release_label):
            raise ValidationException(
                f"Type '{application_type}' is not supported for release label '{release_label}' or release label does not exist"
            )

        application = FakeApplication(
            name=name,
            release_label=release_label,
            application_type=application_type,
            account_id=self.account_id,
            region_name=self.region_name,
            client_token=client_token,
            initial_capacity=initial_capacity,
            maximum_capacity=maximum_capacity,
            tags=tags,
            auto_start_configuration=auto_start_configuration,
            auto_stop_configuration=auto_stop_configuration,
            network_configuration=network_configuration,
        )
        self.applications[application.id] = application
        return application

    def delete_application(self, application_id: str) -> None:
        if application_id not in self.applications.keys():
            raise ResourceNotFoundException(application_id)

        if self.applications[application_id].state not in ["CREATED", "STOPPED"]:
            raise ValidationException(
                f"Application {application_id} must be in one of the following statuses [CREATED, STOPPED]. "
                f"Current status: {self.applications[application_id].state}"
            )
        self.applications[application_id].state = "TERMINATED"

    def get_application(self, application_id: str) -> Dict[str, Any]:
        if application_id not in self.applications.keys():
            raise ResourceNotFoundException(application_id)

        return self.applications[application_id].to_dict()

    def list_applications(
        self, next_token: Optional[str], max_results: int, states: Optional[List[str]]
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        applications = [
            application.to_dict() for application in self.applications.values()
        ]
        if states:
            applications = [
                application
                for application in applications
                if application["state"] in states
            ]
        sort_key = "name"
        return paginated_list(applications, sort_key, max_results, next_token)

    def start_application(self, application_id: str) -> None:
        if application_id not in self.applications.keys():
            raise ResourceNotFoundException(application_id)
        self.applications[application_id].state = "STARTED"

    def stop_application(self, application_id: str) -> None:
        if application_id not in self.applications.keys():
            raise ResourceNotFoundException(application_id)
        self.applications[application_id].state = "STOPPED"

    def update_application(
        self,
        application_id: str,
        initial_capacity: Optional[str],
        maximum_capacity: Optional[str],
        auto_start_configuration: Optional[str],
        auto_stop_configuration: Optional[str],
        network_configuration: Optional[str],
    ) -> Dict[str, Any]:
        if application_id not in self.applications.keys():
            raise ResourceNotFoundException(application_id)

        if self.applications[application_id].state not in ["CREATED", "STOPPED"]:
            raise ValidationException(
                f"Application {application_id} must be in one of the following statuses [CREATED, STOPPED]. "
                f"Current status: {self.applications[application_id].state}"
            )

        if initial_capacity:
            self.applications[application_id].initial_capacity = initial_capacity

        if maximum_capacity:
            self.applications[application_id].maximum_capacity = maximum_capacity

        if auto_start_configuration:
            self.applications[
                application_id
            ].auto_start_configuration = auto_start_configuration

        if auto_stop_configuration:
            self.applications[
                application_id
            ].auto_stop_configuration = auto_stop_configuration

        if network_configuration:
            self.applications[
                application_id
            ].network_configuration = network_configuration

        self.applications[
            application_id
        ].updated_at = iso_8601_datetime_without_milliseconds(
            datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
        )

        return self.applications[application_id].to_dict()

    def start_job_run(
        self,
        application_id: str,
        client_token: str,
        execution_role_arn: str,
        job_driver: Optional[Dict[str, Dict[str, Union[str, List[str]]]]],
        configuration_overrides: Optional[Dict[str, Union[List[Any], Dict[str, Any]]]],
        tags: Optional[Dict[str, str]],
        execution_timeout_minutes: Optional[int],
        name: Optional[str],
    ) -> FakeJobRun:
        role_account_id = execution_role_arn.split(":")[4]
        if role_account_id != self.account_id:
            raise AccessDeniedException("Cross-account pass role is not allowed.")

        if execution_timeout_minutes and execution_timeout_minutes < 5:
            raise ValidationException("RunTimeout must be at least 5 minutes.")

        application_resp = self.get_application(application_id)
        job_run = FakeJobRun(
            application_id=application_id,
            client_token=client_token,
            execution_role_arn=execution_role_arn,
            account_id=self.account_id,
            region_name=self.region_name,
            release_label=application_resp["releaseLabel"],
            application_type=application_resp["type"],
            job_driver=job_driver,
            configuration_overrides=configuration_overrides,
            tags=tags,
            network_configuration=application_resp.get("networkConfiguration"),
            execution_timeout_minutes=execution_timeout_minutes,
            name=name,
        )

        if application_resp["state"] == "TERMINATED":
            raise ValidationException(
                f"Application {application_id} is terminated. Cannot start job run."
            )

        if application_id not in self.job_runs:
            self.job_runs[application_id] = []
        self.job_runs[application_id].append(job_run)

        return job_run

    def get_job_run(self, application_id: str, job_run_id: str) -> FakeJobRun:
        if application_id not in self.job_runs.keys():
            raise ResourceNotFoundException(application_id, "Application")
        job_run_ids = [job_run.id for job_run in self.job_runs[application_id]]
        if job_run_id not in job_run_ids:
            raise ResourceNotFoundException(job_run_id, "JobRun")

        filtered_job_runs = [
            job_run
            for job_run in self.job_runs[application_id]
            if job_run.id == job_run_id
        ]
        assert len(filtered_job_runs) == 1
        job_run: FakeJobRun = filtered_job_runs[0]

        return job_run

    def cancel_job_run(self, application_id: str, job_run_id: str) -> Tuple[str, str]:
        # implement here
        if application_id not in self.job_runs.keys():
            raise ResourceNotFoundException(application_id, "Application")
        job_run_ids = [job_run.id for job_run in self.job_runs[application_id]]
        if job_run_id not in job_run_ids:
            raise ResourceNotFoundException(job_run_id, "JobRun")

        self.job_runs[application_id][job_run_ids.index(job_run_id)].state = "CANCELLED"

        return application_id, job_run_id

    def list_job_runs(
        self,
        application_id: str,
        max_results: int,
        next_token: Optional[str],
        created_at_after: Optional[str],
        created_at_before: Optional[str],
        states: Optional[List[str]],
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        if application_id not in self.job_runs.keys():
            raise ResourceNotFoundException(application_id, "Application")
        job_runs = self.job_runs[application_id]
        if states:
            job_runs = [job_run for job_run in job_runs if job_run.state in states]
        if created_at_after:
            job_runs = [
                job_run
                for job_run in job_runs
                if datetime.strptime(job_run.created_at, "%Y-%m-%dT%H:%M:%SZ")
                > datetime.strptime(created_at_after, "%Y-%m-%dT%H:%M:%SZ")
            ]
        if created_at_before:
            job_runs = [
                job_run
                for job_run in job_runs
                if datetime.strptime(job_run.created_at, "%Y-%m-%dT%H:%M:%SZ")
                < datetime.strptime(created_at_before, "%Y-%m-%dT%H:%M:%SZ")
            ]

        job_run_dicts = [job_run.to_dict("list") for job_run in job_runs]

        sort_key = "createdAt"
        return paginated_list(job_run_dicts, sort_key, max_results, next_token)


emrserverless_backends = BackendDict(EMRServerlessBackend, "emr-serverless")
