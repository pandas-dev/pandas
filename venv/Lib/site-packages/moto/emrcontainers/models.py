"""EMRContainersBackend class with methods for supported APIs."""

import re
from datetime import datetime
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import utcnow
from moto.utilities.utils import get_partition

from ..config.exceptions import ValidationException
from .exceptions import ResourceNotFoundException
from .utils import random_cluster_id, random_job_id

VIRTUAL_CLUSTER_ARN_TEMPLATE = "arn:{partition}:emr-containers:{region}:{account_id}:/virtualclusters/{virtual_cluster_id}"

JOB_ARN_TEMPLATE = "arn:{partition}:emr-containers:{region}:{account_id}:/virtualclusters/{virtual_cluster_id}/jobruns/{job_id}"

# Defaults used for creating a Virtual cluster
VIRTUAL_CLUSTER_STATUS = "RUNNING"
JOB_STATUS = "RUNNING"


class VirtualCluster(BaseModel):
    def __init__(
        self,
        name: str,
        container_provider: dict[str, Any],
        client_token: str,
        account_id: str,
        region_name: str,
        aws_partition: str,
        tags: dict[str, str] | None = None,
        virtual_cluster_id: str | None = None,
    ):
        self.id = virtual_cluster_id or random_cluster_id()

        self.name = name
        self.client_token = client_token
        self.arn = VIRTUAL_CLUSTER_ARN_TEMPLATE.format(
            partition=aws_partition,
            region=region_name,
            account_id=account_id,
            virtual_cluster_id=self.id,
        )
        self.state = VIRTUAL_CLUSTER_STATUS
        self.container_provider = container_provider
        self.container_provider_id = container_provider["id"]
        self.namespace = container_provider["info"]["eksInfo"]["namespace"]
        self.created_at = utcnow().replace(hour=0, minute=0, second=0, microsecond=0)
        self.tags = tags


class JobRun(BaseModel):
    def __init__(
        self,
        name: str,
        virtual_cluster_id: str,
        client_token: str,
        execution_role_arn: str,
        release_label: str,
        job_driver: str,
        configuration_overrides: dict[str, Any],
        account_id: str,
        region_name: str,
        aws_partition: str,
        tags: dict[str, str] | None,
    ):
        self.id = random_job_id()
        self.name = name
        self.virtual_cluster_id = virtual_cluster_id
        self.arn = JOB_ARN_TEMPLATE.format(
            partition=aws_partition,
            region=region_name,
            account_id=account_id,
            virtual_cluster_id=self.virtual_cluster_id,
            job_id=self.id,
        )
        self.state = JOB_STATUS
        self.client_token = client_token
        self.execution_role_arn = execution_role_arn
        self.release_label = release_label
        self.job_driver = job_driver
        self.configuration_overrides = configuration_overrides
        self.created_at = utcnow().replace(hour=0, minute=0, second=0, microsecond=0)
        self.created_by = None
        self.finished_at: datetime | None = None
        self.state_details: str | None = None
        self.failure_reason = None
        self.tags = tags


class EMRContainersBackend(BaseBackend):
    """Implementation of EMRContainers APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.virtual_clusters: dict[str, VirtualCluster] = {}
        self.virtual_cluster_count = 0
        self.job_runs: dict[str, JobRun] = {}
        self.job_count = 0
        self.partition = get_partition(region_name)

    def create_virtual_cluster(
        self,
        name: str,
        container_provider: dict[str, Any],
        client_token: str,
        tags: dict[str, str] | None = None,
    ) -> VirtualCluster:
        occupied_namespaces = [
            virtual_cluster.namespace
            for virtual_cluster in self.virtual_clusters.values()
        ]

        if container_provider["info"]["eksInfo"]["namespace"] in occupied_namespaces:
            raise ValidationException(
                "A virtual cluster already exists in the given namespace"
            )

        virtual_cluster = VirtualCluster(
            name=name,
            container_provider=container_provider,
            client_token=client_token,
            tags=tags,
            account_id=self.account_id,
            region_name=self.region_name,
            aws_partition=self.partition,
        )

        self.virtual_clusters[virtual_cluster.id] = virtual_cluster
        self.virtual_cluster_count += 1
        return virtual_cluster

    def delete_virtual_cluster(self, cluster_id: str) -> VirtualCluster:
        if cluster_id not in self.virtual_clusters:
            raise ValidationException("VirtualCluster does not exist")

        self.virtual_clusters[cluster_id].state = "TERMINATED"
        return self.virtual_clusters[cluster_id]

    def describe_virtual_cluster(self, cluster_id: str) -> VirtualCluster:
        if cluster_id not in self.virtual_clusters:
            raise ValidationException(f"Virtual cluster {cluster_id} doesn't exist.")

        return self.virtual_clusters[cluster_id]

    def list_virtual_clusters(
        self,
        container_provider_id: str,
        container_provider_type: str,
        created_after: datetime,
        created_before: datetime,
        states: list[str] | None,
    ) -> list[VirtualCluster]:
        virtual_clusters = list(self.virtual_clusters.values())

        if container_provider_id:
            virtual_clusters = [
                virtual_cluster
                for virtual_cluster in virtual_clusters
                if virtual_cluster.container_provider["id"] == container_provider_id
            ]

        if container_provider_type:
            virtual_clusters = [
                virtual_cluster
                for virtual_cluster in virtual_clusters
                if virtual_cluster.container_provider["type"] == container_provider_type
            ]

        if created_after:
            virtual_clusters = [
                virtual_cluster
                for virtual_cluster in virtual_clusters
                if virtual_cluster.created_at >= created_after
            ]

        if created_before:
            virtual_clusters = [
                virtual_cluster
                for virtual_cluster in virtual_clusters
                if virtual_cluster.created_at <= created_before
            ]

        if states:
            virtual_clusters = [
                virtual_cluster
                for virtual_cluster in virtual_clusters
                if virtual_cluster.state in states
            ]
        return virtual_clusters

    def start_job_run(
        self,
        name: str,
        virtual_cluster_id: str,
        client_token: str,
        execution_role_arn: str,
        release_label: str,
        job_driver: str,
        configuration_overrides: dict[str, Any],
        tags: dict[str, str],
    ) -> JobRun:
        if virtual_cluster_id not in self.virtual_clusters.keys():
            raise ResourceNotFoundException(
                f"Virtual cluster {virtual_cluster_id} doesn't exist."
            )

        if not re.match(
            r"emr-[0-9]{1}\.[0-9]{1,2}\.0-(latest|[0-9]{8})", release_label
        ):
            raise ResourceNotFoundException(f"Release {release_label} doesn't exist.")

        job_run = JobRun(
            name=name,
            virtual_cluster_id=virtual_cluster_id,
            client_token=client_token,
            execution_role_arn=execution_role_arn,
            release_label=release_label,
            job_driver=job_driver,
            configuration_overrides=configuration_overrides,
            tags=tags,
            account_id=self.account_id,
            region_name=self.region_name,
            aws_partition=self.partition,
        )

        self.job_runs[job_run.id] = job_run
        self.job_count += 1
        return job_run

    def cancel_job_run(self, job_run_id: str, virtual_cluster_id: str) -> JobRun:
        if not re.match(r"[a-z,A-Z,0-9]{19}", job_run_id):
            raise ValidationException("Invalid job run short id")

        if job_run_id not in self.job_runs.keys():
            raise ResourceNotFoundException(f"Job run {job_run_id} doesn't exist.")

        if virtual_cluster_id != self.job_runs[job_run_id].virtual_cluster_id:
            raise ResourceNotFoundException(f"Job run {job_run_id} doesn't exist.")

        if self.job_runs[job_run_id].state in [
            "FAILED",
            "CANCELLED",
            "CANCEL_PENDING",
            "COMPLETED",
        ]:
            raise ValidationException(
                f"Job run {job_run_id} is not in a cancellable state"
            )

        job_run = self.job_runs[job_run_id]
        job_run.state = "CANCELLED"
        job_run.finished_at = utcnow().replace(
            hour=0, minute=1, second=0, microsecond=0
        )
        job_run.state_details = "JobRun CANCELLED successfully."

        return job_run

    def list_job_runs(
        self,
        virtual_cluster_id: str,
        created_before: datetime,
        created_after: datetime,
        name: str,
        states: list[str] | None,
    ) -> list[JobRun]:
        jobs = list(self.job_runs.values())

        jobs = [job for job in jobs if job.virtual_cluster_id == virtual_cluster_id]

        if created_after:
            jobs = [job for job in jobs if job.created_at >= created_after]

        if created_before:
            jobs = [job for job in jobs if job.created_at <= created_before]

        if states:
            jobs = [job for job in jobs if job.state in states]

        if name:
            jobs = [job for job in jobs if job.name in name]

        return jobs

    def describe_job_run(self, job_run_id: str, virtual_cluster_id: str) -> JobRun:
        if not re.match(r"[a-z,A-Z,0-9]{19}", job_run_id):
            raise ValidationException("Invalid job run short id")

        if job_run_id not in self.job_runs.keys():
            raise ResourceNotFoundException(f"Job run {job_run_id} doesn't exist.")

        if virtual_cluster_id != self.job_runs[job_run_id].virtual_cluster_id:
            raise ResourceNotFoundException(f"Job run {job_run_id} doesn't exist.")

        return self.job_runs[job_run_id]


emrcontainers_backends = BackendDict(EMRContainersBackend, "emr-containers")
