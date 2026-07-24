"""Handles incoming emrcontainers requests, invokes methods, returns responses."""

from moto.core.responses import ActionResult, BaseResponse, PaginatedResult

from .models import EMRContainersBackend, emrcontainers_backends

DEFAULT_CONTAINER_PROVIDER_TYPE = "EKS"


class EMRContainersResponse(BaseResponse):
    """Handler for EMRContainers requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="emr-containers")
        self.automated_parameter_parsing = True

    @property
    def emrcontainers_backend(self) -> EMRContainersBackend:
        """Return backend instance specific for this region."""
        return emrcontainers_backends[self.current_account][self.region]

    def create_virtual_cluster(self) -> ActionResult:
        name = self._get_param("name")
        container_provider = self._get_param("containerProvider")
        client_token = self._get_param("clientToken")
        tags = self._get_param("tags")
        virtual_cluster = self.emrcontainers_backend.create_virtual_cluster(
            name=name,
            container_provider=container_provider,
            client_token=client_token,
            tags=tags,
        )
        result = {
            "id": virtual_cluster.id,
            "name": virtual_cluster.name,
            "arn": virtual_cluster.arn,
        }
        return ActionResult(result)

    def delete_virtual_cluster(self) -> ActionResult:
        cluster_id = self._get_param("id")
        virtual_cluster = self.emrcontainers_backend.delete_virtual_cluster(
            cluster_id=cluster_id
        )
        return ActionResult({"id": virtual_cluster.id})

    def describe_virtual_cluster(self) -> ActionResult:
        cluster_id = self._get_param("id")
        virtual_cluster = self.emrcontainers_backend.describe_virtual_cluster(
            cluster_id=cluster_id
        )
        result = {"virtualCluster": virtual_cluster}
        return ActionResult(result)

    def list_virtual_clusters(self) -> ActionResult:
        container_provider_id = self._get_param("containerProviderId")
        container_provider_type = self._get_param(
            "containerProviderType", DEFAULT_CONTAINER_PROVIDER_TYPE
        )
        created_after = self._get_param("createdAfter")
        created_before = self._get_param("createdBefore")
        states = self.querystring.get("states", [])
        virtual_clusters = self.emrcontainers_backend.list_virtual_clusters(
            container_provider_id=container_provider_id,
            container_provider_type=container_provider_type,
            created_after=created_after,
            created_before=created_before,
            states=states,
        )
        result = {"virtualClusters": virtual_clusters}
        return PaginatedResult(result)

    def start_job_run(self) -> ActionResult:
        name = self._get_param("name")
        virtual_cluster_id = self._get_param("virtualClusterId")
        client_token = self._get_param("clientToken")
        execution_role_arn = self._get_param("executionRoleArn")
        release_label = self._get_param("releaseLabel")
        job_driver = self._get_param("jobDriver")
        configuration_overrides = self._get_param("configurationOverrides")
        tags = self._get_param("tags")
        job_run = self.emrcontainers_backend.start_job_run(
            name=name,
            virtual_cluster_id=virtual_cluster_id,
            client_token=client_token,
            execution_role_arn=execution_role_arn,
            release_label=release_label,
            job_driver=job_driver,
            configuration_overrides=configuration_overrides,
            tags=tags,
        )
        result = {
            "id": job_run.id,
            "name": job_run.name,
            "arn": job_run.arn,
            "virtualClusterId": job_run.virtual_cluster_id,
        }
        return ActionResult(result)

    def cancel_job_run(self) -> ActionResult:
        job_run_id = self._get_param("id")
        virtual_cluster_id = self._get_param("virtualClusterId")
        job_run = self.emrcontainers_backend.cancel_job_run(
            job_run_id=job_run_id, virtual_cluster_id=virtual_cluster_id
        )
        result = {
            "id": job_run.id,
            "virtualClusterId": job_run.virtual_cluster_id,
        }
        return ActionResult(result)

    def list_job_runs(self) -> ActionResult:
        virtual_cluster_id = self._get_param("virtualClusterId")
        created_before = self._get_param("createdBefore")
        created_after = self._get_param("createdAfter")
        name = self._get_param("name")
        states = self.querystring.get("states", [])
        job_runs = self.emrcontainers_backend.list_job_runs(
            virtual_cluster_id=virtual_cluster_id,
            created_before=created_before,
            created_after=created_after,
            name=name,
            states=states,
        )
        result = {"jobRuns": job_runs}
        return PaginatedResult(result)

    def describe_job_run(self) -> ActionResult:
        job_run_id = self._get_param("id")
        virtual_cluster_id = self._get_param("virtualClusterId")
        job_run = self.emrcontainers_backend.describe_job_run(
            job_run_id=job_run_id, virtual_cluster_id=virtual_cluster_id
        )
        result = {"jobRun": job_run}
        return ActionResult(result)
