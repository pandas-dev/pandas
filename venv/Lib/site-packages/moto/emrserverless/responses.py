"""Handles incoming emrserverless requests, invokes methods, returns responses."""

from moto.core.responses import ActionResult, BaseResponse, EmptyResult, PaginatedResult

from .models import EMRServerlessBackend, emrserverless_backends


class EMRServerlessResponse(BaseResponse):
    """Handler for EMRServerless requests and responses."""

    def __init__(self) -> None:
        super().__init__("emr-serverless")
        self.automated_parameter_parsing = True

    @property
    def emrserverless_backend(self) -> EMRServerlessBackend:
        """Return backend instance specific for this region."""
        return emrserverless_backends[self.current_account][self.region]

    def create_application(self) -> ActionResult:
        name = self._get_param("name")
        release_label = self._get_param("releaseLabel")
        application_type = self._get_param("type")
        client_token = self._get_param("clientToken")
        initial_capacity = self._get_param("initialCapacity")
        maximum_capacity = self._get_param("maximumCapacity")
        tags = self._get_param("tags")
        auto_start_configuration = self._get_param("autoStartConfiguration")
        auto_stop_configuration = self._get_param("autoStopConfiguration")
        network_configuration = self._get_param("networkConfiguration")
        application = self.emrserverless_backend.create_application(
            name=name,
            release_label=release_label,
            application_type=application_type,
            client_token=client_token,
            initial_capacity=initial_capacity,
            maximum_capacity=maximum_capacity,
            tags=tags,
            auto_start_configuration=auto_start_configuration,
            auto_stop_configuration=auto_stop_configuration,
            network_configuration=network_configuration,
        )
        result = {
            "applicationId": application.id,
            "name": application.name,
            "arn": application.arn,
        }
        return ActionResult(result)

    def delete_application(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        self.emrserverless_backend.delete_application(application_id=application_id)
        return EmptyResult()

    def get_application(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        application = self.emrserverless_backend.get_application(
            application_id=application_id
        )
        response = {"application": application}
        return ActionResult(response)

    def list_applications(self) -> ActionResult:
        states = self._get_param("states", [])
        applications = self.emrserverless_backend.list_applications(states=states)
        result = {
            "applications": [
                {
                    "id": application.id,
                    "name": application.name,
                    "arn": application.arn,
                    "releaseLabel": application.release_label,
                    "type": application.type,
                    "state": application.state,
                    "stateDetails": application.state_details,
                    "createdAt": application.created_at,
                    "updatedAt": application.updated_at,
                }
                for application in applications
            ]
        }
        return PaginatedResult(result)

    def start_application(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        self.emrserverless_backend.start_application(application_id=application_id)
        return EmptyResult()

    def stop_application(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        self.emrserverless_backend.stop_application(application_id=application_id)
        return EmptyResult()

    def update_application(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        initial_capacity = self._get_param("initialCapacity")
        maximum_capacity = self._get_param("maximumCapacity")
        auto_start_configuration = self._get_param("autoStartConfiguration")
        auto_stop_configuration = self._get_param("autoStopConfiguration")
        network_configuration = self._get_param("networkConfiguration")
        application = self.emrserverless_backend.update_application(
            application_id=application_id,
            initial_capacity=initial_capacity,
            maximum_capacity=maximum_capacity,
            auto_start_configuration=auto_start_configuration,
            auto_stop_configuration=auto_stop_configuration,
            network_configuration=network_configuration,
        )
        response = {"application": application}
        return ActionResult(response)

    def start_job_run(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        client_token = self._get_param("clientToken")
        execution_role_arn = self._get_param("executionRoleArn")
        job_driver = self._get_param("jobDriver")
        configuration_overrides = self._get_param("configurationOverrides")
        tags = self._get_param("tags")
        execution_timeout_minutes = self._get_int_param("executionTimeoutMinutes")
        name = self._get_param("name")
        job_run = self.emrserverless_backend.start_job_run(
            application_id=application_id,
            client_token=client_token,
            execution_role_arn=execution_role_arn,
            job_driver=job_driver,
            configuration_overrides=configuration_overrides,
            tags=tags,
            execution_timeout_minutes=execution_timeout_minutes,
            name=name,
        )
        return ActionResult(
            {
                "applicationId": job_run.application_id,
                "jobRunId": job_run.id,
                "arn": job_run.arn,
            }
        )

    def get_job_run(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        job_run_id = self._get_param("jobRunId")
        job_run = self.emrserverless_backend.get_job_run(
            application_id=application_id, job_run_id=job_run_id
        )
        return ActionResult({"jobRun": job_run})

    def cancel_job_run(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        job_run_id = self._get_param("jobRunId")
        application_id, job_run_id = self.emrserverless_backend.cancel_job_run(
            application_id=application_id,
            job_run_id=job_run_id,
        )
        return ActionResult({"applicationId": application_id, "jobRunId": job_run_id})

    def list_job_runs(self) -> ActionResult:
        application_id = self._get_param("applicationId")
        created_at_after = self._get_param("createdAtAfter")
        created_at_before = self._get_param("createdAtBefore")
        states = self._get_param("states", [])
        job_runs = self.emrserverless_backend.list_job_runs(
            application_id=application_id,
            created_at_after=created_at_after,
            created_at_before=created_at_before,
            states=states,
        )
        result = {
            "jobRuns": [
                {
                    "applicationId": job_run.application_id,
                    "id": job_run.id,
                    "name": job_run.name,
                    "arn": job_run.arn,
                    "createdBy": job_run.created_by,
                    "createdAt": job_run.created_at,
                    "updatedAt": job_run.updated_at,
                    "executionRole": job_run.execution_role_arn,
                    "state": job_run.state,
                    "stateDetails": job_run.state_details,
                    "releaseLabel": job_run.release_label,
                    "type": job_run.application_type,
                }
                for job_run in job_runs
            ]
        }
        return PaginatedResult(result)
