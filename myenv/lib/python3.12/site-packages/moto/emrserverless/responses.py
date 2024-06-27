"""Handles incoming emrserverless requests, invokes methods, returns responses."""

import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import EMRServerlessBackend, emrserverless_backends

DEFAULT_MAX_RESULTS = 100
DEFAULT_NEXT_TOKEN = ""

"""
These are the available methods:
    can_paginate()
    cancel_job_run() -> DONE
    close()
    create_application() -> DONE
    delete_application() -> DONE
    get_application() -> DONE
    get_job_run() -> DONE
    get_paginator()
    get_waiter()
    list_applications() -> DONE
    list_job_runs() -> DONE
    list_tags_for_resource()
    start_application() -> DONE
    start_job_run() -> DONE
    stop_application() -> DONE
    tag_resource()
    untag_resource()
    update_application()
"""


class EMRServerlessResponse(BaseResponse):
    """Handler for EMRServerless requests and responses."""

    def __init__(self) -> None:
        super().__init__("emr-serverless")

    @property
    def emrserverless_backend(self) -> EMRServerlessBackend:
        """Return backend instance specific for this region."""
        return emrserverless_backends[self.current_account][self.region]

    def create_application(self) -> TYPE_RESPONSE:
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
        return 200, {}, json.dumps(dict(application))

    def delete_application(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")

        self.emrserverless_backend.delete_application(application_id=application_id)
        return 200, {}, ""

    def get_application(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")

        application = self.emrserverless_backend.get_application(
            application_id=application_id
        )
        response = {"application": application}
        return 200, {}, json.dumps(response)

    def list_applications(self) -> TYPE_RESPONSE:
        states = self.querystring.get("states", [])
        max_results = self._get_int_param("maxResults", DEFAULT_MAX_RESULTS)
        next_token = self._get_param("nextToken", DEFAULT_NEXT_TOKEN)

        applications, next_token = self.emrserverless_backend.list_applications(
            next_token=next_token,
            max_results=max_results,
            states=states,
        )
        response = {"applications": applications, "nextToken": next_token}
        return 200, {}, json.dumps(response)

    def start_application(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")

        self.emrserverless_backend.start_application(application_id=application_id)
        return 200, {}, ""

    def stop_application(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")

        self.emrserverless_backend.stop_application(application_id=application_id)
        return 200, {}, ""

    def update_application(self) -> TYPE_RESPONSE:
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
        return 200, {}, json.dumps(response)

    def start_job_run(self) -> TYPE_RESPONSE:
        # params = self._get_params()
        application_id = self._get_param(
            "applicationId"
        )  # params["application_id"] #.get()
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
        return (
            200,
            {},
            json.dumps(
                {
                    "applicationId": job_run.application_id,
                    "jobRunId": job_run.id,
                    "arn": job_run.arn,
                }
            ),
        )

    def get_job_run(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")
        job_run_id = self._get_param("jobRunId")
        job_run = self.emrserverless_backend.get_job_run(
            application_id=application_id, job_run_id=job_run_id
        )
        # TODO: adjust response
        return 200, {}, json.dumps({"jobRun": job_run.to_dict("get")})

    def cancel_job_run(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")
        job_run_id = self._get_param("jobRunId")
        application_id, job_run_id = self.emrserverless_backend.cancel_job_run(
            application_id=application_id,
            job_run_id=job_run_id,
        )

        return (
            200,
            {},
            json.dumps(dict(applicationId=application_id, jobRunId=job_run_id)),
        )

    def list_job_runs(self) -> TYPE_RESPONSE:
        application_id = self._get_param("applicationId")
        next_token = self._get_param("nextToken", DEFAULT_NEXT_TOKEN)
        max_results = self._get_int_param("maxResults", DEFAULT_MAX_RESULTS)
        created_at_after = self._get_param("createdAtAfter")
        created_at_before = self._get_param("createdAtBefore")
        states = self.querystring.get("states", [])
        job_runs, next_token = self.emrserverless_backend.list_job_runs(
            application_id=application_id,
            next_token=next_token,
            max_results=max_results,
            created_at_after=created_at_after,
            created_at_before=created_at_before,
            states=states,
        )
        return 200, {}, json.dumps({"jobRuns": job_runs, "nextToken": next_token})
