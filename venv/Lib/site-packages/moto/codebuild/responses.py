import json
import re
from typing import Any, Dict, List

from moto.core.responses import BaseResponse
from moto.utilities.utils import get_partition

from .exceptions import (
    InvalidInputException,
    ResourceAlreadyExistsException,
    ResourceNotFoundException,
)
from .models import CodeBuildBackend, codebuild_backends


def _validate_required_params_source(source: Dict[str, Any]) -> None:
    if source["type"] not in [
        "BITBUCKET",
        "CODECOMMIT",
        "CODEPIPELINE",
        "GITHUB",
        "GITHUB_ENTERPRISE",
        "NO_SOURCE",
        "S3",
    ]:
        raise InvalidInputException("Invalid type provided: Project source type")

    if "location" not in source:
        raise InvalidInputException("Project source location is required")

    if source["location"] == "":
        raise InvalidInputException("Project source location is required")


def _validate_required_params_service_role(
    account_id: str, region_name: str, service_role: str
) -> None:
    if not service_role.startswith(
        f"arn:{get_partition(region_name)}:iam::{account_id}:role/"
    ):
        raise InvalidInputException(
            "Invalid service role: Service role account ID does not match caller's account"
        )


def _validate_required_params_artifacts(artifacts: Dict[str, Any]) -> None:
    if artifacts["type"] not in ["CODEPIPELINE", "S3", "NO_ARTIFACTS"]:
        raise InvalidInputException("Invalid type provided: Artifact type")

    if artifacts["type"] == "NO_ARTIFACTS":
        if "location" in artifacts:
            raise InvalidInputException(
                "Invalid artifacts: artifact type NO_ARTIFACTS should have null location"
            )
    elif "location" not in artifacts or artifacts["location"] == "":
        raise InvalidInputException("Project source location is required")


def _validate_required_params_environment(environment: Dict[str, Any]) -> None:
    if environment["type"] not in [
        "WINDOWS_CONTAINER",
        "LINUX_CONTAINER",
        "LINUX_GPU_CONTAINER",
        "ARM_CONTAINER",
    ]:
        raise InvalidInputException(f"Invalid type provided: {environment['type']}")

    if environment["computeType"] not in [
        "BUILD_GENERAL1_SMALL",
        "BUILD_GENERAL1_MEDIUM",
        "BUILD_GENERAL1_LARGE",
        "BUILD_GENERAL1_2XLARGE",
    ]:
        raise InvalidInputException(
            f"Invalid compute type provided: {environment['computeType']}"
        )


def _validate_required_params_project_name(name: str) -> None:
    if len(name) >= 150:
        raise InvalidInputException(
            "Only alphanumeric characters, dash, and underscore are supported"
        )

    if not re.match(r"^[A-Za-z]{1}.*[^!£$%^&*()+=|?`¬{}@~#:;<>\\/\[\]]$", name):
        raise InvalidInputException(
            "Only alphanumeric characters, dash, and underscore are supported"
        )


def _validate_required_params_id(build_id: str, build_ids: List[str]) -> None:
    if ":" not in build_id:
        raise InvalidInputException("Invalid build ID provided")

    if build_id not in build_ids:
        raise ResourceNotFoundException(f"Build {build_id} does not exist")


class CodeBuildResponse(BaseResponse):
    @property
    def codebuild_backend(self) -> CodeBuildBackend:
        return codebuild_backends[self.current_account][self.region]

    def list_builds_for_project(self) -> str:
        _validate_required_params_project_name(self._get_param("projectName"))

        if (
            self._get_param("projectName")
            not in self.codebuild_backend.codebuild_projects.keys()
        ):
            name = self._get_param("projectName")
            raise ResourceNotFoundException(
                f"The provided project arn:{get_partition(self.region)}:codebuild:{self.region}:{self.current_account}:project/{name} does not exist"
            )

        ids = self.codebuild_backend.list_builds_for_project(
            self._get_param("projectName")
        )

        return json.dumps({"ids": ids})

    def create_project(self) -> str:
        _validate_required_params_source(self._get_param("source"))
        service_role = self._get_param("serviceRole")
        _validate_required_params_service_role(
            self.current_account, self.region, service_role
        )
        _validate_required_params_artifacts(self._get_param("artifacts"))
        _validate_required_params_environment(self._get_param("environment"))
        _validate_required_params_project_name(self._get_param("name"))

        if self._get_param("name") in self.codebuild_backend.codebuild_projects.keys():
            name = self._get_param("name")
            raise ResourceAlreadyExistsException(
                f"Project already exists: arn:{get_partition(self.region)}:codebuild:{self.region}:{self.current_account}:project/{name}"
            )

        project_metadata = self.codebuild_backend.create_project(
            project_name=self._get_param("name"),
            description=self._get_param("description"),
            project_source=self._get_param("source"),
            artifacts=self._get_param("artifacts"),
            environment=self._get_param("environment"),
            service_role=service_role,
            tags=self._get_param("tags"),
            cache=self._get_param("cache"),
            timeout=self._get_param("timeoutInMinutes"),
            queued_timeout=self._get_param("queuedTimeoutInMinutes"),
            source_version=self._get_param("sourceVersion"),
            logs_config=self._get_param("logsConfig"),
            vpc_config=self._get_param("vpcConfig"),
        )

        return json.dumps({"project": project_metadata})

    def list_projects(self) -> str:
        project_metadata = self.codebuild_backend.list_projects()
        return json.dumps({"projects": project_metadata})

    def batch_get_projects(self) -> str:
        names = self._get_param("names")
        project_metadata = self.codebuild_backend.batch_get_projects(names)
        return json.dumps({"projects": project_metadata})

    def start_build(self) -> str:
        _validate_required_params_project_name(self._get_param("projectName"))

        if (
            self._get_param("projectName")
            not in self.codebuild_backend.codebuild_projects.keys()
        ):
            name = self._get_param("projectName")
            raise ResourceNotFoundException(
                f"Project cannot be found: arn:{get_partition(self.region)}:codebuild:{self.region}:{self.current_account}:project/{name}"
            )

        metadata = self.codebuild_backend.start_build(
            self._get_param("projectName"),
            self._get_param("sourceVersion"),
            self._get_param("artifactsOverride"),
        )
        return json.dumps({"build": metadata})

    def batch_get_builds(self) -> str:
        for build_id in self._get_param("ids"):
            if ":" not in build_id:
                raise InvalidInputException("Invalid build ID provided")

        metadata = self.codebuild_backend.batch_get_builds(self._get_param("ids"))
        return json.dumps({"builds": metadata})

    def list_builds(self) -> str:
        ids = self.codebuild_backend.list_builds()
        return json.dumps({"ids": ids})

    def delete_project(self) -> str:
        _validate_required_params_project_name(self._get_param("name"))

        self.codebuild_backend.delete_project(self._get_param("name"))
        return "{}"

    def stop_build(self) -> str:
        _validate_required_params_id(
            self._get_param("id"), self.codebuild_backend.list_builds()
        )

        metadata = self.codebuild_backend.stop_build(self._get_param("id"))
        return json.dumps({"build": metadata})
