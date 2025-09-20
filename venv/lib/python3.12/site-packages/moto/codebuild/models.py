import datetime
from collections import defaultdict
from typing import Any, Dict, List, Optional

from dateutil import parser

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.utils import get_partition


class CodeBuildProjectMetadata(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        project_name: str,
        source_version: Optional[str],
        artifacts: Optional[Dict[str, Any]],
        build_id: str,
        service_role: str,
    ):
        current_date = iso_8601_datetime_with_milliseconds()
        self.build_metadata: Dict[str, Any] = dict()

        self.build_metadata["id"] = build_id
        self.build_metadata["arn"] = (
            f"arn:{get_partition(region_name)}:codebuild:{region_name}:{account_id}:build/{build_id}"
        )

        self.build_metadata["buildNumber"] = mock_random.randint(1, 100)
        self.build_metadata["startTime"] = current_date
        self.build_metadata["currentPhase"] = "QUEUED"
        self.build_metadata["buildStatus"] = "IN_PROGRESS"
        self.build_metadata["sourceVersion"] = (
            source_version if source_version else "refs/heads/main"
        )
        self.build_metadata["projectName"] = project_name

        self.build_metadata["phases"] = [
            {
                "phaseType": "SUBMITTED",
                "phaseStatus": "SUCCEEDED",
                "startTime": current_date,
                "endTime": current_date,
                "durationInSeconds": 0,
            },
            {"phaseType": "QUEUED", "startTime": current_date},
        ]

        self.build_metadata["source"] = {
            "type": "CODECOMMIT",  # should be different based on what you pass in
            "location": "https://git-codecommit.eu-west-2.amazonaws.com/v1/repos/testing",
            "gitCloneDepth": 1,
            "gitSubmodulesConfig": {"fetchSubmodules": False},
            "buildspec": "buildspec/stuff.yaml",  # should present in the codebuild project somewhere
            "insecureSsl": False,
        }

        self.build_metadata["secondarySources"] = []
        self.build_metadata["secondarySourceVersions"] = []
        self.build_metadata["artifacts"] = artifacts
        self.build_metadata["secondaryArtifacts"] = []
        self.build_metadata["cache"] = {"type": "NO_CACHE"}

        self.build_metadata["environment"] = {
            "type": "LINUX_CONTAINER",
            "image": "aws/codebuild/amazonlinux2-x86_64-standard:3.0",
            "computeType": "BUILD_GENERAL1_SMALL",
            "environmentVariables": [],
            "privilegedMode": False,
            "imagePullCredentialsType": "CODEBUILD",
        }

        self.build_metadata["serviceRole"] = service_role

        self.build_metadata["logs"] = {
            "deepLink": "https://console.aws.amazon.com/cloudwatch/home?region=eu-west-2#logEvent:group=null;stream=null",
            "cloudWatchLogsArn": f"arn:{get_partition(region_name)}:logs:{region_name}:{account_id}:log-group:null:log-stream:null",
            "cloudWatchLogs": {"status": "ENABLED"},
            "s3Logs": {"status": "DISABLED", "encryptionDisabled": False},
        }

        self.build_metadata["timeoutInMinutes"] = 45
        self.build_metadata["queuedTimeoutInMinutes"] = 480
        self.build_metadata["buildComplete"] = False
        self.build_metadata["initiator"] = "rootme"
        self.build_metadata["encryptionKey"] = (
            f"arn:{get_partition(region_name)}:kms:{region_name}:{account_id}:alias/aws/s3"
        )


class CodeBuild(BaseModel):
    def __init__(
        self,
        account_id: str,
        region: str,
        project_name: str,
        description: Optional[str],
        project_source: Dict[str, Any],
        artifacts: Dict[str, Any],
        environment: Dict[str, Any],
        serviceRole: str = "some_role",
        tags: Optional[List[Dict[str, str]]] = None,
        cache: Optional[Dict[str, Any]] = None,
        timeout: Optional[int] = 0,
        queued_timeout: Optional[int] = 0,
        source_version: Optional[str] = None,
        logs_config: Optional[Dict[str, Any]] = None,
        vpc_config: Optional[Dict[str, Any]] = None,
    ):
        self.arn = f"arn:{get_partition(region)}:codebuild:{region}:{account_id}:project/{project_name}"
        self.service_role = serviceRole
        self.tags = tags
        current_date = unix_time()
        self.project_metadata: Dict[str, Any] = dict()

        self.project_metadata["name"] = project_name
        if description:
            self.project_metadata["description"] = description
        self.project_metadata["arn"] = self.arn
        self.project_metadata["encryptionKey"] = (
            f"arn:{get_partition(region)}:kms:{region}:{account_id}:alias/aws/s3"
        )
        if serviceRole.startswith("arn:"):
            self.project_metadata["serviceRole"] = serviceRole
        else:
            self.project_metadata["serviceRole"] = (
                f"arn:{get_partition(region)}:iam::{account_id}:role/service-role/{serviceRole}"
            )
        self.project_metadata["lastModifiedDate"] = current_date
        self.project_metadata["created"] = current_date
        self.project_metadata["badge"] = dict()
        self.project_metadata["badge"]["badgeEnabled"] = (
            False  # this false needs to be a json false not a python false
        )
        self.project_metadata["environment"] = environment
        self.project_metadata["artifacts"] = artifacts
        self.project_metadata["source"] = project_source
        self.project_metadata["cache"] = cache or {"type": "NO_CACHE"}
        self.project_metadata["timeoutInMinutes"] = timeout or 0
        self.project_metadata["queuedTimeoutInMinutes"] = queued_timeout or 0
        self.project_metadata["tags"] = tags
        if source_version:
            self.project_metadata["sourceVersion"] = source_version
        if logs_config:
            self.project_metadata["logsConfig"] = logs_config
        if vpc_config:
            self.project_metadata["vpcConfig"] = vpc_config


class CodeBuildBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.codebuild_projects: Dict[str, CodeBuild] = dict()
        self.build_history: Dict[str, List[str]] = dict()
        self.build_metadata: Dict[str, CodeBuildProjectMetadata] = dict()
        self.build_metadata_history: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    def create_project(
        self,
        project_name: str,
        description: Optional[str],
        project_source: Dict[str, Any],
        artifacts: Dict[str, Any],
        environment: Dict[str, Any],
        service_role: str,
        tags: Optional[List[Dict[str, str]]],
        cache: Optional[Dict[str, Any]],
        timeout: Optional[int],
        queued_timeout: Optional[int],
        source_version: Optional[str],
        logs_config: Optional[Dict[str, Any]],
        vpc_config: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        self.codebuild_projects[project_name] = CodeBuild(
            self.account_id,
            self.region_name,
            project_name=project_name,
            description=description,
            project_source=project_source,
            artifacts=artifacts,
            environment=environment,
            serviceRole=service_role,
            tags=tags,
            cache=cache,
            timeout=timeout,
            queued_timeout=queued_timeout,
            source_version=source_version,
            logs_config=logs_config,
            vpc_config=vpc_config,
        )

        # empty build history
        self.build_history[project_name] = list()

        return self.codebuild_projects[project_name].project_metadata

    def list_projects(self) -> List[str]:
        projects = []

        for project in self.codebuild_projects.keys():
            projects.append(project)

        return projects

    def batch_get_projects(self, names: List[str]) -> List[Dict[str, Any]]:
        result = []
        for name in names:
            if name in self.codebuild_projects:
                result.append(self.codebuild_projects[name].project_metadata)
            elif name.startswith("arn:"):
                for project in self.codebuild_projects.values():
                    if name == project.arn:
                        result.append(project.project_metadata)
        return result

    def start_build(
        self,
        project_name: str,
        source_version: Optional[str] = None,
        artifact_override: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        project = self.codebuild_projects[project_name]
        build_id = f"{project_name}:{mock_random.uuid4()}"

        # construct a new build
        self.build_metadata[project_name] = CodeBuildProjectMetadata(
            self.account_id,
            self.region_name,
            project_name,
            source_version,
            artifact_override,
            build_id,
            project.service_role,
        )

        self.build_history[project_name].append(build_id)

        # update build histroy with metadata for build id
        self.build_metadata_history[project_name].append(
            self.build_metadata[project_name].build_metadata
        )

        return self.build_metadata[project_name].build_metadata

    def _set_phases(self, phases: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        current_date = iso_8601_datetime_with_milliseconds()
        # No phaseStatus for QUEUED on first start
        for existing_phase in phases:
            if existing_phase["phaseType"] == "QUEUED":
                existing_phase["phaseStatus"] = "SUCCEEDED"

        statuses = [
            "PROVISIONING",
            "DOWNLOAD_SOURCE",
            "INSTALL",
            "PRE_BUILD",
            "BUILD",
            "POST_BUILD",
            "UPLOAD_ARTIFACTS",
            "FINALIZING",
            "COMPLETED",
        ]

        for status in statuses:
            phase: Dict[str, Any] = dict()
            phase["phaseType"] = status
            phase["phaseStatus"] = "SUCCEEDED"
            phase["startTime"] = current_date
            phase["endTime"] = current_date
            phase["durationInSeconds"] = mock_random.randint(10, 100)
            phases.append(phase)

        return phases

    def batch_get_builds(self, ids: List[str]) -> List[Dict[str, Any]]:
        batch_build_metadata: List[Dict[str, Any]] = []

        for metadata in self.build_metadata_history.values():
            for build in metadata:
                if build["id"] in ids:
                    build["phases"] = self._set_phases(build["phases"])
                    build["endTime"] = iso_8601_datetime_with_milliseconds(
                        parser.parse(build["startTime"])
                        + datetime.timedelta(minutes=mock_random.randint(1, 5))
                    )
                    build["currentPhase"] = "COMPLETED"
                    build["buildStatus"] = "SUCCEEDED"

                    batch_build_metadata.append(build)

        return batch_build_metadata

    def list_builds_for_project(self, project_name: str) -> List[str]:
        try:
            return self.build_history[project_name]
        except KeyError:
            return list()

    def list_builds(self) -> List[str]:
        ids = []

        for build_ids in self.build_history.values():
            ids += build_ids
        return ids

    def delete_project(self, project_name: str) -> None:
        self.build_metadata.pop(project_name, None)
        self.codebuild_projects.pop(project_name, None)

    def stop_build(self, build_id: str) -> Optional[Dict[str, Any]]:  # type: ignore[return]
        for metadata in self.build_metadata_history.values():
            for build in metadata:
                if build["id"] == build_id:
                    # set completion properties with variable completion time
                    build["phases"] = self._set_phases(build["phases"])
                    build["endTime"] = iso_8601_datetime_with_milliseconds(
                        parser.parse(build["startTime"])
                        + datetime.timedelta(minutes=mock_random.randint(1, 5))
                    )
                    build["currentPhase"] = "COMPLETED"
                    build["buildStatus"] = "STOPPED"

                    return build


codebuild_backends = BackendDict(CodeBuildBackend, "codebuild")
