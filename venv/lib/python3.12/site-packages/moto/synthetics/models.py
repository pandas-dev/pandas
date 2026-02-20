"""SyntheticsBackend class with methods for supported APIs."""

import uuid
from typing import Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import iso_8601_datetime_with_milliseconds, utcnow


class Canary(BaseModel):
    """
    Represents a CloudWatch Synthetics Canary resource.
    """

    def __init__(
        self,
        name: str,
        code: dict[str, object],
        artifact_s3_location: str,
        execution_role_arn: str,
        schedule: dict[str, object],
        run_config: dict[str, object],
        success_retention_period_in_days: int,
        failure_retention_period_in_days: int,
        runtime_version: str,
        vpc_config: Optional[dict[str, object]],
        resources_to_replicate_tags: Optional[list[str]],
        provisioned_resource_cleanup: Optional[str],
        browser_configs: Optional[list[dict[str, object]]],
        tags: Optional[dict[str, str]],
        artifact_config: Optional[dict[str, object]],
    ):
        self.name = name
        self.id = str(uuid.uuid4())
        self.code = code
        self.artifact_s3_location = artifact_s3_location
        self.execution_role_arn = execution_role_arn
        self.schedule = schedule
        self.run_config = run_config
        self.success_retention_period_in_days = success_retention_period_in_days
        self.failure_retention_period_in_days = failure_retention_period_in_days
        self.runtime_version = runtime_version
        self.vpc_config = vpc_config
        self.resources_to_replicate_tags = resources_to_replicate_tags
        self.provisioned_resource_cleanup = provisioned_resource_cleanup
        self.browser_configs = browser_configs or []
        self.tags = tags or {}
        self.artifact_config = artifact_config
        self.state = "READY"

        now = utcnow()
        self.timeline = {
            "Created": now,
            "LastModified": now,
            "LastStarted": None,
            "LastStopped": None,
        }

    def to_dict(self) -> dict[str, object]:
        """
        Convert the Canary object to a dictionary representation.
        """
        return {
            "Id": self.id,
            "Name": self.name,
            "Code": self.code,
            "ArtifactS3Location": self.artifact_s3_location,
            "ExecutionRoleArn": self.execution_role_arn,
            "Schedule": self.schedule,
            "RunConfig": self.run_config,
            "SuccessRetentionPeriodInDays": self.success_retention_period_in_days,
            "FailureRetentionPeriodInDays": self.failure_retention_period_in_days,
            "RuntimeVersion": self.runtime_version,
            "VpcConfig": self.vpc_config,
            "ProvisionedResourceCleanup": self.provisioned_resource_cleanup,
            "BrowserConfigs": self.browser_configs,
            "Tags": self.tags,
            "ArtifactConfig": self.artifact_config,
            "Status": {
                "State": self.state,
                "StateReason": "Created by Moto",
                "StateReasonCode": "CREATE_COMPLETE",
            },
            "Timeline": {
                "Created": iso_8601_datetime_with_milliseconds(
                    self.timeline["Created"]
                ),
                "LastModified": iso_8601_datetime_with_milliseconds(
                    self.timeline["LastModified"]
                ),
                "LastStarted": None,
                "LastStopped": None,
            },
        }


class SyntheticsBackend(BaseBackend):
    """Implementation of Synthetics APIs."""

    def __init__(self, region_name: str, account_id: str):
        """
        Initialize the SyntheticsBackend with region and account.
        """
        super().__init__(region_name, account_id)
        self.canaries: dict[str, Canary] = {}

    def create_canary(
        self,
        name: str,
        code: dict[str, object],
        artifact_s3_location: str,
        execution_role_arn: str,
        schedule: dict[str, object],
        run_config: dict[str, object],
        success_retention_period_in_days: int,
        failure_retention_period_in_days: int,
        runtime_version: str,
        vpc_config: Optional[dict[str, object]],
        resources_to_replicate_tags: Optional[list[str]],
        provisioned_resource_cleanup: Optional[str],
        browser_configs: Optional[list[dict[str, object]]],
        tags: Optional[dict[str, str]],
        artifact_config: Optional[dict[str, object]],
    ) -> Canary:
        canary = Canary(
            name,
            code,
            artifact_s3_location,
            execution_role_arn,
            schedule,
            run_config,
            success_retention_period_in_days,
            failure_retention_period_in_days,
            runtime_version,
            vpc_config,
            resources_to_replicate_tags,
            provisioned_resource_cleanup,
            browser_configs,
            tags,
            artifact_config,
        )
        self.canaries[name] = canary
        return canary

    def get_canary(self, name: str, dry_run_id: Optional[str] = None) -> Canary:
        """
        The dry-run_id-parameter is not yet supported
        """
        # dry_run_id is unused, included for API compatibility
        return self.canaries[name]

    def describe_canaries(
        self,
        next_token: Optional[str],
        max_results: Optional[int],
        names: Optional[list[str]],
    ) -> tuple[list[Canary], None]:
        """
        Pagination is not yet supported
        """
        canaries = list(self.canaries.values())
        if names:
            canaries = [c for c in canaries if c.name in names]

        return canaries, None

    def list_tags_for_resource(self, resource_arn: str) -> dict[str, str]:
        # Simplified: assume resource_arn is actually the canary name
        canary = self.canaries.get(resource_arn)
        return canary.tags if canary else {}


# Exported backend dict for Moto Synthetics
synthetics_backends = BackendDict(SyntheticsBackend, "synthetics")
