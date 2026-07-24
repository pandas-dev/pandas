"""
Response handlers for AWS CloudWatch Synthetics API emulation in Moto.
"""

import json

from moto.core.responses import BaseResponse
from moto.synthetics.models import SyntheticsBackend, synthetics_backends


class SyntheticsResponse(BaseResponse):
    """
    Handles API responses for AWS CloudWatch Synthetics operations.
    """

    def __init__(self) -> None:
        """
        Initialize the SyntheticsResponse with the synthetics service name.
        """
        super().__init__(service_name="synthetics")

    @property
    def synthetics_backend(self) -> SyntheticsBackend:
        """
        Returns the backend instance for the current region.
        """
        return synthetics_backends[self.current_account][self.region]

    def create_canary(self) -> str:
        """
        Create a new canary using the provided parameters.
        """
        params = json.loads(self.body)
        canary = self.synthetics_backend.create_canary(
            name=params["Name"],
            code=params.get("Code", {}),
            artifact_s3_location=params.get("ArtifactS3Location", "s3://dummy"),
            execution_role_arn=params.get(
                "ExecutionRoleArn", "arn:aws:iam::123:role/service-role"
            ),
            schedule=params.get("Schedule", {"Expression": "rate(5 minutes)"}),
            run_config=params.get("RunConfig", {"TimeoutInSeconds": 60}),
            success_retention_period_in_days=params.get(
                "SuccessRetentionPeriodInDays", 31
            ),
            failure_retention_period_in_days=params.get(
                "FailureRetentionPeriodInDays", 31
            ),
            runtime_version=params.get("RuntimeVersion", "syn-nodejs-puppeteer-3.8"),
            vpc_config=params.get("VpcConfig"),
            resources_to_replicate_tags=params.get("ResourcesToReplicateTags"),
            provisioned_resource_cleanup=params.get("ProvisionedResourceCleanup"),
            browser_configs=params.get("BrowserConfigs"),
            tags=params.get("Tags", {}),
            artifact_config=params.get("ArtifactConfig"),
        )
        return json.dumps({"Canary": canary.to_dict()})

    def get_canary(self) -> str:
        """
        Retrieve details for a specific canary by name.
        """
        # Extract name from the URL path /canary/MyCanary
        path_parts = self.path.split("/")
        name = path_parts[-1] if len(path_parts) > 1 else None
        if not name:
            raise ValueError("Canary name not found in URL")
        canary = self.synthetics_backend.get_canary(name)
        return json.dumps({"Canary": canary.to_dict()})

    def describe_canaries(self) -> str:
        """
        List all canaries in the backend.
        """
        canaries, _ = self.synthetics_backend.describe_canaries(None, None, None)
        return json.dumps({"Canaries": [c.to_dict() for c in canaries]})

    def list_tags_for_resource(self) -> str:
        """
        List tags for a given resource ARN.
        """
        # Extract ARN from the URL path /tags/{resourceArn}
        path_parts = self.path.split("/")
        arn = path_parts[-1] if len(path_parts) > 1 else None
        if not arn:
            raise ValueError("Resource ARN not found in URL")
        tags = self.synthetics_backend.list_tags_for_resource(arn)
        return json.dumps({"Tags": tags})

    def _get_action(self) -> str:
        """
        Override to provide default action for root endpoint.
        """
        action = super()._get_action()
        if action is None and self.path == "/":
            return "GetHealthCheck"  # Default action for root endpoint
        return action or "Unknown"

    def get_health_check(self) -> str:
        """
        Handle root endpoint requests.
        """
        return "What would you like to do?"
