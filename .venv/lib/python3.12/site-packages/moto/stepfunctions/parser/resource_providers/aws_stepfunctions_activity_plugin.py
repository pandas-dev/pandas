from typing import Optional

from localstack.services.cloudformation.resource_provider import (
    CloudFormationResourceProviderPlugin,
    ResourceProvider,
)


class StepFunctionsActivityProviderPlugin(CloudFormationResourceProviderPlugin):
    name = "AWS::StepFunctions::Activity"

    def __init__(self):
        self.factory: Optional[type[ResourceProvider]] = None

    def load(self):
        from localstack.services.stepfunctions.resource_providers.aws_stepfunctions_activity import (
            StepFunctionsActivityProvider,
        )

        self.factory = StepFunctionsActivityProvider
