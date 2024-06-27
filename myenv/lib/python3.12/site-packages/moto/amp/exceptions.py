import json

from moto.core.exceptions import JsonRESTError


class AmpException(JsonRESTError):
    pass


class ResourceNotFoundException(AmpException):
    def __init__(self, message: str, resource_id: str, resource_type: str):
        super().__init__("ResourceNotFoundException", message)
        self.description = json.dumps(
            {
                "resourceId": resource_id,
                "message": self.message,
                "resourceType": resource_type,
            }
        )


class WorkspaceNotFound(ResourceNotFoundException):
    code = 404

    def __init__(self, workspace_id: str):
        super().__init__(
            "Workspace not found",
            resource_id=workspace_id,
            resource_type="AWS::APS::Workspace",
        )


class RuleGroupNamespaceNotFound(ResourceNotFoundException):
    code = 404

    def __init__(self, name: str):
        super().__init__(
            "RuleGroupNamespace not found",
            resource_id=name,
            resource_type="AWS::APS::RuleGroupNamespace",
        )
