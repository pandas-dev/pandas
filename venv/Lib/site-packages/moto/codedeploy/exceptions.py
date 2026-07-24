"""Exceptions raised by the codedeploy service."""

from moto.core.exceptions import JsonRESTError


class CodeDeployException(JsonRESTError):
    pass


class ApplicationDoesNotExistException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ApplicationDoesNotExistException", message)


class DeploymentDoesNotExistException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("DeploymentDoesNotExistException", message)


class ApplicationAlreadyExistsException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ApplicationAlreadyExistsException", message)


class ApplicationNameRequiredException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ApplicationNameRequiredException", message)


class DeploymentGroupAlreadyExistsException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("DeploymentGroupAlreadyExistsException", message)


class DeploymentGroupNameRequiredException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("DeploymentGroupNameRequiredException", message)


class DeploymentGroupDoesNotExistException(CodeDeployException):
    code = 400

    def __init__(self, message: str):
        super().__init__("DeploymentGroupDoesNotExistException", message)
