from typing import Any

from moto.core.exceptions import JsonRESTError


class ApiGatewayException(JsonRESTError):
    pass


class BadRequestException(ApiGatewayException):
    def __init__(self, message: str):
        super().__init__("BadRequestException", message)


class NotFoundException(ApiGatewayException):
    def __init__(self, message: str):
        super().__init__("NotFoundException", message)


class AccessDeniedException(ApiGatewayException):
    pass


class ConflictException(ApiGatewayException):
    code = 409

    def __init__(self, message: str):
        super().__init__("ConflictException", message)


class AwsProxyNotAllowed(BadRequestException):
    def __init__(self) -> None:
        super().__init__(
            "Integrations of type 'AWS_PROXY' currently only supports Lambda function and Firehose stream invocations."
        )


class CrossAccountNotAllowed(AccessDeniedException):
    def __init__(self) -> None:
        super().__init__(
            "AccessDeniedException", "Cross-account pass role is not allowed."
        )


class RoleNotSpecified(BadRequestException):
    def __init__(self) -> None:
        super().__init__("Role ARN must be specified for AWS integrations")


class IntegrationMethodNotDefined(BadRequestException):
    def __init__(self) -> None:
        super().__init__("Enumeration value for HttpMethod must be non-empty")


class InvalidOpenAPIDocumentException(BadRequestException):
    def __init__(self, cause: Any):
        super().__init__(
            f"Failed to parse the uploaded OpenAPI document due to: {cause.message}"
        )


class InvalidOpenApiDocVersionException(BadRequestException):
    def __init__(self) -> None:
        super().__init__("Only OpenAPI 3.x.x are currently supported")


class InvalidOpenApiModeException(BadRequestException):
    def __init__(self) -> None:
        super().__init__(
            'Enumeration value of OpenAPI import mode must be "overwrite" or "merge"',
        )


class InvalidResourcePathException(BadRequestException):
    def __init__(self) -> None:
        super().__init__(
            "Resource's path part only allow a-zA-Z0-9._- and curly braces at the beginning and the end and an optional plus sign before the closing brace."
        )


class InvalidHttpEndpoint(BadRequestException):
    def __init__(self) -> None:
        super().__init__("Invalid HTTP endpoint specified for URI")


class InvalidArn(BadRequestException):
    def __init__(self) -> None:
        super().__init__("Invalid ARN specified in the request")


class InvalidIntegrationArn(BadRequestException):
    def __init__(self) -> None:
        super().__init__("AWS ARN for integration must contain path or action")


class InvalidRequestInput(BadRequestException):
    def __init__(self) -> None:
        super().__init__("Invalid request input")


class NoIntegrationDefined(NotFoundException):
    def __init__(self) -> None:
        super().__init__("No integration defined for method")


class NoIntegrationResponseDefined(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid Response status code specified")


class NoMethodDefined(BadRequestException):
    def __init__(self) -> None:
        super().__init__("The REST API doesn't contain any methods")


class AuthorizerNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid Authorizer identifier specified")


class StageNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid stage identifier specified")


class ApiKeyNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid API Key identifier specified")


class UsagePlanNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid Usage Plan ID specified")


class ApiKeyAlreadyExists(ApiGatewayException):
    code = 409

    def __init__(self) -> None:
        super().__init__("ConflictException", "API Key already exists")


class InvalidDomainName(BadRequestException):
    code = 404

    def __init__(self) -> None:
        super().__init__("No Domain Name specified")


class DomainNameNotFound(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid domain name identifier specified")


class InvalidRestApiId(BadRequestException):
    code = 404

    def __init__(self) -> None:
        super().__init__("No Rest API Id specified")


class InvalidModelName(BadRequestException):
    code = 404

    def __init__(self) -> None:
        super().__init__("No Model Name specified")


class RestAPINotFound(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid Rest API Id specified")


class RequestValidatorNotFound(BadRequestException):
    code = 400

    def __init__(self) -> None:
        super().__init__("Invalid Request Validator Id specified")


class ModelNotFound(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid Model Name specified")


class ApiKeyValueMinLength(BadRequestException):
    code = 400

    def __init__(self) -> None:
        super().__init__("API Key value should be at least 20 characters")


class MethodNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid Method identifier specified")


class InvalidBasePathException(BadRequestException):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "API Gateway V1 doesn't support the slash character (/) in base path mappings. "
            "To create a multi-level base path mapping, use API Gateway V2."
        )


class DeploymentNotFoundException(NotFoundException):
    def __init__(self) -> None:
        super().__init__("Invalid Deployment identifier specified")


class InvalidRestApiIdForBasePathMappingException(BadRequestException):
    code = 400

    def __init__(self) -> None:
        super().__init__("Invalid REST API identifier specified")


class InvalidStageException(BadRequestException):
    code = 400

    def __init__(self) -> None:
        super().__init__("Invalid stage identifier specified")


class BasePathConflictException(ConflictException):
    def __init__(self) -> None:
        super().__init__("Base path already exists for this domain name")


class BasePathNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid base path mapping identifier specified")


class ResourceIdNotFoundException(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("Invalid resource identifier specified")


class VpcLinkNotFound(NotFoundException):
    code = 404

    def __init__(self) -> None:
        super().__init__("VPCLink not found")


class ValidationException(ApiGatewayException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class StageStillActive(BadRequestException):
    def __init__(self) -> None:
        super().__init__(
            "Active stages pointing to this deployment must be moved or deleted"
        )


class GatewayResponseNotFound(NotFoundException):
    def __init__(self) -> None:
        super().__init__("GatewayResponse not found")
