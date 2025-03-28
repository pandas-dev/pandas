import re
import string
import time
from collections import defaultdict
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple, Union
from urllib.parse import urlparse

import requests
import responses

try:
    # Recommended as of 0.7.x
    from openapi_spec_validator import validate  # type: ignore
except ImportError:
    # Only used in < 0.7.x
    # (Also exists in 0.7.0, but throws a warning)
    from openapi_spec_validator import validate_spec as validate  # type: ignore
from openapi_spec_validator.validation.exceptions import OpenAPIValidationError

from moto.apigateway.exceptions import MethodNotFoundException
from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel, CloudFormationModel
from moto.core.utils import path_url
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from ..core.models import responses_mock
from .exceptions import (
    ApiKeyAlreadyExists,
    ApiKeyNotFoundException,
    ApiKeyValueMinLength,
    AuthorizerNotFoundException,
    AwsProxyNotAllowed,
    BadRequestException,
    BasePathConflictException,
    BasePathNotFoundException,
    ConflictException,
    CrossAccountNotAllowed,
    DeploymentNotFoundException,
    DomainNameNotFound,
    GatewayResponseNotFound,
    IntegrationMethodNotDefined,
    InvalidArn,
    InvalidBasePathException,
    InvalidDomainName,
    InvalidHttpEndpoint,
    InvalidIntegrationArn,
    InvalidModelName,
    InvalidOpenAPIDocumentException,
    InvalidOpenApiDocVersionException,
    InvalidOpenApiModeException,
    InvalidResourcePathException,
    InvalidRestApiId,
    InvalidRestApiIdForBasePathMappingException,
    InvalidStageException,
    ModelNotFound,
    NoIntegrationDefined,
    NoIntegrationResponseDefined,
    NoMethodDefined,
    RequestValidatorNotFound,
    ResourceIdNotFoundException,
    RestAPINotFound,
    RoleNotSpecified,
    StageNotFoundException,
    StageStillActive,
    UsagePlanNotFoundException,
    ValidationException,
    VpcLinkNotFound,
)
from .utils import (
    ApigwApiKeyIdentifier,
    ApigwAuthorizerIdentifier,
    ApigwDeploymentIdentifier,
    ApigwModelIdentifier,
    ApigwRequestValidatorIdentifier,
    ApigwResourceIdentifier,
    ApigwRestApiIdentifier,
    ApigwUsagePlanIdentifier,
    ApigwVpcLinkIdentifier,
    create_id,
    to_path,
)

STAGE_URL = "https://{api_id}.execute-api.{region_name}.amazonaws.com/{stage_name}"
PATCH_OPERATIONS = ["add", "remove", "replace", "move", "copy", "test"]


class Deployment(CloudFormationModel):
    def __init__(self, deployment_id: str, name: str, description: str = ""):
        self.id = deployment_id
        self.stage_name = name
        self.description = description
        self.created_date = int(time.time())

    def to_json(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "stageName": self.stage_name,
            "description": self.description,
            "createdDate": self.created_date,
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Deployment"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::ApiGateway::Deployment"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Deployment":
        properties = cloudformation_json["Properties"]
        rest_api_id = properties["RestApiId"]
        name = properties.get("StageName")
        desc = properties.get("Description", "")
        backend: "APIGatewayBackend" = apigateway_backends[account_id][region_name]
        return backend.create_deployment(
            function_id=rest_api_id, name=name, description=desc
        )


class IntegrationResponse(BaseModel):
    def __init__(
        self,
        status_code: Union[str, int],
        selection_pattern: Optional[str] = None,
        response_templates: Optional[Dict[str, Any]] = None,
        response_parameters: Optional[Dict[str, str]] = None,
        content_handling: Optional[Any] = None,
    ):
        if response_templates is None:
            # response_templates = {"application/json": None}  # Note: removed for compatibility with TF
            response_templates = {}
        for key in response_templates.keys():
            response_templates[key] = (
                response_templates[key] or None
            )  # required for compatibility with TF
        self.response_templates = response_templates
        self.status_code = status_code
        self.selection_pattern = selection_pattern
        self.response_parameters = response_parameters
        self.content_handling = content_handling

    def to_json(self) -> Dict[str, Any]:
        resp = {
            "responseTemplates": self.response_templates,
            "statusCode": self.status_code,
        }
        if self.selection_pattern:
            resp["selectionPattern"] = self.selection_pattern
        if self.content_handling:
            resp["contentHandling"] = self.content_handling
        if self.response_parameters:
            resp["responseParameters"] = self.response_parameters
        return resp


class Integration(BaseModel):
    def __init__(
        self,
        integration_type: str,
        uri: str,
        http_method: str,
        request_templates: Optional[Dict[str, Any]] = None,
        passthrough_behavior: Optional[str] = "WHEN_NO_MATCH",
        cache_key_parameters: Optional[List[str]] = None,
        tls_config: Optional[Dict[str, Any]] = None,
        cache_namespace: Optional[str] = None,
        timeout_in_millis: Optional[str] = None,
        request_parameters: Optional[Dict[str, Any]] = None,
        content_handling: Optional[str] = None,
        credentials: Optional[str] = None,
        connection_type: Optional[str] = None,
    ):
        self.integration_type = integration_type
        self.uri = uri
        self.http_method = http_method if integration_type != "MOCK" else None
        self.passthrough_behaviour = passthrough_behavior
        self.cache_key_parameters: List[str] = cache_key_parameters or []
        self.request_templates = request_templates
        self.tls_config = tls_config
        self.cache_namespace = cache_namespace
        self.timeout_in_millis = timeout_in_millis
        self.request_parameters = request_parameters
        self.content_handling = content_handling
        self.credentials = credentials
        self.connection_type = connection_type
        self.integration_responses: Optional[Dict[str, IntegrationResponse]] = None

    def to_json(self) -> Dict[str, Any]:
        int_responses: Optional[Dict[str, Any]] = None
        if self.integration_responses is not None:
            int_responses = {
                k: v.to_json() for k, v in self.integration_responses.items()
            }
        return {
            "type": self.integration_type,
            "uri": self.uri,
            "httpMethod": self.http_method,
            "passthroughBehavior": self.passthrough_behaviour,
            "cacheKeyParameters": self.cache_key_parameters,
            "requestTemplates": self.request_templates,
            "integrationResponses": int_responses,
            "tlsConfig": self.tls_config,
            "cacheNamespace": self.cache_namespace,
            "timeoutInMillis": self.timeout_in_millis,
            "requestParameters": self.request_parameters,
            "contentHandling": self.content_handling,
            "credentials": self.credentials,
            "connectionType": self.connection_type,
        }

    def create_integration_response(
        self,
        status_code: str,
        selection_pattern: str,
        response_templates: Dict[str, str],
        response_parameters: Dict[str, str],
        content_handling: str,
    ) -> IntegrationResponse:
        integration_response = IntegrationResponse(
            status_code,
            selection_pattern,
            response_templates or None,
            response_parameters,
            content_handling,
        )
        if self.integration_responses is None:
            self.integration_responses = {}
        self.integration_responses[status_code] = integration_response
        return integration_response

    def get_integration_response(self, status_code: str) -> IntegrationResponse:
        result = (self.integration_responses or {}).get(status_code)
        if not result:
            raise NoIntegrationResponseDefined()
        return result

    def delete_integration_response(self, status_code: str) -> IntegrationResponse:
        return (self.integration_responses or {}).pop(status_code, None)  # type: ignore[arg-type]


class MethodResponse(BaseModel):
    def __init__(
        self,
        status_code: str,
        response_models: Dict[str, str],
        response_parameters: Dict[str, Dict[str, str]],
    ):
        self.status_code = status_code
        self.response_models = response_models
        self.response_parameters = response_parameters

    def to_json(self) -> Dict[str, Any]:
        return {
            "statusCode": self.status_code,
            "responseModels": self.response_models,
            "responseParameters": self.response_parameters,
        }


class Method(CloudFormationModel):
    def __init__(
        self, method_type: str, authorization_type: Optional[str], **kwargs: Any
    ):
        self.http_method = method_type
        self.authorization_type = authorization_type
        self.authorizer_id = kwargs.get("authorizer_id")
        self.authorization_scopes = kwargs.get("authorization_scopes")
        self.api_key_required = kwargs.get("api_key_required") or False
        self.request_parameters = kwargs.get("request_parameters")
        self.request_models = kwargs.get("request_models")
        self.method_integration: Optional[Integration] = None
        self.operation_name = kwargs.get("operation_name")
        self.request_validator_id = kwargs.get("request_validator_id")
        self.method_responses: Dict[str, MethodResponse] = {}

    def to_json(self) -> Dict[str, Any]:
        return {
            "httpMethod": self.http_method,
            "authorizationType": self.authorization_type,
            "authorizerId": self.authorizer_id,
            "authorizationScopes": self.authorization_scopes,
            "apiKeyRequired": self.api_key_required,
            "requestParameters": self.request_parameters,
            "requestModels": self.request_models,
            "methodIntegration": self.method_integration.to_json()
            if self.method_integration
            else None,
            "operationName": self.operation_name,
            "requestValidatorId": self.request_validator_id,
            "methodResponses": {
                k: v.to_json() for k, v in self.method_responses.items()
            },
        }

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Method"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::ApiGateway::Method"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Method":
        properties = cloudformation_json["Properties"]
        rest_api_id = properties["RestApiId"]
        resource_id = properties["ResourceId"]
        method_type = properties["HttpMethod"]
        auth_type = properties["AuthorizationType"]
        key_req = properties.get("ApiKeyRequired")
        backend = apigateway_backends[account_id][region_name]
        m = backend.put_method(
            function_id=rest_api_id,
            resource_id=resource_id,
            method_type=method_type,
            authorization_type=auth_type,
            api_key_required=key_req,
        )
        int_method = properties["Integration"]["IntegrationHttpMethod"]
        int_type = properties["Integration"]["Type"]
        int_uri = properties["Integration"]["Uri"]
        backend.put_integration(
            function_id=rest_api_id,
            resource_id=resource_id,
            method_type=method_type,
            integration_type=int_type,
            uri=int_uri,
            integration_method=int_method,
        )
        return m

    def create_response(
        self,
        response_code: str,
        response_models: Dict[str, str],
        response_parameters: Dict[str, Dict[str, str]],
    ) -> MethodResponse:
        method_response = MethodResponse(
            response_code, response_models, response_parameters
        )
        self.method_responses[response_code] = method_response
        return method_response

    def get_response(self, response_code: str) -> Optional[MethodResponse]:
        return self.method_responses.get(response_code)

    def delete_response(self, response_code: str) -> Optional[MethodResponse]:
        return self.method_responses.pop(response_code, None)


class Resource(CloudFormationModel):
    def __init__(
        self,
        resource_id: str,
        account_id: str,
        region_name: str,
        api_id: str,
        path_part: str,
        parent_id: Optional[str],
    ):
        self.id = resource_id
        self.account_id = account_id
        self.region_name = region_name
        self.api_id = api_id
        self.path_part = path_part
        self.parent_id = parent_id
        self.resource_methods: Dict[str, Method] = {}
        from .integration_parsers import IntegrationParser
        from .integration_parsers.aws_parser import TypeAwsParser
        from .integration_parsers.http_parser import TypeHttpParser
        from .integration_parsers.unknown_parser import TypeUnknownParser

        self.integration_parsers: Dict[str, IntegrationParser] = defaultdict(
            TypeUnknownParser
        )
        self.integration_parsers["HTTP"] = TypeHttpParser()
        self.integration_parsers["AWS"] = TypeAwsParser()

    def to_dict(self) -> Dict[str, Any]:
        response: Dict[str, Any] = {
            "path": self.get_path(),
            "id": self.id,
        }
        if self.resource_methods:
            response["resourceMethods"] = {
                k: v.to_json() for k, v in self.resource_methods.items()
            }
        if self.parent_id:
            response["parentId"] = self.parent_id
            response["pathPart"] = self.path_part
        return response

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_name_type() -> str:
        return "Resource"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::ApiGateway::Resource"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "Resource":
        properties = cloudformation_json["Properties"]
        api_id = properties["RestApiId"]
        parent = properties["ParentId"]
        path = properties["PathPart"]

        backend = apigateway_backends[account_id][region_name]
        if parent == api_id:
            # A Root path (/) is automatically created. Any new paths should use this as their parent
            resources = backend.get_resources(function_id=api_id)
            root_id = [resource for resource in resources if resource.path_part == "/"][
                0
            ].id
            parent = root_id
        return backend.create_resource(
            function_id=api_id, parent_resource_id=parent, path_part=path
        )

    def get_path(self) -> str:
        return self.get_parent_path() + self.path_part

    def get_parent_path(self) -> str:
        if self.parent_id:
            backend = apigateway_backends[self.account_id][self.region_name]
            parent = backend.get_resource(self.api_id, self.parent_id)
            parent_path = parent.get_path()
            if parent_path != "/":  # Root parent
                parent_path += "/"
            return parent_path
        else:
            return ""

    def get_response(
        self, request: requests.PreparedRequest
    ) -> Tuple[int, Union[str, bytes]]:
        integration = self.get_integration(str(request.method))
        integration_type = integration.integration_type  # type: ignore[union-attr]

        status, result = self.integration_parsers[integration_type].invoke(
            request,
            integration,  # type: ignore[arg-type]
        )

        return status, result

    def add_method(
        self,
        method_type: str,
        authorization_type: Optional[str],
        api_key_required: Optional[bool],
        request_parameters: Any = None,
        request_models: Any = None,
        operation_name: Optional[str] = None,
        authorizer_id: Optional[str] = None,
        authorization_scopes: Any = None,
        request_validator_id: Any = None,
    ) -> Method:
        if authorization_scopes and not isinstance(authorization_scopes, list):
            authorization_scopes = [authorization_scopes]
        method = Method(
            method_type=method_type,
            authorization_type=authorization_type,
            api_key_required=api_key_required,
            request_parameters=request_parameters,
            request_models=request_models,
            operation_name=operation_name,
            authorizer_id=authorizer_id,
            authorization_scopes=authorization_scopes,
            request_validator_id=request_validator_id,
        )
        self.resource_methods[method_type] = method
        return method

    def get_method(self, method_type: str) -> Method:
        method = self.resource_methods.get(method_type)
        if not method:
            raise MethodNotFoundException()
        return method

    def delete_method(self, method_type: str) -> None:
        self.resource_methods.pop(method_type, None)

    def add_integration(
        self,
        method_type: str,
        integration_type: str,
        uri: str,
        request_templates: Optional[Dict[str, Any]] = None,
        passthrough_behavior: Optional[str] = None,
        integration_method: Optional[str] = None,
        tls_config: Optional[Dict[str, Any]] = None,
        cache_namespace: Optional[str] = None,
        timeout_in_millis: Optional[str] = None,
        request_parameters: Optional[Dict[str, Any]] = None,
        content_handling: Optional[str] = None,
        credentials: Optional[str] = None,
        connection_type: Optional[str] = None,
    ) -> Integration:
        integration_method = integration_method or method_type
        integration = Integration(
            integration_type,
            uri,
            integration_method,
            request_templates=request_templates,
            passthrough_behavior=passthrough_behavior,
            tls_config=tls_config,
            cache_namespace=cache_namespace,
            timeout_in_millis=timeout_in_millis,
            request_parameters=request_parameters,
            content_handling=content_handling,
            credentials=credentials,
            connection_type=connection_type,
        )
        self.resource_methods[method_type].method_integration = integration
        return integration

    def get_integration(self, method_type: str) -> Optional[Integration]:
        method = self.resource_methods.get(method_type)
        return method.method_integration if method else None

    def delete_integration(self, method_type: str) -> Integration:
        integration = self.resource_methods[method_type].method_integration
        self.resource_methods[method_type].method_integration = None
        return integration  # type: ignore[return-value]


class Authorizer(BaseModel):
    def __init__(
        self,
        authorizer_id: Optional[str],
        name: Optional[str],
        authorizer_type: Optional[str],
        **kwargs: Any,
    ):
        self.id = authorizer_id
        self.name = name
        self.type = authorizer_type
        self.provider_arns = kwargs.get("provider_arns")
        self.auth_type = kwargs.get("auth_type")
        self.authorizer_uri = kwargs.get("authorizer_uri")
        self.authorizer_credentials = kwargs.get("authorizer_credentials")
        self.identity_source = kwargs.get("identity_source")
        self.identity_validation_expression = kwargs.get(
            "identity_validation_expression"
        )
        self.authorizer_result_ttl = kwargs.get("authorizer_result_ttl")

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "id": self.id,
            "name": self.name,
            "type": self.type,
            "authorizerResultTtlInSeconds": self.authorizer_result_ttl,
        }
        if self.provider_arns:
            dct["providerARNs"] = self.provider_arns
        if self.auth_type:
            dct["authType"] = self.auth_type
        if self.authorizer_uri:
            dct["authorizerUri"] = self.authorizer_uri
        if self.authorizer_credentials:
            dct["authorizerCredentials"] = self.authorizer_credentials
        if self.identity_source:
            dct["identitySource"] = self.identity_source
        if self.identity_validation_expression:
            dct["identityValidationExpression"] = self.identity_validation_expression
        return dct

    def apply_operations(self, patch_operations: List[Dict[str, Any]]) -> "Authorizer":
        for op in patch_operations:
            if "/authorizerUri" in op["path"]:
                self.authorizer_uri = op["value"]
            elif "/authorizerCredentials" in op["path"]:
                self.authorizer_credentials = op["value"]
            elif "/authorizerResultTtlInSeconds" in op["path"]:
                self.authorizer_result_ttl = int(op["value"])
            elif "/authType" in op["path"]:
                self.auth_type = op["value"]
            elif "/identitySource" in op["path"]:
                self.identity_source = op["value"]
            elif "/identityValidationExpression" in op["path"]:
                self.identity_validation_expression = op["value"]
            elif "/name" in op["path"]:
                self.name = op["value"]
            elif "/providerARNs" in op["path"]:
                # TODO: add and remove
                raise Exception(f'Patch operation for "{op["path"]}" not implemented')
            elif "/type" in op["path"]:
                self.type = op["value"]
            else:
                raise BadRequestException(
                    f'Patch operation "{op["op"]}" not implemented'
                )
        return self


class Stage(BaseModel):
    def __init__(
        self,
        name: Optional[str] = None,
        deployment_id: Optional[str] = None,
        variables: Optional[Dict[str, Any]] = None,
        description: str = "",
        cacheClusterEnabled: Optional[bool] = False,
        cacheClusterSize: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        tracing_enabled: Optional[bool] = None,
    ):
        self.name = name
        self.deployment_id = deployment_id
        self.method_settings: Dict[str, Any] = {}
        self.variables = variables or {}
        self.description = description
        self.cache_cluster_enabled = cacheClusterEnabled
        self.cache_cluster_status = "AVAILABLE" if cacheClusterEnabled else None
        self.cache_cluster_size = (
            str(cacheClusterSize) if cacheClusterSize is not None else None
        )
        self.tags = tags
        self.tracing_enabled = tracing_enabled
        self.access_log_settings: Optional[Dict[str, Any]] = None
        self.web_acl_arn: Optional[str] = None

    def to_json(self) -> Dict[str, Any]:
        dct: Dict[str, Any] = {
            "stageName": self.name,
            "deploymentId": self.deployment_id,
            "methodSettings": self.method_settings,
            "variables": self.variables,
            "description": self.description,
            "cacheClusterEnabled": self.cache_cluster_enabled,
            "accessLogSettings": self.access_log_settings,
        }
        if self.cache_cluster_status is not None:
            dct["cacheClusterStatus"] = self.cache_cluster_status
        if self.cache_cluster_enabled:
            if self.cache_cluster_size is not None:
                dct["cacheClusterSize"] = self.cache_cluster_size
            else:
                dct["cacheClusterSize"] = "0.5"
        if self.tags:
            dct["tags"] = self.tags
        if self.tracing_enabled is not None:
            dct["tracingEnabled"] = self.tracing_enabled
        if self.web_acl_arn is not None:
            dct["webAclArn"] = self.web_acl_arn
        return dct

    def apply_operations(self, patch_operations: List[Dict[str, Any]]) -> "Stage":
        for op in patch_operations:
            if "variables/" in op["path"]:
                self._apply_operation_to_variables(op)
            elif "/cacheClusterEnabled" in op["path"]:
                self.cache_cluster_enabled = self._str2bool(op["value"])
                if self.cache_cluster_enabled:
                    self.cache_cluster_status = "AVAILABLE"
                else:
                    self.cache_cluster_status = "NOT_AVAILABLE"
            elif "/cacheClusterSize" in op["path"]:
                self.cache_cluster_size = str(op["value"])
            elif "/description" in op["path"]:
                self.description = op["value"]
            elif "/deploymentId" in op["path"]:
                self.deployment_id = op["value"]
            elif op["op"] == "replace":
                if op["path"] == "/tracingEnabled":
                    self.tracing_enabled = self._str2bool(op["value"])
                elif op["path"].startswith("/accessLogSettings/"):
                    self.access_log_settings = self.access_log_settings or {}
                    self.access_log_settings[op["path"].split("/")[-1]] = op["value"]
                else:
                    # (e.g., path could be '/*/*/logging/loglevel')
                    split_path = op["path"].split("/", 3)
                    if len(split_path) != 4:
                        continue
                    self._patch_method_setting(
                        "/".join(split_path[1:3]), split_path[3], op["value"]
                    )
            elif op["op"] == "remove":
                if op["path"] == "/accessLogSettings":
                    self.access_log_settings = None
            else:
                raise ValidationException(
                    "Member must satisfy enum value set: [add, remove, move, test, replace, copy]"
                )
        return self

    def _patch_method_setting(
        self, resource_path_and_method: str, key: str, value: str
    ) -> None:
        updated_key = self._method_settings_translations(key)
        if updated_key is not None:
            if resource_path_and_method not in self.method_settings:
                self.method_settings[resource_path_and_method] = (
                    self._get_default_method_settings()
                )
            self.method_settings[resource_path_and_method][updated_key] = (
                self._convert_to_type(updated_key, value)
            )

    def _get_default_method_settings(self) -> Dict[str, Any]:
        return {
            "throttlingRateLimit": 1000.0,
            "dataTraceEnabled": False,
            "metricsEnabled": False,
            "unauthorizedCacheControlHeaderStrategy": "SUCCEED_WITH_RESPONSE_HEADER",
            "cacheTtlInSeconds": 300,
            "cacheDataEncrypted": True,
            "cachingEnabled": False,
            "throttlingBurstLimit": 2000,
            "requireAuthorizationForCacheControl": True,
        }

    def _method_settings_translations(self, key: str) -> Optional[str]:
        mappings = {
            "metrics/enabled": "metricsEnabled",
            "logging/loglevel": "loggingLevel",
            "logging/dataTrace": "dataTraceEnabled",
            "throttling/burstLimit": "throttlingBurstLimit",
            "throttling/rateLimit": "throttlingRateLimit",
            "caching/enabled": "cachingEnabled",
            "caching/ttlInSeconds": "cacheTtlInSeconds",
            "caching/dataEncrypted": "cacheDataEncrypted",
            "caching/requireAuthorizationForCacheControl": "requireAuthorizationForCacheControl",
            "caching/unauthorizedCacheControlHeaderStrategy": "unauthorizedCacheControlHeaderStrategy",
        }

        return mappings.get(key)

    def _str2bool(self, v: str) -> bool:
        return v.lower() == "true"

    def _convert_to_type(self, key: str, val: str) -> Union[str, int, float]:
        type_mappings = {
            "metricsEnabled": "bool",
            "loggingLevel": "str",
            "dataTraceEnabled": "bool",
            "throttlingBurstLimit": "int",
            "throttlingRateLimit": "float",
            "cachingEnabled": "bool",
            "cacheTtlInSeconds": "int",
            "cacheDataEncrypted": "bool",
            "requireAuthorizationForCacheControl": "bool",
            "unauthorizedCacheControlHeaderStrategy": "str",
        }

        if key in type_mappings:
            type_value = type_mappings[key]

            if type_value == "bool":
                return self._str2bool(val)
            elif type_value == "int":
                return int(val)
            elif type_value == "float":
                return float(val)
            else:
                return str(val)
        else:
            return str(val)

    def _apply_operation_to_variables(self, op: Dict[str, Any]) -> None:
        key = op["path"][op["path"].rindex("variables/") + 10 :]
        if op["op"] == "remove":
            self.variables.pop(key, None)
        elif op["op"] == "replace":
            self.variables[key] = op["value"]
        else:
            raise Exception(f'Patch operation "{op["op"]}" not implemented')


class ApiKey(BaseModel):
    def __init__(
        self,
        api_key_id: str,
        name: Optional[str] = None,
        description: Optional[str] = None,
        enabled: bool = False,
        generateDistinctId: bool = False,  # pylint: disable=unused-argument
        value: Optional[str] = None,
        stageKeys: Optional[Any] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        customerId: Optional[str] = None,
    ):
        self.id = api_key_id
        self.value = value or "".join(
            random.sample(string.ascii_letters + string.digits, 40)
        )
        self.name = name
        self.customer_id = customerId
        self.description = description
        self.enabled = enabled
        self.created_date = self.last_updated_date = int(time.time())
        self.stage_keys = stageKeys or []
        self.tags = tags

    def to_json(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "value": self.value,
            "name": self.name,
            "customerId": self.customer_id,
            "description": self.description,
            "enabled": self.enabled,
            "createdDate": self.created_date,
            "lastUpdatedDate": self.last_updated_date,
            "stageKeys": self.stage_keys,
            "tags": self.tags,
        }

    def update_operations(self, patch_operations: List[Dict[str, Any]]) -> "ApiKey":
        for op in patch_operations:
            if op["op"] == "replace":
                if "/name" in op["path"]:
                    self.name = op["value"]
                elif "/customerId" in op["path"]:
                    self.customer_id = op["value"]
                elif "/description" in op["path"]:
                    self.description = op["value"]
                elif "/enabled" in op["path"]:
                    self.enabled = self._str2bool(op["value"])
            else:
                raise Exception(f'Patch operation "{op["op"]}" not implemented')
        return self

    def _str2bool(self, v: str) -> bool:
        return v.lower() == "true"


class UsagePlan(BaseModel):
    def __init__(
        self,
        usage_plan_id: str,
        name: Optional[str] = None,
        description: Optional[str] = None,
        apiStages: Any = None,
        throttle: Optional[Dict[str, Any]] = None,
        quota: Optional[Dict[str, Any]] = None,
        productCode: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
    ):
        self.id = usage_plan_id
        self.name = name
        self.description = description
        self.api_stages = apiStages or []
        self.throttle = throttle or {}
        self.quota = quota or {}
        self.product_code = productCode
        self.tags = tags

    def to_json(self) -> Dict[str, Any]:
        resp = {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "apiStages": self.api_stages,
            "productCode": self.product_code,
            "tags": self.tags,
        }
        if self.throttle:
            resp["throttle"] = self.throttle
        if self.quota:
            resp["quota"] = self.quota
        return resp

    def apply_patch_operations(self, patch_operations: List[Dict[str, Any]]) -> None:
        for op in patch_operations:
            path = op["path"]
            if op["op"] == "add":
                value = op["value"]
                if path == "/apiStages":
                    self.api_stages.append(
                        {"apiId": value.split(":")[0], "stage": value.split(":")[1]}
                    )
            if op["op"] == "replace":
                value = op["value"]
                if "/name" in path:
                    self.name = value
                if "/description" in path:
                    self.description = value
            if op["op"] in ["add", "replace"]:
                value = op["value"]
                if "/productCode" in path:
                    self.product_code = value
                if "/quota/limit" in path:
                    self.quota["limit"] = value
                if "/quota/period" in path:
                    self.quota["period"] = value
                if path == "/quota/offset":
                    self.quota["offset"] = value
                if "/throttle/rateLimit" in path:
                    self.throttle["rateLimit"] = int(value)
                if "/throttle/burstLimit" in path:
                    self.throttle["burstLimit"] = int(value)
            if op["op"] == "remove":
                if path == "/apiStages":
                    value = op["value"]
                    self.api_stages.remove(
                        {"apiId": value.split(":")[0], "stage": value.split(":")[1]}
                    )
                if path == "/productCode":
                    self.product_code = None
                if path == "/quota":
                    self.quota.clear()
                if path == "/throttle":
                    self.throttle.clear()


class RequestValidator(BaseModel):
    PROP_ID = "id"
    PROP_NAME = "name"
    PROP_VALIDATE_REQUEST_BODY = "validateRequestBody"
    PROP_VALIDATE_REQUEST_PARAMETERS = "validateRequestParameters"

    # operations
    OP_PATH = "path"
    OP_VALUE = "value"
    OP_REPLACE = "replace"
    OP_OP = "op"

    def __init__(
        self,
        _id: str,
        name: str,
        validateRequestBody: Optional[bool],
        validateRequestParameters: Any,
    ):
        self.id = _id
        self.name = name
        self.validate_request_body = validateRequestBody
        self.validate_request_parameters = validateRequestParameters

    def apply_patch_operations(self, operations: List[Dict[str, Any]]) -> None:
        for operation in operations:
            path = operation[RequestValidator.OP_PATH]
            value = operation[RequestValidator.OP_VALUE]
            if operation[RequestValidator.OP_OP] == RequestValidator.OP_REPLACE:
                if to_path(RequestValidator.PROP_NAME) in path:
                    self.name = value
                if to_path(RequestValidator.PROP_VALIDATE_REQUEST_BODY) in path:
                    self.validate_request_body = value.lower() in ("true")
                if to_path(RequestValidator.PROP_VALIDATE_REQUEST_PARAMETERS) in path:
                    self.validate_request_parameters = value.lower() in ("true")

    def to_dict(self) -> Dict[str, Any]:
        return {
            RequestValidator.PROP_ID: self.id,
            RequestValidator.PROP_NAME: self.name,
            RequestValidator.PROP_VALIDATE_REQUEST_BODY: self.validate_request_body,
            RequestValidator.PROP_VALIDATE_REQUEST_PARAMETERS: self.validate_request_parameters,
        }


class UsagePlanKey(BaseModel):
    def __init__(self, plan_id: str, plan_type: str, name: Optional[str], value: str):
        self.id = plan_id
        self.name = name
        self.type = plan_type
        self.value = value

    def to_json(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "type": self.type,
            "value": self.value,
        }


class VpcLink(BaseModel):
    def __init__(
        self,
        vpc_link_id: str,
        name: str,
        description: str,
        target_arns: List[str],
        tags: List[Dict[str, str]],
    ):
        self.id = vpc_link_id
        self.name = name
        self.description = description
        self.target_arns = target_arns
        self.tags = tags
        self.status = "AVAILABLE"

    def to_json(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "targetArns": self.target_arns,
            "tags": self.tags,
            "status": self.status,
        }


class RestAPI(CloudFormationModel):
    PROP_ID = "id"
    PROP_NAME = "name"
    PROP_DESCRIPTION = "description"
    PROP_VERSION = "version"
    PROP_BINARY_MEDIA_TYPES = "binaryMediaTypes"
    PROP_CREATED_DATE = "createdDate"
    PROP_API_KEY_SOURCE = "apiKeySource"
    PROP_ENDPOINT_CONFIGURATION = "endpointConfiguration"
    PROP_TAGS = "tags"
    PROP_POLICY = "policy"
    PROP_DISABLE_EXECUTE_API_ENDPOINT = "disableExecuteApiEndpoint"
    PROP_MINIMUM_COMPRESSION_SIZE = "minimumCompressionSize"
    PROP_ROOT_RESOURCE_ID = "rootResourceId"

    # operations
    OPERATION_ADD = "add"
    OPERATION_REPLACE = "replace"
    OPERATION_REMOVE = "remove"
    OPERATION_PATH = "path"
    OPERATION_VALUE = "value"
    OPERATION_OP = "op"

    def __init__(
        self,
        api_id: str,
        account_id: str,
        region_name: str,
        name: str,
        description: str,
        **kwargs: Any,
    ):
        self.id = api_id
        self.account_id = account_id
        self.region_name = region_name
        self.name = name
        self.description = description
        self.version = kwargs.get(RestAPI.PROP_VERSION) or "V1"
        self.binaryMediaTypes = kwargs.get(RestAPI.PROP_BINARY_MEDIA_TYPES) or []
        self.create_date = int(time.time())
        self.api_key_source = kwargs.get("api_key_source") or "HEADER"
        self.policy = kwargs.get(RestAPI.PROP_POLICY) or None
        self.endpoint_configuration = kwargs.get("endpoint_configuration") or {
            "types": ["EDGE"]
        }
        self.tags = kwargs.get(RestAPI.PROP_TAGS) or {}
        self.disableExecuteApiEndpoint = (
            kwargs.get("disable_execute_api_endpoint") or False
        )
        self.minimum_compression_size = kwargs.get("minimum_compression_size")
        self.deployments: Dict[str, Deployment] = {}
        self.authorizers: Dict[str, Authorizer] = {}
        self.gateway_responses: Dict[str, GatewayResponse] = {}
        self.stages: Dict[str, Stage] = {}
        self.resources: Dict[str, Resource] = {}
        self.models: Dict[str, Model] = {}
        self.request_validators: Dict[str, RequestValidator] = {}
        self.default = self.add_child("/")  # Add default child
        self.root_resource_id = self.default.id

    def __repr__(self) -> str:
        return str(self.id)

    def to_dict(self) -> Dict[str, Any]:
        return {
            self.PROP_ID: self.id,
            self.PROP_NAME: self.name,
            self.PROP_DESCRIPTION: self.description,
            self.PROP_VERSION: self.version,
            self.PROP_BINARY_MEDIA_TYPES: self.binaryMediaTypes,
            self.PROP_CREATED_DATE: self.create_date,
            self.PROP_API_KEY_SOURCE: self.api_key_source,
            self.PROP_ENDPOINT_CONFIGURATION: self.endpoint_configuration,
            self.PROP_TAGS: self.tags,
            self.PROP_POLICY: self.policy,
            self.PROP_DISABLE_EXECUTE_API_ENDPOINT: self.disableExecuteApiEndpoint,
            self.PROP_MINIMUM_COMPRESSION_SIZE: self.minimum_compression_size,
            self.PROP_ROOT_RESOURCE_ID: self.root_resource_id,
        }

    def apply_patch_operations(self, patch_operations: List[Dict[str, Any]]) -> None:
        for op in patch_operations:
            path = op[self.OPERATION_PATH]
            value = ""
            if self.OPERATION_VALUE in op:
                value = op[self.OPERATION_VALUE]
            operaton = op[self.OPERATION_OP]
            if operaton == self.OPERATION_REPLACE:
                if to_path(self.PROP_NAME) in path:
                    self.name = value
                if to_path(self.PROP_DESCRIPTION) in path:
                    self.description = value
                if to_path(self.PROP_API_KEY_SOURCE) in path:
                    self.api_key_source = value
                if to_path(self.PROP_BINARY_MEDIA_TYPES) in path:
                    self.binaryMediaTypes = [value]
                if to_path(self.PROP_DISABLE_EXECUTE_API_ENDPOINT) in path:
                    self.disableExecuteApiEndpoint = bool(value)
                if to_path(self.PROP_POLICY) in path:
                    self.policy = value
            elif operaton == self.OPERATION_ADD:
                if to_path(self.PROP_BINARY_MEDIA_TYPES) in path:
                    self.binaryMediaTypes.append(value)
            elif operaton == self.OPERATION_REMOVE:
                if to_path(self.PROP_BINARY_MEDIA_TYPES) in path:
                    self.binaryMediaTypes.remove(value)
                if to_path(self.PROP_DESCRIPTION) in path:
                    self.description = ""

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:
        return attr in ["RootResourceId"]

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        from moto.cloudformation.exceptions import UnformattedGetAttTemplateException

        if attribute_name == "RootResourceId":
            for res_id, res_obj in self.resources.items():
                if res_obj.path_part == "/" and not res_obj.parent_id:
                    return res_id
            raise Exception(f"Unable to find root resource for API {self}")
        raise UnformattedGetAttTemplateException()

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_name_type() -> str:
        return "RestApi"

    @staticmethod
    def cloudformation_type() -> str:
        return "AWS::ApiGateway::RestApi"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "RestAPI":
        properties = cloudformation_json["Properties"]
        name = properties["Name"]
        desc = properties.get("Description", "")
        config = properties.get("EndpointConfiguration", None)
        backend = apigateway_backends[account_id][region_name]
        return backend.create_rest_api(
            name=name, description=desc, endpoint_configuration=config
        )

    def add_child(self, path: str, parent_id: Optional[str] = None) -> Resource:
        child_id = ApigwResourceIdentifier(
            self.account_id, self.region_name, parent_id or "", path
        ).generate()
        child = Resource(
            resource_id=child_id,
            account_id=self.account_id,
            region_name=self.region_name,
            api_id=self.id,
            path_part=path,
            parent_id=parent_id,
        )
        self.resources[child_id] = child
        return child

    def add_model(
        self,
        name: str,
        description: str,
        schema: str,
        content_type: str,
    ) -> "Model":
        model_id = ApigwModelIdentifier(
            self.account_id, self.region_name, name
        ).generate()
        new_model = Model(
            model_id=model_id,
            name=name,
            description=description,
            schema=schema,
            content_type=content_type,
        )

        self.models[name] = new_model
        return new_model

    def get_resource_for_path(self, path_after_stage_name: str) -> Resource:  # type: ignore[return]
        for resource in self.resources.values():
            if resource.get_path() == path_after_stage_name:
                return resource
        # TODO deal with no matching resource

    def resource_callback(
        self, request: Any
    ) -> Tuple[int, Dict[str, str], Union[str, bytes]]:
        path = path_url(request.url)
        path_after_stage_name = "/" + "/".join(path.split("/")[2:])

        resource = self.get_resource_for_path(path_after_stage_name)
        status_code, response = resource.get_response(request)
        return status_code, {}, response

    def update_integration_mocks(self, stage_name: str) -> None:
        stage_url_lower = STAGE_URL.format(
            api_id=self.id.lower(), region_name=self.region_name, stage_name=stage_name
        )
        stage_url_upper = STAGE_URL.format(
            api_id=self.id.upper(), region_name=self.region_name, stage_name=stage_name
        )

        for resource in self.resources.values():
            path = resource.get_path()
            path = "" if path == "/" else path

            for http_method in resource.resource_methods.keys():
                for url in [stage_url_lower, stage_url_upper]:
                    callback_response = responses.CallbackResponse(
                        url=url + path,
                        method=http_method,
                        callback=self.resource_callback,
                        content_type="text/plain",
                    )
                    responses_mock.add(callback_response)

    def create_authorizer(
        self,
        authorizer_id: str,
        name: str,
        authorizer_type: str,
        provider_arns: Optional[List[str]],
        auth_type: Optional[str],
        authorizer_uri: Optional[str],
        authorizer_credentials: Optional[str],
        identity_source: Optional[str],
        identiy_validation_expression: Optional[str],
        authorizer_result_ttl: Optional[int],
    ) -> Authorizer:
        authorizer = Authorizer(
            authorizer_id=authorizer_id,
            name=name,
            authorizer_type=authorizer_type,
            provider_arns=provider_arns,
            auth_type=auth_type,
            authorizer_uri=authorizer_uri,
            authorizer_credentials=authorizer_credentials,
            identity_source=identity_source,
            identiy_validation_expression=identiy_validation_expression,
            authorizer_result_ttl=authorizer_result_ttl,
        )
        self.authorizers[authorizer_id] = authorizer
        return authorizer

    def create_stage(
        self,
        name: str,
        deployment_id: str,
        variables: Any,
        description: str,
        cacheClusterEnabled: Optional[bool],
        cacheClusterSize: Optional[str],
        tags: Optional[Dict[str, str]],
        tracing_enabled: Optional[bool],
    ) -> Stage:
        if name in self.stages:
            raise ConflictException("Stage already exists")
        if variables is None:
            variables = {}
        stage = Stage(
            name=name,
            deployment_id=deployment_id,
            variables=variables,
            description=description,
            cacheClusterSize=cacheClusterSize,
            cacheClusterEnabled=cacheClusterEnabled,
            tags=tags,
            tracing_enabled=tracing_enabled,
        )
        self.stages[name] = stage
        self.update_integration_mocks(name)
        return stage

    def create_deployment(
        self, name: str, description: str, stage_variables: Any = None
    ) -> Deployment:
        if stage_variables is None:
            stage_variables = {}
        # Since there are no unique values to a deployment, we will use the stage name for the deployment.
        # We are also passing a list of deployment ids to generate to prevent overwriting deployments.
        deployment_id = ApigwDeploymentIdentifier(
            self.account_id, self.region_name, stage_name=name
        ).generate(list(self.deployments.keys()))
        deployment = Deployment(deployment_id, name, description)
        self.deployments[deployment_id] = deployment
        if name:
            self.stages[name] = Stage(
                name=name, deployment_id=deployment_id, variables=stage_variables
            )
        self.update_integration_mocks(name)

        return deployment

    def get_deployment(self, deployment_id: str) -> Deployment:
        return self.deployments[deployment_id]

    def get_authorizers(self) -> List[Authorizer]:
        return list(self.authorizers.values())

    def get_stages(self) -> List[Stage]:
        return list(self.stages.values())

    def get_deployments(self) -> List[Deployment]:
        return list(self.deployments.values())

    def delete_deployment(self, deployment_id: str) -> Deployment:
        if deployment_id not in self.deployments:
            raise DeploymentNotFoundException()
        deployment = self.deployments[deployment_id]
        if deployment.stage_name and deployment.stage_name in self.stages:
            # Stage is still active
            raise StageStillActive()

        return self.deployments.pop(deployment_id)

    def create_request_validator(
        self,
        name: str,
        validateRequestBody: Optional[bool],
        validateRequestParameters: Any,
    ) -> RequestValidator:
        validator_id = ApigwRequestValidatorIdentifier(
            self.account_id, self.region_name, name
        ).generate()
        request_validator = RequestValidator(
            _id=validator_id,
            name=name,
            validateRequestBody=validateRequestBody,
            validateRequestParameters=validateRequestParameters,
        )
        self.request_validators[validator_id] = request_validator
        return request_validator

    def get_request_validators(self) -> List[RequestValidator]:
        return list(self.request_validators.values())

    def get_request_validator(self, validator_id: str) -> RequestValidator:
        reqeust_validator = self.request_validators.get(validator_id)
        if reqeust_validator is None:
            raise RequestValidatorNotFound()
        return reqeust_validator

    def delete_request_validator(self, validator_id: str) -> RequestValidator:
        return self.request_validators.pop(validator_id)

    def update_request_validator(
        self, validator_id: str, patch_operations: List[Dict[str, Any]]
    ) -> RequestValidator:
        self.request_validators[validator_id].apply_patch_operations(patch_operations)
        return self.request_validators[validator_id]

    def put_gateway_response(
        self,
        response_type: str,
        status_code: int,
        response_parameters: Dict[str, Any],
        response_templates: Dict[str, str],
    ) -> "GatewayResponse":
        response = GatewayResponse(
            response_type=response_type,
            status_code=status_code,
            response_parameters=response_parameters,
            response_templates=response_templates,
        )
        self.gateway_responses[response_type] = response
        return response

    def get_gateway_response(self, response_type: str) -> "GatewayResponse":
        if response_type not in self.gateway_responses:
            raise GatewayResponseNotFound()
        return self.gateway_responses[response_type]

    def get_gateway_responses(self) -> List["GatewayResponse"]:
        return list(self.gateway_responses.values())

    def delete_gateway_response(self, response_type: str) -> None:
        self.gateway_responses.pop(response_type, None)


class DomainName(BaseModel):
    def __init__(self, domain_name: str, **kwargs: Any):
        self.domain_name = domain_name
        region = kwargs.get("region_name") or "us-east-1"
        self.regional_domain_name = (
            f"d-{create_id()}.execute-api.{region}.amazonaws.com"
        )
        self.distribution_domain_name = f"d{create_id()}.cloudfront.net"
        self.domain_name_status = "AVAILABLE"
        self.status_message = "Domain Name Available"
        self.regional_hosted_zone_id = "Z2FDTNDATAQYW2"
        self.distribution_hosted_zone_id = "Z2FDTNDATAQYW2"
        self.certificate_upload_date = int(time.time())
        self.certificate_name = kwargs.get("certificate_name")
        self.certificate_arn = kwargs.get("certificate_arn")
        self.certificate_body = kwargs.get("certificate_body")
        self.tags = kwargs.get("tags")
        self.security_policy = kwargs.get("security_policy")
        self.certificate_chain = kwargs.get("certificate_chain")
        self.regional_certificate_name = kwargs.get("regional_certificate_name")
        self.certificate_private_key = kwargs.get("certificate_private_key")
        self.regional_certificate_arn = kwargs.get("regional_certificate_arn")
        self.endpoint_configuration = kwargs.get("endpoint_configuration")

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "domainName": self.domain_name,
            "regionalDomainName": self.regional_domain_name,
            "distributionDomainName": self.distribution_domain_name,
            "domainNameStatus": self.domain_name_status,
            "domainNameStatusMessage": self.status_message,
            "regionalHostedZoneId": self.regional_hosted_zone_id,
            "distributionHostedZoneId": self.distribution_hosted_zone_id,
            "certificateUploadDate": self.certificate_upload_date,
        }
        if self.certificate_name:
            dct["certificateName"] = self.certificate_name
        if self.certificate_arn:
            dct["certificateArn"] = self.certificate_arn
        if self.certificate_body:
            dct["certificateBody"] = self.certificate_body
        if self.tags:
            dct["tags"] = self.tags
        if self.security_policy:
            dct["securityPolicy"] = self.security_policy
        if self.certificate_chain:
            dct["certificateChain"] = self.certificate_chain
        if self.regional_certificate_name:
            dct["regionalCertificateName"] = self.regional_certificate_name
        if self.certificate_private_key:
            dct["certificatePrivateKey"] = self.certificate_private_key
        if self.regional_certificate_arn:
            dct["regionalCertificateArn"] = self.regional_certificate_arn
        if self.endpoint_configuration:
            dct["endpointConfiguration"] = self.endpoint_configuration
        return dct


class Model(BaseModel):
    def __init__(self, model_id: str, name: str, **kwargs: Any):
        self.id = model_id
        self.name = name
        self.description = kwargs.get("description")
        self.schema = kwargs.get("schema")
        self.content_type = kwargs.get("content_type")

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "id": self.id,
            "name": self.name,
        }
        if self.description:
            dct["description"] = self.description
        if self.schema:
            dct["schema"] = self.schema
        if self.content_type:
            dct["contentType"] = self.content_type
        return dct


class BasePathMapping(BaseModel):
    # operations
    OPERATION_REPLACE = "replace"
    OPERATION_PATH = "path"
    OPERATION_VALUE = "value"
    OPERATION_OP = "op"

    def __init__(self, domain_name: str, rest_api_id: str, **kwargs: Any):
        self.domain_name = domain_name
        self.rest_api_id = rest_api_id
        self.base_path = kwargs.get("basePath") or "(none)"
        self.stage = kwargs.get("stage")

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "domain_name": self.domain_name,
            "restApiId": self.rest_api_id,
            "basePath": self.base_path,
        }
        if self.stage is not None:
            dct["stage"] = self.stage
        return dct

    def apply_patch_operations(self, patch_operations: List[Dict[str, Any]]) -> None:
        for op in patch_operations:
            path = op["path"]
            value = op["value"]
            operation = op["op"]
            if operation == self.OPERATION_REPLACE:
                if "/basePath" in path:
                    self.base_path = value
                if "/restapiId" in path:
                    self.rest_api_id = value
                if "/stage" in path:
                    self.stage = value


class GatewayResponse(BaseModel):
    def __init__(
        self,
        response_type: str,
        status_code: int,
        response_parameters: Dict[str, Any],
        response_templates: Dict[str, str],
    ):
        self.response_type = response_type
        self.default_response = False
        self.status_code = status_code
        self.response_parameters = response_parameters
        self.response_templates = response_templates

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "responseType": self.response_type,
            "defaultResponse": self.default_response,
        }
        if self.status_code is not None:
            dct["statusCode"] = self.status_code
        if self.response_parameters is not None:
            dct["responseParameters"] = self.response_parameters
        if self.response_templates is not None:
            dct["responseTemplates"] = self.response_templates
        return dct


class Account(BaseModel):
    def __init__(self) -> None:
        self.cloudwatch_role_arn: Optional[str] = None
        self.throttle_settings: Dict[str, Any] = {
            "burstLimit": 5000,
            "rateLimit": 10000.0,
        }
        self.features: Optional[List[str]] = None
        self.api_key_version: str = "1"

    def apply_patch_operations(
        self, patch_operations: List[Dict[str, Any]]
    ) -> "Account":
        for op in patch_operations:
            if "/cloudwatchRoleArn" in op["path"]:
                self.cloudwatch_role_arn = op["value"]
            elif "/features" in op["path"]:
                if op["op"] == "add":
                    if self.features is None:
                        self.features = [op["value"]]
                    else:
                        self.features.append(op["value"])
                elif op["op"] == "remove":
                    if op["value"] == "UsagePlans":
                        raise BadRequestException(
                            "Usage Plans cannot be disabled once enabled"
                        )
                    if self.features is not None:
                        self.features.remove(op["value"])
                else:
                    raise NotImplementedError(
                        f'Patch operation "{op["op"]}" for "/features" not implemented'
                    )
            else:
                raise NotImplementedError(
                    f'Patch operation "{op["op"]}" for "{op["path"]}" not implemented'
                )
        return self

    def to_json(self) -> Dict[str, Any]:
        return {
            "cloudwatchRoleArn": self.cloudwatch_role_arn,
            "throttleSettings": self.throttle_settings,
            "features": self.features,
            "apiKeyVersion": self.api_key_version,
        }


class APIGatewayBackend(BaseBackend):
    """
    API Gateway mock.

    The public URLs of an API integration are mocked as well, i.e. the following would be supported in Moto:

    .. sourcecode:: python

        client.put_integration(
            restApiId=api_id,
            ...,
            uri="http://httpbin.org/robots.txt",
            integrationHttpMethod="GET"
        )
        deploy_url = f"https://{api_id}.execute-api.us-east-1.amazonaws.com/dev"
        assert requests.get(deploy_url).content == b"a fake response"

    Limitations:
     - Integrations of type HTTP are supported
     - Integrations of type AWS with service DynamoDB are supported
     - Other types (AWS_PROXY, MOCK, etc) are ignored
     - Other services are not yet supported
     - The BasePath of an API is ignored
     - TemplateMapping is not yet supported for requests/responses
     - This only works when using the decorators, not in ServerMode
    """

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.account: Account = Account()
        self.apis: Dict[str, RestAPI] = {}
        self.keys: Dict[str, ApiKey] = {}
        self.usage_plans: Dict[str, UsagePlan] = {}
        self.usage_plan_keys: Dict[str, Dict[str, UsagePlanKey]] = {}
        self.domain_names: Dict[str, DomainName] = {}
        self.models: Dict[str, Model] = {}
        self.base_path_mappings: Dict[str, Dict[str, BasePathMapping]] = {}
        self.vpc_links: Dict[str, VpcLink] = {}

    def create_rest_api(
        self,
        name: str,
        description: str,
        api_key_source: Optional[str] = None,
        endpoint_configuration: Optional[str] = None,
        tags: Optional[List[Dict[str, str]]] = None,
        policy: Optional[str] = None,
        minimum_compression_size: Optional[int] = None,
        disable_execute_api_endpoint: Optional[bool] = None,
    ) -> RestAPI:
        api_id = ApigwRestApiIdentifier(
            self.account_id, self.region_name, name
        ).generate(tags=tags)
        rest_api = RestAPI(
            api_id,
            self.account_id,
            self.region_name,
            name,
            description,
            api_key_source=api_key_source,
            endpoint_configuration=endpoint_configuration,
            tags=tags,
            policy=policy,
            minimum_compression_size=minimum_compression_size,
            disable_execute_api_endpoint=disable_execute_api_endpoint,
        )
        self.apis[api_id] = rest_api
        return rest_api

    def import_rest_api(
        self, api_doc: Dict[str, Any], fail_on_warnings: bool
    ) -> RestAPI:
        """
        Only a subset of the OpenAPI spec 3.x is currently implemented.
        """
        if fail_on_warnings:
            try:
                validate(api_doc)  # type: ignore[arg-type]
            except OpenAPIValidationError as e:
                raise InvalidOpenAPIDocumentException(e)
        name = api_doc["info"]["title"]
        description = api_doc["info"]["description"]
        api = self.create_rest_api(name=name, description=description)
        self.put_rest_api(api.id, api_doc, fail_on_warnings=fail_on_warnings)
        return api

    def export_api(self, rest_api_id: str, export_type: str) -> Dict[str, Any]:
        """
        Not all fields are implemented yet.
        The export-type is currently ignored - we will only return the 'swagger'-format
        """
        try:
            api = self.get_rest_api(rest_api_id)
        except RestAPINotFound:
            raise StageNotFoundException
        if export_type not in ["swagger", "oas30"]:
            raise BadRequestException(f"No API exporter for type '{export_type}'")
        now = datetime.now().strftime("%Y-%m-%dT%H:%m:%S")
        resp: Dict[str, Any] = {
            "swagger": "2.0",
            "info": {"version": now, "title": api.name},
            "host": f"{api.id}.execute-api.{self.region_name}.amazonaws.com",
            "basePath": "/",
            "schemes": ["https"],
            "paths": {},
            "definitions": {"Empty": {"type": "object", "title": "Empty Schema"}},
        }
        for res in api.resources.values():
            path = res.get_path()
            resp["paths"][path] = {}
            for method_type, method in res.resource_methods.items():
                resp["paths"][path][method_type] = {
                    "produces": ["application/json"],
                    "responses": {},
                }
                for code, _ in method.method_responses.items():
                    resp["paths"][path][method_type]["responses"][code] = {
                        "description": f"{code} response",
                        "schema": {"$ref": "#/definitions/Empty"},
                    }
        return resp

    def get_rest_api(self, function_id: str) -> RestAPI:
        rest_api = self.apis.get(function_id)
        if rest_api is None:
            raise RestAPINotFound()
        return rest_api

    def put_rest_api(
        self,
        function_id: str,
        api_doc: Dict[str, Any],
        fail_on_warnings: bool,
        mode: str = "merge",
    ) -> RestAPI:
        """
        Only a subset of the OpenAPI spec 3.x is currently implemented.
        """
        if mode not in ["merge", "overwrite"]:
            raise InvalidOpenApiModeException()

        if api_doc.get("swagger") is not None or (
            api_doc.get("openapi") is not None and api_doc["openapi"][0] != "3"
        ):
            raise InvalidOpenApiDocVersionException()

        if fail_on_warnings:
            try:
                validate(api_doc)  # type: ignore[arg-type]
            except OpenAPIValidationError as e:
                raise InvalidOpenAPIDocumentException(e)

        if mode == "overwrite":
            api = self.get_rest_api(function_id)
            api.resources = {}
            api.default = api.add_child("/")  # Add default child

        for path, resource_doc in sorted(api_doc["paths"].items(), key=lambda x: x[0]):
            # We may want to create a path like /store/inventory
            # Ensure that /store exists first, so we can use it as a parent
            ancestors = path.split("/")[
                1:-1
            ]  # skip first (empty), skip last (child) - only process ancestors
            direct_parent = ""
            parent_id = self.apis[function_id].get_resource_for_path("/").id
            for a in ancestors:
                res = self.apis[function_id].get_resource_for_path(
                    direct_parent + "/" + a
                )
                if res is None:
                    res = self.create_resource(
                        function_id=function_id,
                        parent_resource_id=parent_id,
                        path_part=a,
                    )
                parent_id = res.id
                direct_parent = direct_parent + "/" + a

            # Now that we know all ancestors are created, create the resource itself
            parent_path_part = path[0 : path.rfind("/")] or "/"
            parent_resource_id = (
                self.apis[function_id].get_resource_for_path(parent_path_part).id
            )
            resource = self.create_resource(
                function_id=function_id,
                parent_resource_id=parent_resource_id,
                path_part=path[path.rfind("/") + 1 :],
            )

            for method_type, method_doc in resource_doc.items():
                method_type = method_type.upper()
                if method_doc.get("x-amazon-apigateway-integration") is None:
                    self.put_method(function_id, resource.id, method_type, None)
                    method_responses = method_doc.get("responses", {}).items()
                    for response_code, _ in method_responses:
                        self.put_method_response(
                            function_id,
                            resource.id,
                            method_type,
                            response_code,
                            response_models=None,
                            response_parameters=None,
                        )

        return self.get_rest_api(function_id)

    def update_rest_api(
        self, function_id: str, patch_operations: List[Dict[str, Any]]
    ) -> RestAPI:
        rest_api = self.apis.get(function_id)
        if rest_api is None:
            raise RestAPINotFound()
        self.apis[function_id].apply_patch_operations(patch_operations)
        return self.apis[function_id]

    def list_apis(self) -> List[RestAPI]:
        return list(self.apis.values())

    def delete_rest_api(self, function_id: str) -> RestAPI:
        rest_api = self.apis.pop(function_id)
        return rest_api

    def get_resources(self, function_id: str) -> List[Resource]:
        api = self.get_rest_api(function_id)
        return list(api.resources.values())

    def get_resource(self, function_id: str, resource_id: str) -> Resource:
        api = self.get_rest_api(function_id)
        if resource_id not in api.resources:
            raise ResourceIdNotFoundException
        return api.resources[resource_id]

    def create_resource(
        self, function_id: str, parent_resource_id: str, path_part: str
    ) -> Resource:
        api = self.get_rest_api(function_id)
        if not path_part:
            # We're attempting to create the default resource, which already exists.
            return api.default
        if not re.match("^\\{?[a-zA-Z0-9._\\-\\:]+\\+?\\}?$", path_part):
            raise InvalidResourcePathException()
        return api.add_child(path=path_part, parent_id=parent_resource_id)

    def delete_resource(self, function_id: str, resource_id: str) -> Resource:
        api = self.get_rest_api(function_id)
        return api.resources.pop(resource_id)

    def get_method(
        self, function_id: str, resource_id: str, method_type: str
    ) -> Method:
        resource = self.get_resource(function_id, resource_id)
        return resource.get_method(method_type)

    def put_method(
        self,
        function_id: str,
        resource_id: str,
        method_type: str,
        authorization_type: Optional[str],
        api_key_required: Optional[bool] = None,
        request_parameters: Optional[Dict[str, Any]] = None,
        request_models: Optional[Dict[str, Any]] = None,
        operation_name: Optional[str] = None,
        authorizer_id: Optional[str] = None,
        authorization_scopes: Optional[str] = None,
        request_validator_id: Optional[str] = None,
    ) -> Method:
        resource = self.get_resource(function_id, resource_id)
        method = resource.add_method(
            method_type,
            authorization_type,
            api_key_required=api_key_required,
            request_parameters=request_parameters,
            request_models=request_models,
            operation_name=operation_name,
            authorizer_id=authorizer_id,
            authorization_scopes=authorization_scopes,
            request_validator_id=request_validator_id,
        )
        return method

    def delete_method(
        self, function_id: str, resource_id: str, method_type: str
    ) -> None:
        resource = self.get_resource(function_id, resource_id)
        resource.delete_method(method_type)

    def get_authorizer(self, restapi_id: str, authorizer_id: str) -> Authorizer:
        api = self.get_rest_api(restapi_id)
        authorizer = api.authorizers.get(authorizer_id)
        if authorizer is None:
            raise AuthorizerNotFoundException()
        else:
            return authorizer

    def get_authorizers(self, restapi_id: str) -> List[Authorizer]:
        api = self.get_rest_api(restapi_id)
        return api.get_authorizers()

    def create_authorizer(
        self, restapi_id: str, name: str, authorizer_type: str, **kwargs: Any
    ) -> Authorizer:
        api = self.get_rest_api(restapi_id)
        authorizer_id = ApigwAuthorizerIdentifier(
            self.account_id, self.region_name, name
        ).generate()
        return api.create_authorizer(
            authorizer_id,
            name,
            authorizer_type,
            provider_arns=kwargs.get("provider_arns"),
            auth_type=kwargs.get("auth_type"),
            authorizer_uri=kwargs.get("authorizer_uri"),
            authorizer_credentials=kwargs.get("authorizer_credentials"),
            identity_source=kwargs.get("identity_source"),
            identiy_validation_expression=kwargs.get("identiy_validation_expression"),
            authorizer_result_ttl=kwargs.get("authorizer_result_ttl"),
        )

    def update_authorizer(
        self, restapi_id: str, authorizer_id: str, patch_operations: Any
    ) -> Authorizer:
        authorizer = self.get_authorizer(restapi_id, authorizer_id)
        return authorizer.apply_operations(patch_operations)

    def delete_authorizer(self, restapi_id: str, authorizer_id: str) -> None:
        api = self.get_rest_api(restapi_id)
        del api.authorizers[authorizer_id]

    def get_stage(self, function_id: str, stage_name: str) -> Stage:
        api = self.get_rest_api(function_id)
        stage = api.stages.get(stage_name)
        if stage is None:
            raise StageNotFoundException()
        return stage

    def get_stages(self, function_id: str) -> List[Stage]:
        api = self.get_rest_api(function_id)
        return api.get_stages()

    def create_stage(
        self,
        function_id: str,
        stage_name: str,
        deploymentId: str,
        variables: Optional[Any] = None,
        description: str = "",
        cacheClusterEnabled: Optional[bool] = None,
        cacheClusterSize: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        tracing_enabled: Optional[bool] = None,
    ) -> Stage:
        if variables is None:
            variables = {}
        api = self.get_rest_api(function_id)
        return api.create_stage(
            stage_name,
            deploymentId,
            variables=variables,
            description=description,
            cacheClusterEnabled=cacheClusterEnabled,
            cacheClusterSize=cacheClusterSize,
            tags=tags,
            tracing_enabled=tracing_enabled,
        )

    def update_stage(
        self, function_id: str, stage_name: str, patch_operations: Any
    ) -> Stage:
        stage = self.get_stage(function_id, stage_name)
        if not stage:
            api = self.get_rest_api(function_id)
            stage = api.stages[stage_name] = Stage()
        return stage.apply_operations(patch_operations)

    def delete_stage(self, function_id: str, stage_name: str) -> None:
        api = self.get_rest_api(function_id)
        deleted = api.stages.pop(stage_name, None)
        if not deleted:
            raise StageNotFoundException()

    def get_method_response(
        self, function_id: str, resource_id: str, method_type: str, response_code: str
    ) -> Optional[MethodResponse]:
        method = self.get_method(function_id, resource_id, method_type)
        return method.get_response(response_code)

    def put_method_response(
        self,
        function_id: str,
        resource_id: str,
        method_type: str,
        response_code: str,
        response_models: Any,
        response_parameters: Any,
    ) -> MethodResponse:
        method = self.get_method(function_id, resource_id, method_type)
        return method.create_response(
            response_code, response_models, response_parameters
        )

    def delete_method_response(
        self, function_id: str, resource_id: str, method_type: str, response_code: str
    ) -> Optional[MethodResponse]:
        method = self.get_method(function_id, resource_id, method_type)
        return method.delete_response(response_code)

    def put_integration(
        self,
        function_id: str,
        resource_id: str,
        method_type: str,
        integration_type: str,
        uri: str,
        integration_method: str,
        credentials: Optional[str] = None,
        request_templates: Optional[Dict[str, Any]] = None,
        passthrough_behavior: Optional[str] = None,
        tls_config: Optional[Dict[str, Any]] = None,
        cache_namespace: Optional[str] = None,
        timeout_in_millis: Optional[str] = None,
        request_parameters: Optional[Dict[str, Any]] = None,
        content_handling: Optional[str] = None,
        connection_type: Optional[str] = None,
    ) -> Integration:
        resource = self.get_resource(function_id, resource_id)
        if credentials and not re.match(
            f"^arn:{get_partition(self.region_name)}:iam::" + str(self.account_id),
            credentials,
        ):
            raise CrossAccountNotAllowed()
        if not integration_method and integration_type in [
            "HTTP",
            "HTTP_PROXY",
            "AWS",
            "AWS_PROXY",
        ]:
            raise IntegrationMethodNotDefined()
        if integration_type in ["AWS_PROXY"] and re.match(
            ARN_PARTITION_REGEX + ":apigateway:[a-zA-Z0-9-]+:s3", uri
        ):
            raise AwsProxyNotAllowed()
        if (
            integration_type in ["AWS"]
            and re.match(ARN_PARTITION_REGEX + ":apigateway:[a-zA-Z0-9-]+:s3", uri)
            and not credentials
        ):
            raise RoleNotSpecified()
        if integration_type in ["HTTP", "HTTP_PROXY"] and not self._uri_validator(uri):
            raise InvalidHttpEndpoint()
        if integration_type in ["AWS", "AWS_PROXY"] and not re.match(
            ARN_PARTITION_REGEX + ":", uri
        ):
            raise InvalidArn()
        if integration_type in ["AWS", "AWS_PROXY"] and not re.match(
            ARN_PARTITION_REGEX
            + ":apigateway:[a-zA-Z0-9-]+:[a-zA-Z0-9-.]+:(path|action)/",
            uri,
        ):
            raise InvalidIntegrationArn()
        integration = resource.add_integration(
            method_type,
            integration_type,
            uri,
            integration_method=integration_method,
            request_templates=request_templates,
            passthrough_behavior=passthrough_behavior,
            tls_config=tls_config,
            cache_namespace=cache_namespace,
            timeout_in_millis=timeout_in_millis,
            request_parameters=request_parameters,
            content_handling=content_handling,
            credentials=credentials,
            connection_type=connection_type,
        )
        return integration

    def get_integration(
        self, function_id: str, resource_id: str, method_type: str
    ) -> Integration:
        resource = self.get_resource(function_id, resource_id)
        return resource.get_integration(method_type)  # type: ignore[return-value]

    def delete_integration(
        self, function_id: str, resource_id: str, method_type: str
    ) -> Integration:
        resource = self.get_resource(function_id, resource_id)
        return resource.delete_integration(method_type)

    def put_integration_response(
        self,
        function_id: str,
        resource_id: str,
        method_type: str,
        status_code: str,
        selection_pattern: str,
        response_templates: Dict[str, str],
        response_parameters: Dict[str, str],
        content_handling: str,
    ) -> IntegrationResponse:
        integration = self.get_integration(function_id, resource_id, method_type)
        if integration:
            return integration.create_integration_response(
                status_code,
                selection_pattern,
                response_templates,
                response_parameters,
                content_handling,
            )
        raise NoIntegrationResponseDefined()

    def get_integration_response(
        self, function_id: str, resource_id: str, method_type: str, status_code: str
    ) -> IntegrationResponse:
        integration = self.get_integration(function_id, resource_id, method_type)
        integration_response = integration.get_integration_response(status_code)
        return integration_response

    def delete_integration_response(
        self, function_id: str, resource_id: str, method_type: str, status_code: str
    ) -> IntegrationResponse:
        integration = self.get_integration(function_id, resource_id, method_type)
        integration_response = integration.delete_integration_response(status_code)
        return integration_response

    def create_deployment(
        self,
        function_id: str,
        name: str,
        description: str = "",
        stage_variables: Any = None,
    ) -> Deployment:
        if stage_variables is None:
            stage_variables = {}
        api = self.get_rest_api(function_id)
        nested_methods = [
            list(res.resource_methods.values())
            for res in self.get_resources(function_id)
        ]
        methods = [m for sublist in nested_methods for m in sublist]
        if not any(methods):
            raise NoMethodDefined()
        method_integrations = [
            method.method_integration for method in methods if method.method_integration
        ]
        if not any(method_integrations):
            raise NoIntegrationDefined()
        deployment = api.create_deployment(name, description, stage_variables)
        return deployment

    def get_deployment(self, function_id: str, deployment_id: str) -> Deployment:
        api = self.get_rest_api(function_id)
        return api.get_deployment(deployment_id)

    def get_deployments(self, function_id: str) -> List[Deployment]:
        api = self.get_rest_api(function_id)
        return api.get_deployments()

    def delete_deployment(self, function_id: str, deployment_id: str) -> Deployment:
        api = self.get_rest_api(function_id)
        return api.delete_deployment(deployment_id)

    def create_api_key(self, payload: Dict[str, Any]) -> ApiKey:
        if payload.get("value"):
            if len(payload.get("value", [])) < 20:
                raise ApiKeyValueMinLength()
            for api_key in self.get_api_keys():
                if api_key.value == payload["value"]:
                    raise ApiKeyAlreadyExists()
        api_key_id = ApigwApiKeyIdentifier(
            self.account_id,
            self.region_name,
            # The value of an api key must be unique on aws
            payload.get("value", ""),
        ).generate()
        key = ApiKey(api_key_id=api_key_id, **payload)
        self.keys[key.id] = key
        return key

    def get_api_keys(self, name: Optional[str] = None) -> List[ApiKey]:
        return [
            key
            for key in self.keys.values()
            if not name or (key.name and key.name.startswith(name))
        ]

    def get_api_key(self, api_key_id: str) -> ApiKey:
        if api_key_id not in self.keys:
            raise ApiKeyNotFoundException()
        return self.keys[api_key_id]

    def update_api_key(self, api_key_id: str, patch_operations: Any) -> ApiKey:
        key = self.keys[api_key_id]
        return key.update_operations(patch_operations)

    def delete_api_key(self, api_key_id: str) -> None:
        self.keys.pop(api_key_id)

    def create_usage_plan(self, payload: Any) -> UsagePlan:
        usage_plan_id = ApigwUsagePlanIdentifier(
            self.account_id, self.region_name, payload["name"]
        ).generate()
        plan = UsagePlan(usage_plan_id=usage_plan_id, **payload)
        self.usage_plans[plan.id] = plan
        return plan

    def get_usage_plans(self, api_key_id: Optional[str] = None) -> List[UsagePlan]:
        plans = list(self.usage_plans.values())
        if api_key_id is not None:
            plans = [
                plan
                for plan in plans
                if dict(self.usage_plan_keys.get(plan.id, {})).get(api_key_id)
            ]
        return plans

    def get_usage_plan(self, usage_plan_id: str) -> UsagePlan:
        if usage_plan_id not in self.usage_plans:
            raise UsagePlanNotFoundException()
        return self.usage_plans[usage_plan_id]

    def update_usage_plan(self, usage_plan_id: str, patch_operations: Any) -> UsagePlan:
        """
        The following PatchOperations are currently supported:
        add    : Everything except /apiStages/{apidId:stageName}/throttle/ and children
        replace: Everything except /apiStages/{apidId:stageName}/throttle/ and children
        remove : Everything except /apiStages/{apidId:stageName}/throttle/ and children
        copy   : Nothing yet
        """
        if usage_plan_id not in self.usage_plans:
            raise UsagePlanNotFoundException()
        self.usage_plans[usage_plan_id].apply_patch_operations(patch_operations)
        return self.usage_plans[usage_plan_id]

    def delete_usage_plan(self, usage_plan_id: str) -> None:
        self.usage_plans.pop(usage_plan_id)

    def create_usage_plan_key(
        self, usage_plan_id: str, payload: Dict[str, Any]
    ) -> UsagePlanKey:
        if usage_plan_id not in self.usage_plan_keys:
            self.usage_plan_keys[usage_plan_id] = {}

        key_id = payload["keyId"]
        if key_id not in self.keys:
            raise ApiKeyNotFoundException()

        api_key = self.keys[key_id]

        usage_plan_key = UsagePlanKey(
            plan_id=key_id,
            plan_type=payload["keyType"],
            name=api_key.name,
            value=api_key.value,
        )
        self.usage_plan_keys[usage_plan_id][usage_plan_key.id] = usage_plan_key
        return usage_plan_key

    def get_usage_plan_keys(
        self, usage_plan_id: str, name: Optional[str] = None
    ) -> List[UsagePlanKey]:
        if usage_plan_id not in self.usage_plan_keys:
            return []

        plan_keys = self.usage_plan_keys[usage_plan_id].values()
        if name:
            return [
                key
                for key in plan_keys
                if not name or (key.name and key.name.startswith(name))
            ]
        return list(plan_keys)

    def get_usage_plan_key(self, usage_plan_id: str, key_id: str) -> UsagePlanKey:
        # first check if is a valid api key
        if key_id not in self.keys:
            raise ApiKeyNotFoundException()

        # then check if is a valid api key and that the key is in the plan
        if (
            usage_plan_id not in self.usage_plan_keys
            or key_id not in self.usage_plan_keys[usage_plan_id]
        ):
            raise UsagePlanNotFoundException()

        return self.usage_plan_keys[usage_plan_id][key_id]

    def delete_usage_plan_key(self, usage_plan_id: str, key_id: str) -> None:
        self.usage_plan_keys[usage_plan_id].pop(key_id)

    def _uri_validator(self, uri: str) -> bool:
        try:
            result = urlparse(uri)
            return all([result.scheme, result.netloc, result.path or "/"])
        except Exception:
            return False

    def create_domain_name(
        self,
        domain_name: str,
        certificate_name: str,
        tags: List[Dict[str, str]],
        certificate_arn: str,
        certificate_body: str,
        certificate_private_key: str,
        certificate_chain: str,
        regional_certificate_name: str,
        regional_certificate_arn: str,
        endpoint_configuration: Any,
        security_policy: str,
    ) -> DomainName:
        if not domain_name:
            raise InvalidDomainName()

        new_domain_name = DomainName(
            domain_name=domain_name,
            certificate_name=certificate_name,
            certificate_private_key=certificate_private_key,
            certificate_arn=certificate_arn,
            certificate_body=certificate_body,
            certificate_chain=certificate_chain,
            regional_certificate_name=regional_certificate_name,
            regional_certificate_arn=regional_certificate_arn,
            endpoint_configuration=endpoint_configuration,
            tags=tags,
            security_policy=security_policy,
            region_name=self.region_name,
        )

        self.domain_names[domain_name] = new_domain_name
        return new_domain_name

    def get_domain_names(self) -> List[DomainName]:
        return list(self.domain_names.values())

    def get_domain_name(self, domain_name: str) -> DomainName:
        domain_info = self.domain_names.get(domain_name)
        if domain_info is None:
            raise DomainNameNotFound()
        else:
            return domain_info

    def delete_domain_name(self, domain_name: str) -> None:
        domain_info = self.domain_names.pop(domain_name, None)
        if domain_info is None:
            raise DomainNameNotFound()

    def create_model(
        self,
        rest_api_id: str,
        name: str,
        content_type: str,
        description: str,
        schema: str,
    ) -> Model:
        if not rest_api_id:
            raise InvalidRestApiId
        if not name:
            raise InvalidModelName

        api = self.get_rest_api(rest_api_id)
        new_model = api.add_model(
            name=name,
            description=description,
            schema=schema,
            content_type=content_type,
        )

        return new_model

    def get_models(self, rest_api_id: str) -> List[Model]:
        if not rest_api_id:
            raise InvalidRestApiId
        api = self.get_rest_api(rest_api_id)
        models = api.models.values()
        return list(models)

    def get_model(self, rest_api_id: str, model_name: str) -> Model:
        if not rest_api_id:
            raise InvalidRestApiId
        api = self.get_rest_api(rest_api_id)
        model = api.models.get(model_name)
        if model is None:
            raise ModelNotFound
        else:
            return model

    def get_request_validators(self, restapi_id: str) -> List[RequestValidator]:
        restApi = self.get_rest_api(restapi_id)
        return restApi.get_request_validators()

    def create_request_validator(
        self, restapi_id: str, name: str, body: Optional[bool], params: Any
    ) -> RequestValidator:
        restApi = self.get_rest_api(restapi_id)
        return restApi.create_request_validator(
            name=name, validateRequestBody=body, validateRequestParameters=params
        )

    def get_request_validator(
        self, restapi_id: str, validator_id: str
    ) -> RequestValidator:
        restApi = self.get_rest_api(restapi_id)
        return restApi.get_request_validator(validator_id)

    def delete_request_validator(self, restapi_id: str, validator_id: str) -> None:
        restApi = self.get_rest_api(restapi_id)
        restApi.delete_request_validator(validator_id)

    def update_request_validator(
        self, restapi_id: str, validator_id: str, patch_operations: Any
    ) -> RequestValidator:
        restApi = self.get_rest_api(restapi_id)
        return restApi.update_request_validator(validator_id, patch_operations)

    def create_base_path_mapping(
        self, domain_name: str, rest_api_id: str, base_path: str, stage: str
    ) -> BasePathMapping:
        if domain_name not in self.domain_names:
            raise DomainNameNotFound()

        if base_path and "/" in base_path:
            raise InvalidBasePathException()

        if rest_api_id not in self.apis:
            raise InvalidRestApiIdForBasePathMappingException()

        if stage and self.apis[rest_api_id].stages.get(stage) is None:
            raise InvalidStageException()

        new_base_path_mapping = BasePathMapping(
            domain_name=domain_name,
            rest_api_id=rest_api_id,
            basePath=base_path,
            stage=stage,
        )

        new_base_path = new_base_path_mapping.base_path
        if self.base_path_mappings.get(domain_name) is None:
            self.base_path_mappings[domain_name] = {}
        else:
            if (
                self.base_path_mappings[domain_name].get(new_base_path)
                and new_base_path != "(none)"
            ):
                raise BasePathConflictException()
        self.base_path_mappings[domain_name][new_base_path] = new_base_path_mapping
        return new_base_path_mapping

    def get_base_path_mappings(self, domain_name: str) -> List[BasePathMapping]:
        if domain_name not in self.domain_names:
            raise DomainNameNotFound()

        return list(self.base_path_mappings[domain_name].values())

    def get_base_path_mapping(
        self, domain_name: str, base_path: str
    ) -> BasePathMapping:
        if domain_name not in self.domain_names:
            raise DomainNameNotFound()

        if base_path not in self.base_path_mappings[domain_name]:
            raise BasePathNotFoundException()

        return self.base_path_mappings[domain_name][base_path]

    def delete_base_path_mapping(self, domain_name: str, base_path: str) -> None:
        if domain_name not in self.domain_names:
            raise DomainNameNotFound()

        if base_path not in self.base_path_mappings[domain_name]:
            raise BasePathNotFoundException()

        self.base_path_mappings[domain_name].pop(base_path)

    def update_base_path_mapping(
        self, domain_name: str, base_path: str, patch_operations: Any
    ) -> BasePathMapping:
        if domain_name not in self.domain_names:
            raise DomainNameNotFound()

        if base_path not in self.base_path_mappings[domain_name]:
            raise BasePathNotFoundException()

        base_path_mapping = self.get_base_path_mapping(domain_name, base_path)

        rest_api_ids = [
            op["value"] for op in patch_operations if op["path"] == "/restapiId"
        ]
        if len(rest_api_ids) == 0:
            modified_rest_api_id = base_path_mapping.rest_api_id
        else:
            modified_rest_api_id = rest_api_ids[-1]

        stages = [op["value"] for op in patch_operations if op["path"] == "/stage"]
        if len(stages) == 0:
            modified_stage = base_path_mapping.stage
        else:
            modified_stage = stages[-1]

        base_paths = [
            op["value"] for op in patch_operations if op["path"] == "/basePath"
        ]
        if len(base_paths) == 0:
            modified_base_path = base_path_mapping.base_path
        else:
            modified_base_path = base_paths[-1]

        rest_api = self.apis.get(modified_rest_api_id)
        if rest_api is None:
            raise InvalidRestApiIdForBasePathMappingException()
        if modified_stage and rest_api.stages.get(modified_stage) is None:
            raise InvalidStageException()

        base_path_mapping.apply_patch_operations(patch_operations)

        if base_path != modified_base_path:
            self.base_path_mappings[domain_name].pop(base_path)
            self.base_path_mappings[domain_name][modified_base_path] = base_path_mapping

        return base_path_mapping

    def create_vpc_link(
        self,
        name: str,
        description: str,
        target_arns: List[str],
        tags: List[Dict[str, str]],
    ) -> VpcLink:
        vpc_link_id = ApigwVpcLinkIdentifier(
            self.account_id, self.region_name, name
        ).generate()
        vpc_link = VpcLink(
            vpc_link_id,
            name,
            description=description,
            target_arns=target_arns,
            tags=tags,
        )
        self.vpc_links[vpc_link.id] = vpc_link
        return vpc_link

    def delete_vpc_link(self, vpc_link_id: str) -> None:
        self.vpc_links.pop(vpc_link_id, None)

    def get_vpc_link(self, vpc_link_id: str) -> VpcLink:
        if vpc_link_id not in self.vpc_links:
            raise VpcLinkNotFound
        return self.vpc_links[vpc_link_id]

    def get_vpc_links(self) -> List[VpcLink]:
        """
        Pagination has not yet been implemented
        """
        return list(self.vpc_links.values())

    def put_gateway_response(
        self,
        rest_api_id: str,
        response_type: str,
        status_code: int,
        response_parameters: Any,
        response_templates: Any,
    ) -> GatewayResponse:
        api = self.get_rest_api(rest_api_id)
        response = api.put_gateway_response(
            response_type,
            status_code=status_code,
            response_parameters=response_parameters,
            response_templates=response_templates,
        )
        return response

    def get_gateway_response(
        self, rest_api_id: str, response_type: str
    ) -> GatewayResponse:
        api = self.get_rest_api(rest_api_id)
        return api.get_gateway_response(response_type)

    def get_gateway_responses(self, rest_api_id: str) -> List[GatewayResponse]:
        """
        Pagination is not yet implemented
        """
        api = self.get_rest_api(rest_api_id)
        return api.get_gateway_responses()

    def delete_gateway_response(self, rest_api_id: str, response_type: str) -> None:
        api = self.get_rest_api(rest_api_id)
        api.delete_gateway_response(response_type)

    def update_account(self, patch_operations: List[Dict[str, Any]]) -> Account:
        account = self.account.apply_patch_operations(patch_operations)
        return account

    def get_account(self) -> Account:
        return self.account


apigateway_backends = BackendDict(APIGatewayBackend, "apigateway")
