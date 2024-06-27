"""ApiGatewayV2Backend class with methods for supported APIs."""

import hashlib
import string
from typing import Any, Dict, List, Optional, Union

import yaml

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random as random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    ApiMappingNotFound,
    ApiNotFound,
    AuthorizerNotFound,
    BadRequestException,
    DomainNameAlreadyExists,
    DomainNameNotFound,
    IntegrationNotFound,
    IntegrationResponseNotFound,
    ModelNotFound,
    RouteNotFound,
    RouteResponseNotFound,
    StageNotFound,
    VpcLinkNotFound,
)


class Stage(BaseModel):
    def __init__(self, api: "Api", config: Dict[str, Any]):
        self.config = config
        self.name = config["stageName"]
        if api.protocol_type == "HTTP":
            self.default_route_settings = config.get(
                "defaultRouteSettings", {"detailedMetricsEnabled": False}
            )
        elif api.protocol_type == "WEBSOCKET":
            self.default_route_settings = config.get(
                "defaultRouteSettings",
                {
                    "dataTraceEnabled": False,
                    "detailedMetricsEnabled": False,
                    "loggingLevel": "OFF",
                },
            )
        self.access_log_settings = config.get("accessLogSettings")
        self.auto_deploy = config.get("autoDeploy")
        self.client_certificate_id = config.get("clientCertificateId")
        self.description = config.get("description")
        self.route_settings = config.get("routeSettings", {})
        self.stage_variables = config.get("stageVariables", {})
        self.tags = config.get("tags", {})
        self.created = self.updated = unix_time()

    def to_json(self) -> Dict[str, Any]:
        dct = {
            "stageName": self.name,
            "defaultRouteSettings": self.default_route_settings,
            "createdDate": self.created,
            "lastUpdatedDate": self.updated,
            "routeSettings": self.route_settings,
            "stageVariables": self.stage_variables,
            "tags": self.tags,
        }
        if self.access_log_settings:
            dct["accessLogSettings"] = self.access_log_settings
        if self.auto_deploy is not None:
            dct["autoDeploy"] = self.auto_deploy
        if self.client_certificate_id:
            dct["clientCertificateId"] = self.client_certificate_id
        if self.description:
            dct["description"] = self.description
        return dct


class Authorizer(BaseModel):
    def __init__(
        self,
        auth_creds_arn: str,
        auth_payload_format_version: str,
        auth_result_ttl: str,
        authorizer_type: str,
        authorizer_uri: str,
        enable_simple_response: str,
        identity_source: str,
        identity_validation_expr: str,
        jwt_config: str,
        name: str,
    ):
        self.id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.auth_creds_arn = auth_creds_arn
        self.auth_payload_format_version = auth_payload_format_version
        self.auth_result_ttl = auth_result_ttl
        self.authorizer_type = authorizer_type
        self.authorizer_uri = authorizer_uri
        self.enable_simple_response = enable_simple_response
        self.identity_source = identity_source
        self.identity_validation_expr = identity_validation_expr
        self.jwt_config = jwt_config
        self.name = name

    def update(
        self,
        auth_creds_arn: str,
        auth_payload_format_version: str,
        auth_result_ttl: str,
        authorizer_type: str,
        authorizer_uri: str,
        enable_simple_response: str,
        identity_source: str,
        identity_validation_expr: str,
        jwt_config: str,
        name: str,
    ) -> None:
        if auth_creds_arn is not None:
            self.auth_creds_arn = auth_creds_arn
        if auth_payload_format_version is not None:
            self.auth_payload_format_version = auth_payload_format_version
        if auth_result_ttl is not None:
            self.auth_result_ttl = auth_result_ttl
        if authorizer_type is not None:
            self.authorizer_type = authorizer_type
        if authorizer_uri is not None:
            self.authorizer_uri = authorizer_uri
        if enable_simple_response is not None:
            self.enable_simple_response = enable_simple_response
        if identity_source is not None:
            self.identity_source = identity_source
        if identity_validation_expr is not None:
            self.identity_validation_expr = identity_validation_expr
        if jwt_config is not None:
            self.jwt_config = jwt_config
        if name is not None:
            self.name = name

    def to_json(self) -> Dict[str, Any]:
        return {
            "authorizerId": self.id,
            "authorizerCredentialsArn": self.auth_creds_arn,
            "authorizerPayloadFormatVersion": self.auth_payload_format_version,
            "authorizerResultTtlInSeconds": self.auth_result_ttl,
            "authorizerType": self.authorizer_type,
            "authorizerUri": self.authorizer_uri,
            "enableSimpleResponses": self.enable_simple_response,
            "identitySource": self.identity_source,
            "identityValidationExpression": self.identity_validation_expr,
            "jwtConfiguration": self.jwt_config,
            "name": self.name,
        }


class Integration(BaseModel):
    def __init__(
        self,
        connection_id: Optional[str],
        connection_type: str,
        content_handling_strategy: Optional[str],
        credentials_arn: Optional[str],
        description: str,
        integration_method: str,
        integration_type: str,
        integration_uri: str,
        passthrough_behavior: Optional[str],
        payload_format_version: Optional[str],
        integration_subtype: Optional[str],
        request_parameters: Optional[Dict[str, str]],
        request_templates: Optional[Dict[str, str]],
        response_parameters: Optional[Dict[str, Dict[str, str]]],
        template_selection_expression: Optional[str],
        timeout_in_millis: Optional[str],
        tls_config: Optional[Dict[str, str]],
    ):
        self.id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.connection_id = connection_id
        self.connection_type = connection_type
        self.content_handling_strategy = content_handling_strategy
        self.credentials_arn = credentials_arn
        self.description = description
        self.integration_method = integration_method
        self.integration_response_selection_expression = None
        self.integration_type = integration_type
        self.integration_subtype = integration_subtype
        self.integration_uri = integration_uri
        self.passthrough_behavior = passthrough_behavior
        self.payload_format_version = payload_format_version
        self.request_parameters = request_parameters
        self.request_templates = request_templates
        self.response_parameters = response_parameters
        self.template_selection_expression = template_selection_expression
        self.timeout_in_millis = int(timeout_in_millis) if timeout_in_millis else None
        self.tls_config = tls_config

        if self.integration_type in ["MOCK", "HTTP"]:
            self.integration_response_selection_expression = (
                "${integration.response.statuscode}"
            )
        elif self.integration_type in ["AWS"]:
            self.integration_response_selection_expression = (
                "${integration.response.body.errorMessage}"
            )
        if (
            self.integration_type in ["AWS", "MOCK", "HTTP"]
            and self.passthrough_behavior is None
        ):
            self.passthrough_behavior = "WHEN_NO_MATCH"
        if self.integration_uri is not None and self.integration_method is None:
            self.integration_method = "POST"
        if self.integration_type in ["AWS", "MOCK"]:
            self.timeout_in_millis = self.timeout_in_millis or 29000
        else:
            self.timeout_in_millis = self.timeout_in_millis or 30000

        self.responses: Dict[str, IntegrationResponse] = dict()

    def create_response(
        self,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> "IntegrationResponse":
        response = IntegrationResponse(
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )
        self.responses[response.id] = response
        return response

    def delete_response(self, integration_response_id: str) -> None:
        self.responses.pop(integration_response_id)

    def get_response(self, integration_response_id: str) -> "IntegrationResponse":
        if integration_response_id not in self.responses:
            raise IntegrationResponseNotFound(integration_response_id)
        return self.responses[integration_response_id]

    def get_responses(self) -> List["IntegrationResponse"]:
        return list(self.responses.values())

    def update_response(
        self,
        integration_response_id: str,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> "IntegrationResponse":
        int_response = self.responses[integration_response_id]
        int_response.update(
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )
        return int_response

    def update(
        self,
        connection_id: str,
        connection_type: str,
        content_handling_strategy: str,
        credentials_arn: str,
        description: str,
        integration_method: str,
        integration_type: str,
        integration_uri: str,
        passthrough_behavior: str,
        payload_format_version: str,
        integration_subtype: str,
        request_parameters: Dict[str, str],
        request_templates: Dict[str, str],
        response_parameters: Dict[str, Dict[str, str]],
        template_selection_expression: str,
        timeout_in_millis: Optional[int],
        tls_config: Dict[str, str],
    ) -> None:
        if connection_id is not None:
            self.connection_id = connection_id
        if connection_type is not None:
            self.connection_type = connection_type
        if content_handling_strategy is not None:
            self.content_handling_strategy = content_handling_strategy
        if credentials_arn is not None:
            self.credentials_arn = credentials_arn
        if description is not None:
            self.description = description
        if integration_method is not None:
            self.integration_method = integration_method
        if integration_type is not None:
            self.integration_type = integration_type
        if integration_uri is not None:
            self.integration_uri = integration_uri
        if passthrough_behavior is not None:
            self.passthrough_behavior = passthrough_behavior
        if payload_format_version is not None:
            self.payload_format_version = payload_format_version
        if integration_subtype is not None:
            self.integration_subtype = integration_subtype
        if request_parameters is not None:
            # Skip parameters with an empty value
            req_params = {
                key: value for (key, value) in request_parameters.items() if value
            }
            self.request_parameters = req_params
        if request_templates is not None:
            self.request_templates = request_templates
        if response_parameters is not None:
            self.response_parameters = response_parameters
        if template_selection_expression is not None:
            self.template_selection_expression = template_selection_expression
        if timeout_in_millis is not None:
            self.timeout_in_millis = timeout_in_millis
        if tls_config is not None:
            self.tls_config = tls_config

    def to_json(self) -> Dict[str, Any]:
        return {
            "connectionId": self.connection_id,
            "connectionType": self.connection_type,
            "contentHandlingStrategy": self.content_handling_strategy,
            "credentialsArn": self.credentials_arn,
            "description": self.description,
            "integrationId": self.id,
            "integrationMethod": self.integration_method,
            "integrationResponseSelectionExpression": self.integration_response_selection_expression,
            "integrationType": self.integration_type,
            "integrationSubtype": self.integration_subtype,
            "integrationUri": self.integration_uri,
            "passthroughBehavior": self.passthrough_behavior,
            "payloadFormatVersion": self.payload_format_version,
            "requestParameters": self.request_parameters,
            "requestTemplates": self.request_templates,
            "responseParameters": self.response_parameters,
            "templateSelectionExpression": self.template_selection_expression,
            "timeoutInMillis": self.timeout_in_millis,
            "tlsConfig": self.tls_config,
        }


class IntegrationResponse(BaseModel):
    def __init__(
        self,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ):
        self.id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.content_handling_strategy = content_handling_strategy
        self.integration_response_key = integration_response_key
        self.response_parameters = response_parameters
        self.response_templates = response_templates
        self.template_selection_expression = template_selection_expression

    def update(
        self,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> None:
        if content_handling_strategy is not None:
            self.content_handling_strategy = content_handling_strategy
        if integration_response_key is not None:
            self.integration_response_key = integration_response_key
        if response_parameters is not None:
            self.response_parameters = response_parameters
        if response_templates is not None:
            self.response_templates = response_templates
        if template_selection_expression is not None:
            self.template_selection_expression = template_selection_expression

    def to_json(self) -> Dict[str, str]:
        return {
            "integrationResponseId": self.id,
            "integrationResponseKey": self.integration_response_key,
            "contentHandlingStrategy": self.content_handling_strategy,
            "responseParameters": self.response_parameters,
            "responseTemplates": self.response_templates,
            "templateSelectionExpression": self.template_selection_expression,
        }


class Model(BaseModel):
    def __init__(self, content_type: str, description: str, name: str, schema: str):
        self.id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.content_type = content_type
        self.description = description
        self.name = name
        self.schema = schema

    def update(
        self, content_type: str, description: str, name: str, schema: str
    ) -> None:
        if content_type is not None:
            self.content_type = content_type
        if description is not None:
            self.description = description
        if name is not None:
            self.name = name
        if schema is not None:
            self.schema = schema

    def to_json(self) -> Dict[str, str]:
        return {
            "modelId": self.id,
            "contentType": self.content_type,
            "description": self.description,
            "name": self.name,
            "schema": self.schema,
        }


class RouteResponse(BaseModel):
    def __init__(
        self,
        route_response_key: str,
        model_selection_expression: str,
        response_models: str,
    ):
        self.id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.route_response_key = route_response_key
        self.model_selection_expression = model_selection_expression
        self.response_models = response_models

    def to_json(self) -> Dict[str, str]:
        return {
            "modelSelectionExpression": self.model_selection_expression,
            "responseModels": self.response_models,
            "routeResponseId": self.id,
            "routeResponseKey": self.route_response_key,
        }


class Route(BaseModel):
    def __init__(
        self,
        api_key_required: bool,
        authorization_scopes: List[str],
        authorization_type: Optional[str],
        authorizer_id: Optional[str],
        model_selection_expression: Optional[str],
        operation_name: Optional[str],
        request_models: Optional[Dict[str, str]],
        request_parameters: Optional[Dict[str, Dict[str, bool]]],
        route_key: str,
        route_response_selection_expression: Optional[str],
        target: str,
    ):
        self.route_id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.api_key_required = api_key_required
        self.authorization_scopes = authorization_scopes
        self.authorization_type = authorization_type
        self.authorizer_id = authorizer_id
        self.model_selection_expression = model_selection_expression
        self.operation_name = operation_name
        self.request_models = request_models
        self.request_parameters = request_parameters or {}
        self.route_key = route_key
        self.route_response_selection_expression = route_response_selection_expression
        self.target = target

        self.route_responses: Dict[str, RouteResponse] = dict()

    def create_route_response(
        self,
        route_response_key: str,
        model_selection_expression: str,
        response_models: str,
    ) -> RouteResponse:
        route_response = RouteResponse(
            route_response_key,
            model_selection_expression=model_selection_expression,
            response_models=response_models,
        )
        self.route_responses[route_response.id] = route_response
        return route_response

    def get_route_response(self, route_response_id: str) -> RouteResponse:
        if route_response_id not in self.route_responses:
            raise RouteResponseNotFound(route_response_id)
        return self.route_responses[route_response_id]

    def delete_route_response(self, route_response_id: str) -> None:
        self.route_responses.pop(route_response_id, None)

    def delete_route_request_parameter(self, request_param: str) -> None:
        del self.request_parameters[request_param]

    def update(
        self,
        api_key_required: Optional[bool],
        authorization_scopes: Optional[List[str]],
        authorization_type: str,
        authorizer_id: str,
        model_selection_expression: str,
        operation_name: str,
        request_models: Dict[str, str],
        request_parameters: Dict[str, Dict[str, bool]],
        route_key: str,
        route_response_selection_expression: str,
        target: str,
    ) -> None:
        if api_key_required is not None:
            self.api_key_required = api_key_required
        if authorization_scopes:
            self.authorization_scopes = authorization_scopes
        if authorization_type:
            self.authorization_type = authorization_type
        if authorizer_id is not None:
            self.authorizer_id = authorizer_id
        if model_selection_expression:
            self.model_selection_expression = model_selection_expression
        if operation_name is not None:
            self.operation_name = operation_name
        if request_models:
            self.request_models = request_models
        if request_parameters:
            self.request_parameters = request_parameters
        if route_key:
            self.route_key = route_key
        if route_response_selection_expression is not None:
            self.route_response_selection_expression = (
                route_response_selection_expression
            )
        if target:
            self.target = target

    def to_json(self) -> Dict[str, Any]:
        return {
            "apiKeyRequired": self.api_key_required,
            "authorizationScopes": self.authorization_scopes,
            "authorizationType": self.authorization_type,
            "authorizerId": self.authorizer_id,
            "modelSelectionExpression": self.model_selection_expression,
            "operationName": self.operation_name,
            "requestModels": self.request_models,
            "requestParameters": self.request_parameters,
            "routeId": self.route_id,
            "routeKey": self.route_key,
            "routeResponseSelectionExpression": self.route_response_selection_expression,
            "target": self.target,
        }


class Api(BaseModel):
    def __init__(
        self,
        region: str,
        name: str,
        api_key_selection_expression: str,
        cors_configuration: Optional[str],
        description: str,
        disable_execute_api_endpoint: str,
        disable_schema_validation: str,
        protocol_type: str,
        route_selection_expression: str,
        tags: Dict[str, str],
        version: str,
        backend: "ApiGatewayV2Backend",
    ):
        self.api_id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.api_endpoint = f"https://{self.api_id}.execute-api.{region}.amazonaws.com"
        self.backend = backend
        self.name = name
        self.api_key_selection_expression = (
            api_key_selection_expression or "$request.header.x-api-key"
        )
        self.created_date = unix_time()
        self.cors_configuration = cors_configuration
        self.description = description
        self.disable_execute_api_endpoint = disable_execute_api_endpoint or False
        self.disable_schema_validation = disable_schema_validation
        self.protocol_type = protocol_type
        self.route_selection_expression = (
            route_selection_expression or "$request.method $request.path"
        )
        self.version = version

        self.authorizers: Dict[str, Authorizer] = dict()
        self.integrations: Dict[str, Integration] = dict()
        self.models: Dict[str, Model] = dict()
        self.routes: Dict[str, Route] = dict()
        self.stages: Dict[str, Stage] = dict()

        self.arn = (
            f"arn:{get_partition(region)}:apigateway:{region}::/apis/{self.api_id}"
        )
        self.backend.tag_resource(self.arn, tags)

    def clear(self) -> None:
        self.authorizers.clear()
        self.integrations.clear()
        self.models.clear()
        self.routes.clear()
        self.stages.clear()

    def delete_cors_configuration(self) -> None:
        self.cors_configuration = None

    def create_authorizer(
        self,
        auth_creds_arn: str,
        auth_payload_format_version: str,
        auth_result_ttl: str,
        authorizer_type: str,
        authorizer_uri: str,
        enable_simple_response: str,
        identity_source: str,
        identity_validation_expr: str,
        jwt_config: str,
        name: str,
    ) -> Authorizer:
        authorizer = Authorizer(
            auth_creds_arn=auth_creds_arn,
            auth_payload_format_version=auth_payload_format_version,
            auth_result_ttl=auth_result_ttl,
            authorizer_type=authorizer_type,
            authorizer_uri=authorizer_uri,
            enable_simple_response=enable_simple_response,
            identity_source=identity_source,
            identity_validation_expr=identity_validation_expr,
            jwt_config=jwt_config,
            name=name,
        )
        self.authorizers[authorizer.id] = authorizer
        return authorizer

    def delete_authorizer(self, authorizer_id: str) -> None:
        self.authorizers.pop(authorizer_id, None)

    def get_authorizer(self, authorizer_id: str) -> Authorizer:
        if authorizer_id not in self.authorizers:
            raise AuthorizerNotFound(authorizer_id)
        return self.authorizers[authorizer_id]

    def update_authorizer(
        self,
        authorizer_id: str,
        auth_creds_arn: str,
        auth_payload_format_version: str,
        auth_result_ttl: str,
        authorizer_type: str,
        authorizer_uri: str,
        enable_simple_response: str,
        identity_source: str,
        identity_validation_expr: str,
        jwt_config: str,
        name: str,
    ) -> Authorizer:
        authorizer = self.authorizers[authorizer_id]
        authorizer.update(
            auth_creds_arn=auth_creds_arn,
            auth_payload_format_version=auth_payload_format_version,
            auth_result_ttl=auth_result_ttl,
            authorizer_type=authorizer_type,
            authorizer_uri=authorizer_uri,
            enable_simple_response=enable_simple_response,
            identity_source=identity_source,
            identity_validation_expr=identity_validation_expr,
            jwt_config=jwt_config,
            name=name,
        )
        return authorizer

    def create_model(
        self, content_type: str, description: str, name: str, schema: str
    ) -> Model:
        model = Model(content_type, description, name, schema)
        self.models[model.id] = model
        return model

    def delete_model(self, model_id: str) -> None:
        self.models.pop(model_id, None)

    def get_model(self, model_id: str) -> Model:
        if model_id not in self.models:
            raise ModelNotFound(model_id)
        return self.models[model_id]

    def update_model(
        self, model_id: str, content_type: str, description: str, name: str, schema: str
    ) -> Model:
        model = self.models[model_id]
        model.update(content_type, description, name, schema)
        return model

    def import_api(self, body_str: str, fail_on_warnings: bool) -> None:
        self.clear()
        body = yaml.safe_load(body_str)
        for path, path_details in body.get("paths", {}).items():
            for method, method_details in path_details.items():
                route_key = f"{method.upper()} {path}"
                for int_type, type_details in method_details.items():
                    if int_type == "responses":
                        for status_code, response_details in type_details.items():
                            content = response_details.get("content", {})
                            for content_type in content.values():
                                for ref in content_type.get("schema", {}).values():
                                    if ref not in self.models and fail_on_warnings:
                                        attr = f"paths.'{path}'({method}).{int_type}.{status_code}.content.schema.{ref}"
                                        raise BadRequestException(
                                            f"Warnings found during import:\n\tParse issue: attribute {attr} is missing"
                                        )
                    if int_type == "x-amazon-apigateway-integration":
                        integration = self.create_integration(
                            connection_type="INTERNET",
                            description="AutoCreate from OpenAPI Import",
                            integration_type=type_details.get("type"),
                            integration_method=type_details.get("httpMethod"),
                            payload_format_version=type_details.get(
                                "payloadFormatVersion"
                            ),
                            integration_uri=type_details.get("uri"),
                        )
                        self.create_route(
                            api_key_required=False,
                            authorization_scopes=[],
                            route_key=route_key,
                            target=f"integrations/{integration.id}",
                        )
        if "title" in body.get("info", {}):
            self.name = body["info"]["title"]
        if "version" in body.get("info", {}):
            self.version = str(body["info"]["version"])
        if "x-amazon-apigateway-cors" in body:
            self.cors_configuration = body["x-amazon-apigateway-cors"]

    def update(
        self,
        api_key_selection_expression: str,
        cors_configuration: str,
        description: str,
        disable_schema_validation: str,
        disable_execute_api_endpoint: str,
        name: str,
        route_selection_expression: str,
        version: str,
    ) -> None:
        if api_key_selection_expression is not None:
            self.api_key_selection_expression = api_key_selection_expression
        if cors_configuration is not None:
            self.cors_configuration = cors_configuration
        if description is not None:
            self.description = description
        if disable_execute_api_endpoint is not None:
            self.disable_execute_api_endpoint = disable_execute_api_endpoint
        if disable_schema_validation is not None:
            self.disable_schema_validation = disable_schema_validation
        if name is not None:
            self.name = name
        if route_selection_expression is not None:
            self.route_selection_expression = route_selection_expression
        if version is not None:
            self.version = version

    def create_integration(
        self,
        connection_type: str,
        description: str,
        integration_method: str,
        integration_type: str,
        integration_uri: str,
        connection_id: Optional[str] = None,
        content_handling_strategy: Optional[str] = None,
        credentials_arn: Optional[str] = None,
        passthrough_behavior: Optional[str] = None,
        payload_format_version: Optional[str] = None,
        integration_subtype: Optional[str] = None,
        request_parameters: Optional[Dict[str, str]] = None,
        request_templates: Optional[Dict[str, str]] = None,
        response_parameters: Optional[Dict[str, Dict[str, str]]] = None,
        template_selection_expression: Optional[str] = None,
        timeout_in_millis: Optional[str] = None,
        tls_config: Optional[Dict[str, str]] = None,
    ) -> Integration:
        integration = Integration(
            connection_id=connection_id,
            connection_type=connection_type,
            content_handling_strategy=content_handling_strategy,
            credentials_arn=credentials_arn,
            description=description,
            integration_method=integration_method,
            integration_type=integration_type,
            integration_uri=integration_uri,
            passthrough_behavior=passthrough_behavior,
            payload_format_version=payload_format_version,
            integration_subtype=integration_subtype,
            request_parameters=request_parameters,
            request_templates=request_templates,
            response_parameters=response_parameters,
            template_selection_expression=template_selection_expression,
            timeout_in_millis=timeout_in_millis,
            tls_config=tls_config,
        )
        self.integrations[integration.id] = integration
        return integration

    def delete_integration(self, integration_id: str) -> None:
        self.integrations.pop(integration_id, None)

    def get_integration(self, integration_id: str) -> Integration:
        if integration_id not in self.integrations:
            raise IntegrationNotFound(integration_id)
        return self.integrations[integration_id]

    def get_integrations(self) -> List[Integration]:
        return list(self.integrations.values())

    def update_integration(
        self,
        integration_id: str,
        connection_id: str,
        connection_type: str,
        content_handling_strategy: str,
        credentials_arn: str,
        description: str,
        integration_method: str,
        integration_type: str,
        integration_uri: str,
        passthrough_behavior: str,
        payload_format_version: str,
        integration_subtype: str,
        request_parameters: Dict[str, str],
        request_templates: Dict[str, str],
        response_parameters: Dict[str, Dict[str, str]],
        template_selection_expression: str,
        timeout_in_millis: Optional[int],
        tls_config: Dict[str, str],
    ) -> Integration:
        integration = self.integrations[integration_id]
        integration.update(
            connection_id=connection_id,
            connection_type=connection_type,
            content_handling_strategy=content_handling_strategy,
            credentials_arn=credentials_arn,
            description=description,
            integration_method=integration_method,
            integration_type=integration_type,
            integration_uri=integration_uri,
            passthrough_behavior=passthrough_behavior,
            payload_format_version=payload_format_version,
            integration_subtype=integration_subtype,
            request_parameters=request_parameters,
            request_templates=request_templates,
            response_parameters=response_parameters,
            template_selection_expression=template_selection_expression,
            timeout_in_millis=timeout_in_millis,
            tls_config=tls_config,
        )
        return integration

    def create_integration_response(
        self,
        integration_id: str,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> IntegrationResponse:
        integration = self.get_integration(integration_id)
        return integration.create_response(
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )

    def delete_integration_response(
        self, integration_id: str, integration_response_id: str
    ) -> None:
        integration = self.get_integration(integration_id)
        integration.delete_response(integration_response_id)

    def get_integration_response(
        self, integration_id: str, integration_response_id: str
    ) -> IntegrationResponse:
        integration = self.get_integration(integration_id)
        return integration.get_response(integration_response_id)

    def get_integration_responses(
        self, integration_id: str
    ) -> List[IntegrationResponse]:
        integration = self.get_integration(integration_id)
        return integration.get_responses()

    def update_integration_response(
        self,
        integration_id: str,
        integration_response_id: str,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> IntegrationResponse:
        integration = self.get_integration(integration_id)
        return integration.update_response(
            integration_response_id=integration_response_id,
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )

    def create_route(
        self,
        api_key_required: bool,
        authorization_scopes: List[str],
        route_key: str,
        target: str,
        authorization_type: Optional[str] = None,
        authorizer_id: Optional[str] = None,
        model_selection_expression: Optional[str] = None,
        operation_name: Optional[str] = None,
        request_models: Optional[Dict[str, str]] = None,
        request_parameters: Optional[Dict[str, Dict[str, bool]]] = None,
        route_response_selection_expression: Optional[str] = None,
    ) -> Route:
        route = Route(
            api_key_required=api_key_required,
            authorization_scopes=authorization_scopes,
            authorization_type=authorization_type,
            authorizer_id=authorizer_id,
            model_selection_expression=model_selection_expression,
            operation_name=operation_name,
            request_models=request_models,
            request_parameters=request_parameters,
            route_key=route_key,
            route_response_selection_expression=route_response_selection_expression,
            target=target,
        )
        self.routes[route.route_id] = route
        return route

    def delete_route(self, route_id: str) -> None:
        self.routes.pop(route_id, None)

    def delete_route_request_parameter(self, route_id: str, request_param: str) -> None:
        route = self.get_route(route_id)
        route.delete_route_request_parameter(request_param)

    def get_route(self, route_id: str) -> Route:
        if route_id not in self.routes:
            raise RouteNotFound(route_id)
        return self.routes[route_id]

    def get_routes(self) -> List[Route]:
        return list(self.routes.values())

    def update_route(
        self,
        route_id: str,
        api_key_required: Optional[bool],
        authorization_scopes: List[str],
        authorization_type: str,
        authorizer_id: str,
        model_selection_expression: str,
        operation_name: str,
        request_models: Dict[str, str],
        request_parameters: Dict[str, Dict[str, bool]],
        route_key: str,
        route_response_selection_expression: str,
        target: str,
    ) -> Route:
        route = self.get_route(route_id)
        route.update(
            api_key_required=api_key_required,
            authorization_scopes=authorization_scopes,
            authorization_type=authorization_type,
            authorizer_id=authorizer_id,
            model_selection_expression=model_selection_expression,
            operation_name=operation_name,
            request_models=request_models,
            request_parameters=request_parameters,
            route_key=route_key,
            route_response_selection_expression=route_response_selection_expression,
            target=target,
        )
        return route

    def create_route_response(
        self,
        route_id: str,
        route_response_key: str,
        model_selection_expression: str,
        response_models: str,
    ) -> RouteResponse:
        route = self.get_route(route_id)
        return route.create_route_response(
            route_response_key,
            model_selection_expression=model_selection_expression,
            response_models=response_models,
        )

    def delete_route_response(self, route_id: str, route_response_id: str) -> None:
        route = self.get_route(route_id)
        route.delete_route_response(route_response_id)

    def get_route_response(
        self, route_id: str, route_response_id: str
    ) -> RouteResponse:
        route = self.get_route(route_id)
        return route.get_route_response(route_response_id)

    def create_stage(self, config: Dict[str, Any]) -> Stage:
        stage = Stage(api=self, config=config)
        self.stages[stage.name] = stage
        return stage

    def get_stage(self, stage_name: str) -> Stage:
        if stage_name not in self.stages:
            raise StageNotFound
        return self.stages[stage_name]

    def delete_stage(self, stage_name: str) -> None:
        self.stages.pop(stage_name, None)

    def to_json(self) -> Dict[str, Any]:
        return {
            "apiId": self.api_id,
            "apiEndpoint": self.api_endpoint,
            "apiKeySelectionExpression": self.api_key_selection_expression,
            "createdDate": self.created_date,
            "corsConfiguration": self.cors_configuration,
            "description": self.description,
            "disableExecuteApiEndpoint": self.disable_execute_api_endpoint,
            "disableSchemaValidation": self.disable_schema_validation,
            "name": self.name,
            "protocolType": self.protocol_type,
            "routeSelectionExpression": self.route_selection_expression,
            "tags": self.backend.get_tags(self.arn),
            "version": self.version,
        }


class VpcLink(BaseModel):
    def __init__(
        self,
        name: str,
        sg_ids: List[str],
        subnet_ids: List[str],
        tags: Dict[str, str],
        backend: "ApiGatewayV2Backend",
    ):
        self.created = unix_time()
        self.id = "".join(random.choice(string.ascii_lowercase) for _ in range(8))
        self.name = name
        self.sg_ids = sg_ids
        self.subnet_ids = subnet_ids

        self.arn = f"arn:{get_partition(backend.region_name)}:apigateway:{backend.region_name}::/vpclinks/{self.id}"
        self.backend = backend
        self.backend.tag_resource(self.arn, tags)

    def update(self, name: str) -> None:
        self.name = name

    def to_json(self) -> Dict[str, Any]:
        return {
            "createdDate": self.created,
            "name": self.name,
            "securityGroupIds": self.sg_ids,
            "subnetIds": self.subnet_ids,
            "tags": self.backend.get_tags(self.arn),
            "vpcLinkId": self.id,
            "vpcLinkStatus": "AVAILABLE",
            "vpcLinkVersion": "V2",
        }


class DomainName(BaseModel):
    def __init__(
        self,
        domain_name: str,
        domain_name_configurations: List[Dict[str, str]],
        mutual_tls_authentication: Dict[str, str],
        tags: Dict[str, str],
    ):
        self.api_mapping_selection_expression = "$request.basepath"
        self.domain_name = domain_name
        self.domain_name_configurations = domain_name_configurations
        self.mutual_tls_authentication = mutual_tls_authentication
        self.tags = tags

    def to_json(self) -> Dict[str, Any]:
        return {
            "apiMappingSelectionExpression": self.api_mapping_selection_expression,
            "domainName": self.domain_name,
            "domainNameConfigurations": self.domain_name_configurations,
            "mutualTlsAuthentication": self.mutual_tls_authentication,
            "tags": self.tags,
        }


class ApiMapping(BaseModel):
    def __init__(
        self,
        api_id: str,
        api_mapping_key: str,
        api_mapping_id: str,
        domain_name: str,
        stage: str,
    ) -> None:
        self.api_id = api_id
        self.api_mapping_key = api_mapping_key
        self.api_mapping_id = api_mapping_id
        self.domain_name = domain_name
        self.stage = stage

    def to_json(self) -> Dict[str, Any]:
        return {
            "apiId": self.api_id,
            "apiMappingId": self.api_mapping_id,
            "apiMappingKey": self.api_mapping_key,
            "domainName": self.domain_name,
            "stage": self.stage,
        }


class ApiGatewayV2Backend(BaseBackend):
    """Implementation of ApiGatewayV2 APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.apis: Dict[str, Api] = dict()
        self.vpc_links: Dict[str, VpcLink] = dict()
        self.domain_names: Dict[str, DomainName] = dict()
        self.api_mappings: Dict[str, ApiMapping] = dict()
        self.tagger = TaggingService()

    def create_api(
        self,
        api_key_selection_expression: str,
        cors_configuration: str,
        description: str,
        disable_schema_validation: str,
        disable_execute_api_endpoint: str,
        name: str,
        protocol_type: str,
        route_selection_expression: str,
        tags: Dict[str, str],
        version: str,
    ) -> Api:
        """
        The following parameters are not yet implemented:
        CredentialsArn, RouteKey, Tags, Target
        """
        api = Api(
            region=self.region_name,
            cors_configuration=cors_configuration,
            description=description,
            name=name,
            api_key_selection_expression=api_key_selection_expression,
            disable_execute_api_endpoint=disable_execute_api_endpoint,
            disable_schema_validation=disable_schema_validation,
            protocol_type=protocol_type,
            route_selection_expression=route_selection_expression,
            tags=tags,
            version=version,
            backend=self,
        )
        self.apis[api.api_id] = api
        return api

    def delete_api(self, api_id: str) -> None:
        self.apis.pop(api_id, None)

    def get_api(self, api_id: str) -> Api:
        if api_id not in self.apis:
            raise ApiNotFound(api_id)
        return self.apis[api_id]

    def get_apis(self) -> List[Api]:
        """
        Pagination is not yet implemented
        """
        return list(self.apis.values())

    def update_api(
        self,
        api_id: str,
        api_key_selection_expression: str,
        cors_configuration: str,
        description: str,
        disable_schema_validation: str,
        disable_execute_api_endpoint: str,
        name: str,
        route_selection_expression: str,
        version: str,
    ) -> Api:
        """
        The following parameters have not yet been implemented: CredentialsArn, RouteKey, Target
        """
        api = self.get_api(api_id)
        api.update(
            api_key_selection_expression=api_key_selection_expression,
            cors_configuration=cors_configuration,
            description=description,
            disable_schema_validation=disable_schema_validation,
            disable_execute_api_endpoint=disable_execute_api_endpoint,
            name=name,
            route_selection_expression=route_selection_expression,
            version=version,
        )
        return api

    def reimport_api(self, api_id: str, body: str, fail_on_warnings: bool) -> Api:
        """
        Only YAML is supported at the moment. Full OpenAPI-support is not guaranteed. Only limited validation is implemented
        """
        api = self.get_api(api_id)
        api.import_api(body, fail_on_warnings)
        return api

    def delete_cors_configuration(self, api_id: str) -> None:
        api = self.get_api(api_id)
        api.delete_cors_configuration()

    def create_authorizer(
        self,
        api_id: str,
        auth_creds_arn: str,
        auth_payload_format_version: str,
        auth_result_ttl: str,
        authorizer_uri: str,
        authorizer_type: str,
        enable_simple_response: str,
        identity_source: str,
        identity_validation_expr: str,
        jwt_config: str,
        name: str,
    ) -> Authorizer:
        api = self.get_api(api_id)

        if (
            api.protocol_type == "HTTP"
            and authorizer_type == "REQUEST"
            and not auth_payload_format_version
        ):
            raise BadRequestException(
                "AuthorizerPayloadFormatVersion is a required parameter for REQUEST authorizer"
            )

        authorizer = api.create_authorizer(
            auth_creds_arn=auth_creds_arn,
            auth_payload_format_version=auth_payload_format_version,
            auth_result_ttl=auth_result_ttl,
            authorizer_type=authorizer_type,
            authorizer_uri=authorizer_uri,
            enable_simple_response=enable_simple_response,
            identity_source=identity_source,
            identity_validation_expr=identity_validation_expr,
            jwt_config=jwt_config,
            name=name,
        )
        return authorizer

    def delete_authorizer(self, api_id: str, authorizer_id: str) -> None:
        api = self.get_api(api_id)
        api.delete_authorizer(authorizer_id=authorizer_id)

    def get_authorizer(self, api_id: str, authorizer_id: str) -> Authorizer:
        api = self.get_api(api_id)
        authorizer = api.get_authorizer(authorizer_id=authorizer_id)
        return authorizer

    def update_authorizer(
        self,
        api_id: str,
        authorizer_id: str,
        auth_creds_arn: str,
        auth_payload_format_version: str,
        auth_result_ttl: str,
        authorizer_uri: str,
        authorizer_type: str,
        enable_simple_response: str,
        identity_source: str,
        identity_validation_expr: str,
        jwt_config: str,
        name: str,
    ) -> Authorizer:
        api = self.get_api(api_id)
        authorizer = api.update_authorizer(
            authorizer_id=authorizer_id,
            auth_creds_arn=auth_creds_arn,
            auth_payload_format_version=auth_payload_format_version,
            auth_result_ttl=auth_result_ttl,
            authorizer_type=authorizer_type,
            authorizer_uri=authorizer_uri,
            enable_simple_response=enable_simple_response,
            identity_source=identity_source,
            identity_validation_expr=identity_validation_expr,
            jwt_config=jwt_config,
            name=name,
        )
        return authorizer

    def create_model(
        self, api_id: str, content_type: str, description: str, name: str, schema: str
    ) -> Model:
        api = self.get_api(api_id)
        model = api.create_model(
            content_type=content_type, description=description, name=name, schema=schema
        )
        return model

    def delete_model(self, api_id: str, model_id: str) -> None:
        api = self.get_api(api_id)
        api.delete_model(model_id=model_id)

    def get_model(self, api_id: str, model_id: str) -> Model:
        api = self.get_api(api_id)
        return api.get_model(model_id)

    def update_model(
        self,
        api_id: str,
        model_id: str,
        content_type: str,
        description: str,
        name: str,
        schema: str,
    ) -> Model:
        api = self.get_api(api_id)
        return api.update_model(model_id, content_type, description, name, schema)

    def get_tags(self, resource_id: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(resource_id)

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tags_input = TaggingService.convert_dict_to_tags_input(tags or {})
        self.tagger.tag_resource(resource_arn, tags_input)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def create_route(
        self,
        api_id: str,
        api_key_required: bool,
        authorization_scopes: List[str],
        authorization_type: str,
        authorizer_id: str,
        model_selection_expression: str,
        operation_name: str,
        request_models: Optional[Dict[str, str]],
        request_parameters: Optional[Dict[str, Dict[str, bool]]],
        route_key: str,
        route_response_selection_expression: str,
        target: str,
    ) -> Route:
        api = self.get_api(api_id)
        route = api.create_route(
            api_key_required=api_key_required,
            authorization_scopes=authorization_scopes,
            authorization_type=authorization_type,
            authorizer_id=authorizer_id,
            model_selection_expression=model_selection_expression,
            operation_name=operation_name,
            request_models=request_models,
            request_parameters=request_parameters,
            route_key=route_key,
            route_response_selection_expression=route_response_selection_expression,
            target=target,
        )
        return route

    def delete_route(self, api_id: str, route_id: str) -> None:
        api = self.get_api(api_id)
        api.delete_route(route_id)

    def delete_route_request_parameter(
        self, api_id: str, route_id: str, request_param: str
    ) -> None:
        api = self.get_api(api_id)
        api.delete_route_request_parameter(route_id, request_param)

    def get_route(self, api_id: str, route_id: str) -> Route:
        api = self.get_api(api_id)
        return api.get_route(route_id)

    def get_routes(self, api_id: str) -> List[Route]:
        """
        Pagination is not yet implemented
        """
        api = self.get_api(api_id)
        return api.get_routes()

    def update_route(
        self,
        api_id: str,
        api_key_required: bool,
        authorization_scopes: List[str],
        authorization_type: str,
        authorizer_id: str,
        model_selection_expression: str,
        operation_name: str,
        request_models: Dict[str, str],
        request_parameters: Dict[str, Dict[str, bool]],
        route_id: str,
        route_key: str,
        route_response_selection_expression: str,
        target: str,
    ) -> Route:
        api = self.get_api(api_id)
        route = api.update_route(
            route_id=route_id,
            api_key_required=api_key_required,
            authorization_scopes=authorization_scopes,
            authorization_type=authorization_type,
            authorizer_id=authorizer_id,
            model_selection_expression=model_selection_expression,
            operation_name=operation_name,
            request_models=request_models,
            request_parameters=request_parameters,
            route_key=route_key,
            route_response_selection_expression=route_response_selection_expression,
            target=target,
        )
        return route

    def create_route_response(
        self,
        api_id: str,
        route_id: str,
        route_response_key: str,
        model_selection_expression: str,
        response_models: str,
    ) -> RouteResponse:
        """
        The following parameters are not yet implemented: ResponseModels, ResponseParameters
        """
        api = self.get_api(api_id)
        return api.create_route_response(
            route_id,
            route_response_key,
            model_selection_expression=model_selection_expression,
            response_models=response_models,
        )

    def delete_route_response(
        self, api_id: str, route_id: str, route_response_id: str
    ) -> None:
        api = self.get_api(api_id)
        api.delete_route_response(route_id, route_response_id)

    def get_route_response(
        self, api_id: str, route_id: str, route_response_id: str
    ) -> RouteResponse:
        api = self.get_api(api_id)
        return api.get_route_response(route_id, route_response_id)

    def create_integration(
        self,
        api_id: str,
        connection_id: str,
        connection_type: str,
        content_handling_strategy: str,
        credentials_arn: str,
        description: str,
        integration_method: str,
        integration_subtype: str,
        integration_type: str,
        integration_uri: str,
        passthrough_behavior: str,
        payload_format_version: str,
        request_parameters: Optional[Dict[str, str]],
        request_templates: Optional[Dict[str, str]],
        response_parameters: Optional[Dict[str, Dict[str, str]]],
        template_selection_expression: str,
        timeout_in_millis: str,
        tls_config: Dict[str, str],
    ) -> Integration:
        api = self.get_api(api_id)
        integration = api.create_integration(
            connection_id=connection_id,
            connection_type=connection_type,
            content_handling_strategy=content_handling_strategy,
            credentials_arn=credentials_arn,
            description=description,
            integration_method=integration_method,
            integration_type=integration_type,
            integration_uri=integration_uri,
            passthrough_behavior=passthrough_behavior,
            payload_format_version=payload_format_version,
            integration_subtype=integration_subtype,
            request_parameters=request_parameters,
            request_templates=request_templates,
            response_parameters=response_parameters,
            template_selection_expression=template_selection_expression,
            timeout_in_millis=timeout_in_millis,
            tls_config=tls_config,
        )
        return integration

    def get_integration(self, api_id: str, integration_id: str) -> Integration:
        api = self.get_api(api_id)
        integration = api.get_integration(integration_id)
        return integration

    def get_integrations(self, api_id: str) -> List[Integration]:
        """
        Pagination is not yet implemented
        """
        api = self.get_api(api_id)
        return api.get_integrations()

    def delete_integration(self, api_id: str, integration_id: str) -> None:
        api = self.get_api(api_id)
        api.delete_integration(integration_id)

    def update_integration(
        self,
        api_id: str,
        connection_id: str,
        connection_type: str,
        content_handling_strategy: str,
        credentials_arn: str,
        description: str,
        integration_id: str,
        integration_method: str,
        integration_subtype: str,
        integration_type: str,
        integration_uri: str,
        passthrough_behavior: str,
        payload_format_version: str,
        request_parameters: Dict[str, str],
        request_templates: Dict[str, str],
        response_parameters: Dict[str, Dict[str, str]],
        template_selection_expression: str,
        timeout_in_millis: Optional[int],
        tls_config: Dict[str, str],
    ) -> Integration:
        api = self.get_api(api_id)
        integration = api.update_integration(
            integration_id=integration_id,
            connection_id=connection_id,
            connection_type=connection_type,
            content_handling_strategy=content_handling_strategy,
            credentials_arn=credentials_arn,
            description=description,
            integration_method=integration_method,
            integration_type=integration_type,
            integration_uri=integration_uri,
            passthrough_behavior=passthrough_behavior,
            payload_format_version=payload_format_version,
            integration_subtype=integration_subtype,
            request_parameters=request_parameters,
            request_templates=request_templates,
            response_parameters=response_parameters,
            template_selection_expression=template_selection_expression,
            timeout_in_millis=timeout_in_millis,
            tls_config=tls_config,
        )
        return integration

    def create_integration_response(
        self,
        api_id: str,
        integration_id: str,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> IntegrationResponse:
        api = self.get_api(api_id)
        integration_response = api.create_integration_response(
            integration_id=integration_id,
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )
        return integration_response

    def delete_integration_response(
        self, api_id: str, integration_id: str, integration_response_id: str
    ) -> None:
        api = self.get_api(api_id)
        api.delete_integration_response(
            integration_id, integration_response_id=integration_response_id
        )

    def get_integration_response(
        self, api_id: str, integration_id: str, integration_response_id: str
    ) -> IntegrationResponse:
        api = self.get_api(api_id)
        return api.get_integration_response(
            integration_id, integration_response_id=integration_response_id
        )

    def get_integration_responses(
        self, api_id: str, integration_id: str
    ) -> List[IntegrationResponse]:
        api = self.get_api(api_id)
        return api.get_integration_responses(integration_id)

    def update_integration_response(
        self,
        api_id: str,
        integration_id: str,
        integration_response_id: str,
        content_handling_strategy: str,
        integration_response_key: str,
        response_parameters: str,
        response_templates: str,
        template_selection_expression: str,
    ) -> IntegrationResponse:
        api = self.get_api(api_id)
        integration_response = api.update_integration_response(
            integration_id=integration_id,
            integration_response_id=integration_response_id,
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )
        return integration_response

    def create_vpc_link(
        self, name: str, sg_ids: List[str], subnet_ids: List[str], tags: Dict[str, str]
    ) -> VpcLink:
        vpc_link = VpcLink(
            name, sg_ids=sg_ids, subnet_ids=subnet_ids, tags=tags, backend=self
        )
        self.vpc_links[vpc_link.id] = vpc_link
        return vpc_link

    def get_vpc_link(self, vpc_link_id: str) -> VpcLink:
        if vpc_link_id not in self.vpc_links:
            raise VpcLinkNotFound(vpc_link_id)
        return self.vpc_links[vpc_link_id]

    def delete_vpc_link(self, vpc_link_id: str) -> None:
        self.vpc_links.pop(vpc_link_id, None)

    def get_vpc_links(self) -> List[VpcLink]:
        return list(self.vpc_links.values())

    def update_vpc_link(self, vpc_link_id: str, name: str) -> VpcLink:
        vpc_link = self.get_vpc_link(vpc_link_id)
        vpc_link.update(name)
        return vpc_link

    def create_domain_name(
        self,
        domain_name: str,
        domain_name_configurations: List[Dict[str, str]],
        mutual_tls_authentication: Dict[str, str],
        tags: Dict[str, str],
    ) -> DomainName:
        if domain_name in self.domain_names.keys():
            raise DomainNameAlreadyExists

        domain = DomainName(
            domain_name=domain_name,
            domain_name_configurations=domain_name_configurations,
            mutual_tls_authentication=mutual_tls_authentication,
            tags=tags,
        )
        self.domain_names[domain.domain_name] = domain
        return domain

    def get_domain_name(self, domain_name: Union[str, None]) -> DomainName:
        if domain_name is None or domain_name not in self.domain_names:
            raise DomainNameNotFound
        return self.domain_names[domain_name]

    def get_domain_names(self) -> List[DomainName]:
        """
        Pagination is not yet implemented
        """
        return list(self.domain_names.values())

    def delete_domain_name(self, domain_name: str) -> None:
        if domain_name not in self.domain_names.keys():
            raise DomainNameNotFound

        for mapping_id, mapping in self.api_mappings.items():
            if mapping.domain_name == domain_name:
                del self.api_mappings[mapping_id]

        del self.domain_names[domain_name]

    def _generate_api_maping_id(
        self, api_mapping_key: str, stage: str, domain_name: str
    ) -> str:
        return str(
            hashlib.sha256(
                f"{stage} {domain_name}/{api_mapping_key}".encode("utf-8")
            ).hexdigest()
        )[:5]

    def create_api_mapping(
        self, api_id: str, api_mapping_key: str, domain_name: str, stage: str
    ) -> ApiMapping:
        if domain_name not in self.domain_names.keys():
            raise DomainNameNotFound

        if api_id not in self.apis.keys():
            raise ApiNotFound("The resource specified in the request was not found.")

        if api_mapping_key.startswith("/") or "//" in api_mapping_key:
            raise BadRequestException(
                "API mapping key should not start with a '/' or have consecutive '/'s."
            )

        if api_mapping_key.endswith("/"):
            raise BadRequestException("API mapping key should not end with a '/'.")

        api_mapping_id = self._generate_api_maping_id(
            api_mapping_key=api_mapping_key, stage=stage, domain_name=domain_name
        )

        mapping = ApiMapping(
            domain_name=domain_name,
            api_id=api_id,
            api_mapping_key=api_mapping_key,
            api_mapping_id=api_mapping_id,
            stage=stage,
        )

        self.api_mappings[api_mapping_id] = mapping
        return mapping

    def get_api_mapping(self, api_mapping_id: str, domain_name: str) -> ApiMapping:
        if domain_name not in self.domain_names.keys():
            raise DomainNameNotFound

        if api_mapping_id not in self.api_mappings.keys():
            raise ApiMappingNotFound

        return self.api_mappings[api_mapping_id]

    def get_api_mappings(self, domain_name: str) -> List[ApiMapping]:
        domain_mappings = []
        for mapping in self.api_mappings.values():
            if mapping.domain_name == domain_name:
                domain_mappings.append(mapping)
        return domain_mappings

    def delete_api_mapping(self, api_mapping_id: str, domain_name: str) -> None:
        if api_mapping_id not in self.api_mappings.keys():
            raise ApiMappingNotFound

        if self.api_mappings[api_mapping_id].domain_name != domain_name:
            raise BadRequestException(
                f"given domain name {domain_name} does not match with mapping definition of mapping {api_mapping_id}"
            )

        del self.api_mappings[api_mapping_id]

    def create_stage(self, api_id: str, config: Dict[str, Any]) -> Stage:
        api = self.get_api(api_id)
        return api.create_stage(config)

    def get_stage(self, api_id: str, stage_name: str) -> Stage:
        api = self.get_api(api_id)
        return api.get_stage(stage_name)

    def delete_stage(self, api_id: str, stage_name: str) -> None:
        api = self.get_api(api_id)
        api.delete_stage(stage_name)

    def get_stages(self, api_id: str) -> List[Stage]:
        api = self.get_api(api_id)
        return list(api.stages.values())


apigatewayv2_backends = BackendDict(ApiGatewayV2Backend, "apigatewayv2")
