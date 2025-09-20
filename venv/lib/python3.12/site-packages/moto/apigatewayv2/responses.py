"""Handles incoming apigatewayv2 requests, invokes methods, returns responses."""

import json
from typing import Any
from urllib.parse import unquote

from moto.core.responses import TYPE_RESPONSE, BaseResponse

from .exceptions import UnknownProtocol
from .models import ApiGatewayV2Backend, apigatewayv2_backends


class ApiGatewayV2Response(BaseResponse):
    """Handler for ApiGatewayV2 requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="apigatewayv2")

    @property
    def apigatewayv2_backend(self) -> ApiGatewayV2Backend:
        """Return backend instance specific for this region."""
        return apigatewayv2_backends[self.current_account][self.region]

    def create_api(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)

        api_key_selection_expression = params.get("apiKeySelectionExpression")
        cors_configuration = params.get("corsConfiguration")
        description = params.get("description")
        disable_schema_validation = params.get("disableSchemaValidation")
        disable_execute_api_endpoint = params.get("disableExecuteApiEndpoint")
        name = params.get("name")
        protocol_type = params.get("protocolType")
        route_selection_expression = params.get("routeSelectionExpression")
        tags = params.get("tags")
        version = params.get("version")

        if protocol_type not in ["HTTP", "WEBSOCKET"]:
            raise UnknownProtocol

        api = self.apigatewayv2_backend.create_api(
            api_key_selection_expression=api_key_selection_expression,
            cors_configuration=cors_configuration,
            description=description,
            disable_schema_validation=disable_schema_validation,
            disable_execute_api_endpoint=disable_execute_api_endpoint,
            name=name,
            protocol_type=protocol_type,
            route_selection_expression=route_selection_expression,
            tags=tags,
            version=version,
        )
        return 200, {}, json.dumps(api.to_json())

    def delete_api(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-1]
        self.apigatewayv2_backend.delete_api(api_id=api_id)
        return 200, {}, "{}"

    def get_api(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-1]
        api = self.apigatewayv2_backend.get_api(api_id=api_id)
        return 200, {}, json.dumps(api.to_json())

    @staticmethod
    def get_api_without_id(*args: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        """
        AWS is returning an empty response when apiId is an empty string. This is slightly odd and it seems an
        outlier, therefore it was decided we could have a custom handler for this particular use case instead of
        trying to make it work with the existing url-matcher.
        """
        return 200, {}, "{}"

    def get_apis(self) -> TYPE_RESPONSE:
        apis = self.apigatewayv2_backend.get_apis()
        return 200, {}, json.dumps({"items": [a.to_json() for a in apis]})

    def update_api(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-1]
        params = json.loads(self.body)
        api_key_selection_expression = params.get("apiKeySelectionExpression")
        cors_configuration = params.get("corsConfiguration")
        description = params.get("description")
        disable_schema_validation = params.get("disableSchemaValidation")
        disable_execute_api_endpoint = params.get("disableExecuteApiEndpoint")
        name = params.get("name")
        route_selection_expression = params.get("routeSelectionExpression")
        version = params.get("version")
        api = self.apigatewayv2_backend.update_api(
            api_id=api_id,
            api_key_selection_expression=api_key_selection_expression,
            cors_configuration=cors_configuration,
            description=description,
            disable_schema_validation=disable_schema_validation,
            disable_execute_api_endpoint=disable_execute_api_endpoint,
            name=name,
            route_selection_expression=route_selection_expression,
            version=version,
        )
        return 200, {}, json.dumps(api.to_json())

    def reimport_api(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-1]
        params = json.loads(self.body)
        body = params.get("body")
        fail_on_warnings = (
            str(self._get_param("failOnWarnings", "false")).lower() == "true"
        )

        api = self.apigatewayv2_backend.reimport_api(api_id, body, fail_on_warnings)
        return 201, {}, json.dumps(api.to_json())

    def create_authorizer(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        params = json.loads(self.body)

        auth_creds_arn = params.get("authorizerCredentialsArn")
        auth_payload_format_version = params.get("authorizerPayloadFormatVersion")
        auth_result_ttl = params.get("authorizerResultTtlInSeconds")
        authorizer_type = params.get("authorizerType")
        authorizer_uri = params.get("authorizerUri")
        enable_simple_response = params.get("enableSimpleResponses")
        identity_source = params.get("identitySource")
        identity_validation_expr = params.get("identityValidationExpression")
        jwt_config = params.get("jwtConfiguration")
        name = params.get("name")
        authorizer = self.apigatewayv2_backend.create_authorizer(
            api_id,
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
        return 200, {}, json.dumps(authorizer.to_json())

    def delete_authorizer(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        authorizer_id = self.path.split("/")[-1]

        self.apigatewayv2_backend.delete_authorizer(api_id, authorizer_id)
        return 200, {}, "{}"

    def get_authorizer(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        authorizer_id = self.path.split("/")[-1]

        authorizer = self.apigatewayv2_backend.get_authorizer(api_id, authorizer_id)
        return 200, {}, json.dumps(authorizer.to_json())

    def update_authorizer(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        authorizer_id = self.path.split("/")[-1]
        params = json.loads(self.body)

        auth_creds_arn = params.get("authorizerCredentialsArn")
        auth_payload_format_version = params.get("authorizerPayloadFormatVersion")
        auth_result_ttl = params.get("authorizerResultTtlInSeconds")
        authorizer_type = params.get("authorizerType")
        authorizer_uri = params.get("authorizerUri")
        enable_simple_response = params.get("enableSimpleResponses")
        identity_source = params.get("identitySource")
        identity_validation_expr = params.get("identityValidationExpression")
        jwt_config = params.get("jwtConfiguration")
        name = params.get("name")
        authorizer = self.apigatewayv2_backend.update_authorizer(
            api_id,
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
        return 200, {}, json.dumps(authorizer.to_json())

    def delete_cors_configuration(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        self.apigatewayv2_backend.delete_cors_configuration(api_id)
        return 200, {}, "{}"

    def create_model(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        params = json.loads(self.body)

        content_type = params.get("contentType")
        description = params.get("description")
        name = params.get("name")
        schema = params.get("schema")
        model = self.apigatewayv2_backend.create_model(
            api_id, content_type, description, name, schema
        )
        return 200, {}, json.dumps(model.to_json())

    def delete_model(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        model_id = self.path.split("/")[-1]

        self.apigatewayv2_backend.delete_model(api_id, model_id)
        return 200, {}, "{}"

    def get_model(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        model_id = self.path.split("/")[-1]

        model = self.apigatewayv2_backend.get_model(api_id, model_id)
        return 200, {}, json.dumps(model.to_json())

    def update_model(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        model_id = self.path.split("/")[-1]
        params = json.loads(self.body)

        content_type = params.get("contentType")
        description = params.get("description")
        name = params.get("name")
        schema = params.get("schema")

        model = self.apigatewayv2_backend.update_model(
            api_id,
            model_id,
            content_type=content_type,
            description=description,
            name=name,
            schema=schema,
        )
        return 200, {}, json.dumps(model.to_json())

    def get_tags(self) -> TYPE_RESPONSE:
        resource_arn = unquote(self.path.split("/tags/")[1])
        tags = self.apigatewayv2_backend.get_tags(resource_arn)
        return 200, {}, json.dumps({"tags": tags})

    def tag_resource(self) -> TYPE_RESPONSE:
        resource_arn = unquote(self.path.split("/tags/")[1])
        tags = json.loads(self.body).get("tags", {})
        self.apigatewayv2_backend.tag_resource(resource_arn, tags)
        return 201, {}, "{}"

    def untag_resource(self) -> TYPE_RESPONSE:
        resource_arn = unquote(self.path.split("/tags/")[1])
        tag_keys = self.querystring.get("tagKeys") or []
        self.apigatewayv2_backend.untag_resource(resource_arn, tag_keys)
        return 200, {}, "{}"

    def create_route(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        params = json.loads(self.body)
        api_key_required: bool = params.get("apiKeyRequired", False)
        authorization_scopes = params.get("authorizationScopes")
        authorization_type = params.get("authorizationType", "NONE")
        authorizer_id = params.get("authorizerId")
        model_selection_expression = params.get("modelSelectionExpression")
        operation_name = params.get("operationName")
        request_models = params.get("requestModels")
        request_parameters = params.get("requestParameters")
        route_key = params.get("routeKey")
        route_response_selection_expression = params.get(
            "routeResponseSelectionExpression"
        )
        target = params.get("target")
        route = self.apigatewayv2_backend.create_route(
            api_id=api_id,
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
        return 201, {}, json.dumps(route.to_json())

    def delete_route(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        route_id = self.path.split("/")[-1]
        self.apigatewayv2_backend.delete_route(api_id=api_id, route_id=route_id)
        return 200, {}, "{}"

    def delete_route_request_parameter(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-5]
        route_id = self.path.split("/")[-3]
        request_param = self.path.split("/")[-1]
        self.apigatewayv2_backend.delete_route_request_parameter(
            api_id, route_id, request_param
        )
        return 200, {}, "{}"

    def get_route(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        route_id = self.path.split("/")[-1]
        api = self.apigatewayv2_backend.get_route(api_id=api_id, route_id=route_id)
        return 200, {}, json.dumps(api.to_json())

    def get_routes(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        apis = self.apigatewayv2_backend.get_routes(api_id=api_id)
        return 200, {}, json.dumps({"items": [api.to_json() for api in apis]})

    def update_route(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        route_id = self.path.split("/")[-1]

        params = json.loads(self.body)
        api_key_required = params.get("apiKeyRequired")
        authorization_scopes = params.get("authorizationScopes")
        authorization_type = params.get("authorizationType")
        authorizer_id = params.get("authorizerId")
        model_selection_expression = params.get("modelSelectionExpression")
        operation_name = params.get("operationName")
        request_models = params.get("requestModels")
        request_parameters = params.get("requestParameters")
        route_key = params.get("routeKey")
        route_response_selection_expression = params.get(
            "routeResponseSelectionExpression"
        )
        target = params.get("target")
        api = self.apigatewayv2_backend.update_route(
            api_id=api_id,
            api_key_required=api_key_required,
            authorization_scopes=authorization_scopes,
            authorization_type=authorization_type,
            authorizer_id=authorizer_id,
            model_selection_expression=model_selection_expression,
            operation_name=operation_name,
            request_models=request_models,
            request_parameters=request_parameters,
            route_id=route_id,
            route_key=route_key,
            route_response_selection_expression=route_response_selection_expression,
            target=target,
        )
        return 200, {}, json.dumps(api.to_json())

    def create_route_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-4]
        route_id = self.path.split("/")[-2]
        params = json.loads(self.body)

        response_models = params.get("responseModels")
        route_response_key = params.get("routeResponseKey")
        model_selection_expression = params.get("modelSelectionExpression")
        route_response = self.apigatewayv2_backend.create_route_response(
            api_id,
            route_id,
            route_response_key,
            model_selection_expression=model_selection_expression,
            response_models=response_models,
        )
        return 200, {}, json.dumps(route_response.to_json())

    def delete_route_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-5]
        route_id = self.path.split("/")[-3]
        route_response_id = self.path.split("/")[-1]

        self.apigatewayv2_backend.delete_route_response(
            api_id, route_id, route_response_id
        )
        return 200, {}, "{}"

    def get_route_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-5]
        route_id = self.path.split("/")[-3]
        route_response_id = self.path.split("/")[-1]

        route_response = self.apigatewayv2_backend.get_route_response(
            api_id, route_id, route_response_id
        )
        return 200, {}, json.dumps(route_response.to_json())

    def create_integration(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]

        params = json.loads(self.body)
        connection_id = params.get("connectionId")
        connection_type = params.get("connectionType")
        content_handling_strategy = params.get("contentHandlingStrategy")
        credentials_arn = params.get("credentialsArn")
        description = params.get("description")
        integration_method = params.get("integrationMethod")
        integration_subtype = params.get("integrationSubtype")
        integration_type = params.get("integrationType")
        integration_uri = params.get("integrationUri")
        passthrough_behavior = params.get("passthroughBehavior")
        payload_format_version = params.get("payloadFormatVersion")
        request_parameters = params.get("requestParameters")
        request_templates = params.get("requestTemplates")
        response_parameters = params.get("responseParameters")
        template_selection_expression = params.get("templateSelectionExpression")
        timeout_in_millis = params.get("timeoutInMillis")
        tls_config = params.get("tlsConfig")
        integration = self.apigatewayv2_backend.create_integration(
            api_id=api_id,
            connection_id=connection_id,
            connection_type=connection_type,
            content_handling_strategy=content_handling_strategy,
            credentials_arn=credentials_arn,
            description=description,
            integration_method=integration_method,
            integration_subtype=integration_subtype,
            integration_type=integration_type,
            integration_uri=integration_uri,
            passthrough_behavior=passthrough_behavior,
            payload_format_version=payload_format_version,
            request_parameters=request_parameters,
            request_templates=request_templates,
            response_parameters=response_parameters,
            template_selection_expression=template_selection_expression,
            timeout_in_millis=timeout_in_millis,
            tls_config=tls_config,
        )
        return 200, {}, json.dumps(integration.to_json())

    def get_integration(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        integration_id = self.path.split("/")[-1]

        integration = self.apigatewayv2_backend.get_integration(
            api_id=api_id, integration_id=integration_id
        )
        return 200, {}, json.dumps(integration.to_json())

    def get_integrations(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]

        integrations = self.apigatewayv2_backend.get_integrations(api_id=api_id)
        return 200, {}, json.dumps({"items": [i.to_json() for i in integrations]})

    def delete_integration(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        integration_id = self.path.split("/")[-1]

        self.apigatewayv2_backend.delete_integration(
            api_id=api_id, integration_id=integration_id
        )
        return 200, {}, "{}"

    def update_integration(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        integration_id = self.path.split("/")[-1]

        params = json.loads(self.body)
        connection_id = params.get("connectionId")
        connection_type = params.get("connectionType")
        content_handling_strategy = params.get("contentHandlingStrategy")
        credentials_arn = params.get("credentialsArn")
        description = params.get("description")
        integration_method = params.get("integrationMethod")
        integration_subtype = params.get("integrationSubtype")
        integration_type = params.get("integrationType")
        integration_uri = params.get("integrationUri")
        passthrough_behavior = params.get("passthroughBehavior")
        payload_format_version = params.get("payloadFormatVersion")
        request_parameters = params.get("requestParameters")
        request_templates = params.get("requestTemplates")
        response_parameters = params.get("responseParameters")
        template_selection_expression = params.get("templateSelectionExpression")
        timeout_in_millis = params.get("timeoutInMillis")
        tls_config = params.get("tlsConfig")
        integration = self.apigatewayv2_backend.update_integration(
            api_id=api_id,
            connection_id=connection_id,
            connection_type=connection_type,
            content_handling_strategy=content_handling_strategy,
            credentials_arn=credentials_arn,
            description=description,
            integration_id=integration_id,
            integration_method=integration_method,
            integration_subtype=integration_subtype,
            integration_type=integration_type,
            integration_uri=integration_uri,
            passthrough_behavior=passthrough_behavior,
            payload_format_version=payload_format_version,
            request_parameters=request_parameters,
            request_templates=request_templates,
            response_parameters=response_parameters,
            template_selection_expression=template_selection_expression,
            timeout_in_millis=timeout_in_millis,
            tls_config=tls_config,
        )
        return 200, {}, json.dumps(integration.to_json())

    def create_integration_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-4]
        int_id = self.path.split("/")[-2]

        params = json.loads(self.body)
        content_handling_strategy = params.get("contentHandlingStrategy")
        integration_response_key = params.get("integrationResponseKey")
        response_parameters = params.get("responseParameters")
        response_templates = params.get("responseTemplates")
        template_selection_expression = params.get("templateSelectionExpression")
        integration_response = self.apigatewayv2_backend.create_integration_response(
            api_id=api_id,
            integration_id=int_id,
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )
        return 200, {}, json.dumps(integration_response.to_json())

    def delete_integration_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-5]
        int_id = self.path.split("/")[-3]
        int_res_id = self.path.split("/")[-1]

        self.apigatewayv2_backend.delete_integration_response(
            api_id=api_id, integration_id=int_id, integration_response_id=int_res_id
        )
        return 200, {}, "{}"

    def get_integration_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-5]
        int_id = self.path.split("/")[-3]
        int_res_id = self.path.split("/")[-1]

        int_response = self.apigatewayv2_backend.get_integration_response(
            api_id=api_id, integration_id=int_id, integration_response_id=int_res_id
        )
        return 200, {}, json.dumps(int_response.to_json())

    def get_integration_responses(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-4]
        int_id = self.path.split("/")[-2]

        int_response = self.apigatewayv2_backend.get_integration_responses(
            api_id=api_id, integration_id=int_id
        )
        return 200, {}, json.dumps({"items": [res.to_json() for res in int_response]})

    def update_integration_response(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-5]
        int_id = self.path.split("/")[-3]
        int_res_id = self.path.split("/")[-1]

        params = json.loads(self.body)
        content_handling_strategy = params.get("contentHandlingStrategy")
        integration_response_key = params.get("integrationResponseKey")
        response_parameters = params.get("responseParameters")
        response_templates = params.get("responseTemplates")
        template_selection_expression = params.get("templateSelectionExpression")
        integration_response = self.apigatewayv2_backend.update_integration_response(
            api_id=api_id,
            integration_id=int_id,
            integration_response_id=int_res_id,
            content_handling_strategy=content_handling_strategy,
            integration_response_key=integration_response_key,
            response_parameters=response_parameters,
            response_templates=response_templates,
            template_selection_expression=template_selection_expression,
        )
        return 200, {}, json.dumps(integration_response.to_json())

    def create_vpc_link(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)

        name = params.get("name")
        sg_ids = params.get("securityGroupIds")
        subnet_ids = params.get("subnetIds")
        tags = params.get("tags")
        vpc_link = self.apigatewayv2_backend.create_vpc_link(
            name, sg_ids, subnet_ids, tags
        )
        return 200, {}, json.dumps(vpc_link.to_json())

    def delete_vpc_link(self) -> TYPE_RESPONSE:
        vpc_link_id = self.path.split("/")[-1]
        self.apigatewayv2_backend.delete_vpc_link(vpc_link_id)
        return 200, {}, "{}"

    def get_vpc_link(self) -> TYPE_RESPONSE:
        vpc_link_id = self.path.split("/")[-1]
        vpc_link = self.apigatewayv2_backend.get_vpc_link(vpc_link_id)
        return 200, {}, json.dumps(vpc_link.to_json())

    def get_vpc_links(self) -> TYPE_RESPONSE:
        vpc_links = self.apigatewayv2_backend.get_vpc_links()
        return 200, {}, json.dumps({"items": [link.to_json() for link in vpc_links]})

    def update_vpc_link(self) -> TYPE_RESPONSE:
        vpc_link_id = self.path.split("/")[-1]
        params = json.loads(self.body)
        name = params.get("name")

        vpc_link = self.apigatewayv2_backend.update_vpc_link(vpc_link_id, name=name)
        return 200, {}, json.dumps(vpc_link.to_json())

    def create_domain_name(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)
        domain_name = params.get("domainName")
        domain_name_configurations = params.get("domainNameConfigurations", [{}])
        mutual_tls_authentication = params.get("mutualTlsAuthentication", {})
        tags = params.get("tags", {})
        domain_name = self.apigatewayv2_backend.create_domain_name(
            domain_name=domain_name,
            domain_name_configurations=domain_name_configurations,
            mutual_tls_authentication=mutual_tls_authentication,
            tags=tags,
        )
        return 201, {}, json.dumps(domain_name.to_json())

    def get_domain_name(self) -> TYPE_RESPONSE:
        domain_name_param = self.path.split("/")[-1]
        domain_name = self.apigatewayv2_backend.get_domain_name(
            domain_name=domain_name_param
        )
        return 200, {}, json.dumps(domain_name.to_json())

    def get_domain_names(self) -> TYPE_RESPONSE:
        domain_names = self.apigatewayv2_backend.get_domain_names()
        list_of_dict = [domain_name.to_json() for domain_name in domain_names]
        return 200, {}, json.dumps(dict(items=list_of_dict))

    def create_api_mapping(self) -> TYPE_RESPONSE:
        domain_name = self.path.split("/")[-2]
        params = json.loads(self.body)
        api_id = params.get("apiId")
        api_mapping_key = params.get("apiMappingKey", "")
        stage = params.get("stage")
        mapping = self.apigatewayv2_backend.create_api_mapping(
            api_id=api_id,
            api_mapping_key=api_mapping_key,
            domain_name=domain_name,
            stage=stage,
        )
        return 201, {}, json.dumps(mapping.to_json())

    def get_api_mapping(self) -> TYPE_RESPONSE:
        api_mapping_id = self.path.split("/")[-1]
        domain_name = self.path.split("/")[-3]
        mapping = self.apigatewayv2_backend.get_api_mapping(
            api_mapping_id=api_mapping_id,
            domain_name=domain_name,
        )
        return 200, {}, json.dumps(mapping.to_json())

    def get_api_mappings(self) -> TYPE_RESPONSE:
        domain_name = self.path.split("/")[-2]
        mappings = self.apigatewayv2_backend.get_api_mappings(domain_name=domain_name)
        list_of_dict = [mapping.to_json() for mapping in mappings]
        return 200, {}, json.dumps(dict(items=list_of_dict))

    def delete_domain_name(self) -> TYPE_RESPONSE:
        domain_name = self.path.split("/")[-1]
        self.apigatewayv2_backend.delete_domain_name(
            domain_name=domain_name,
        )
        return 204, {}, ""

    def delete_api_mapping(self) -> TYPE_RESPONSE:
        api_mapping_id = self.path.split("/")[-1]
        domain_name = self.path.split("/")[-3]
        self.apigatewayv2_backend.delete_api_mapping(
            api_mapping_id=api_mapping_id,
            domain_name=domain_name,
        )
        return 204, {}, ""

    def create_stage(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        config = json.loads(self.body)
        stage = self.apigatewayv2_backend.create_stage(api_id, config)
        return 200, {}, json.dumps(stage.to_json())

    def get_stage(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        stage_name = self.path.split("/")[-1]
        stage = self.apigatewayv2_backend.get_stage(api_id, stage_name)
        return 200, {}, json.dumps(stage.to_json())

    def delete_stage(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-3]
        stage_name = self.path.split("/")[-1]
        self.apigatewayv2_backend.delete_stage(api_id, stage_name)
        return 200, {}, "{}"

    def get_stages(self) -> TYPE_RESPONSE:
        api_id = self.path.split("/")[-2]
        stages = self.apigatewayv2_backend.get_stages(api_id)
        return 200, {}, json.dumps({"items": [st.to_json() for st in stages]})
