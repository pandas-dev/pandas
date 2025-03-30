"""Handles incoming appsync requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import AppSyncBackend, appsync_backends


class AppSyncResponse(BaseResponse):
    """Handler for AppSync requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="appsync")

    @property
    def appsync_backend(self) -> AppSyncBackend:
        """Return backend instance specific for this region."""
        return appsync_backends[self.current_account][self.region]

    def create_graphql_api(self) -> str:
        params = json.loads(self.body)
        name = params.get("name")
        log_config = params.get("logConfig")
        authentication_type = params.get("authenticationType")
        user_pool_config = params.get("userPoolConfig")
        open_id_connect_config = params.get("openIDConnectConfig")
        tags = params.get("tags")
        additional_authentication_providers = params.get(
            "additionalAuthenticationProviders"
        )
        xray_enabled = params.get("xrayEnabled", False)
        lambda_authorizer_config = params.get("lambdaAuthorizerConfig")
        visibility = params.get("visibility")
        graphql_api = self.appsync_backend.create_graphql_api(
            name=name,
            log_config=log_config,
            authentication_type=authentication_type,
            user_pool_config=user_pool_config,
            open_id_connect_config=open_id_connect_config,
            additional_authentication_providers=additional_authentication_providers,
            xray_enabled=xray_enabled,
            lambda_authorizer_config=lambda_authorizer_config,
            tags=tags,
            visibility=visibility,
        )
        response = graphql_api.to_json()
        response["tags"] = self.appsync_backend.list_tags_for_resource(graphql_api.arn)
        return json.dumps(dict(graphqlApi=response))

    def get_graphql_api(self) -> str:
        api_id = self.path.split("/")[-1]

        graphql_api = self.appsync_backend.get_graphql_api(api_id=api_id)
        response = graphql_api.to_json()
        response["tags"] = self.appsync_backend.list_tags_for_resource(graphql_api.arn)
        return json.dumps(dict(graphqlApi=response))

    def delete_graphql_api(self) -> str:
        api_id = self.path.split("/")[-1]
        self.appsync_backend.delete_graphql_api(api_id=api_id)
        return "{}"

    def update_graphql_api(self) -> str:
        api_id = self.path.split("/")[-1]

        params = json.loads(self.body)
        name = params.get("name")
        log_config = params.get("logConfig")
        authentication_type = params.get("authenticationType")
        user_pool_config = params.get("userPoolConfig")
        open_id_connect_config = params.get("openIDConnectConfig")
        additional_authentication_providers = params.get(
            "additionalAuthenticationProviders"
        )
        xray_enabled = params.get("xrayEnabled", False)
        lambda_authorizer_config = params.get("lambdaAuthorizerConfig")

        api = self.appsync_backend.update_graphql_api(
            api_id=api_id,
            name=name,
            log_config=log_config,
            authentication_type=authentication_type,
            user_pool_config=user_pool_config,
            open_id_connect_config=open_id_connect_config,
            additional_authentication_providers=additional_authentication_providers,
            xray_enabled=xray_enabled,
            lambda_authorizer_config=lambda_authorizer_config,
        )
        return json.dumps(dict(graphqlApi=api.to_json()))

    def list_graphql_apis(self) -> str:
        graphql_apis = self.appsync_backend.list_graphql_apis()
        return json.dumps(dict(graphqlApis=[api.to_json() for api in graphql_apis]))

    def create_api_key(self) -> str:
        params = json.loads(self.body)
        # /v1/apis/[api_id]/apikeys
        api_id = self.path.split("/")[-2]
        description = params.get("description")
        expires = params.get("expires")
        api_key = self.appsync_backend.create_api_key(
            api_id=api_id, description=description, expires=expires
        )
        return json.dumps(dict(apiKey=api_key.to_json()))

    def delete_api_key(self) -> str:
        api_id = self.path.split("/")[-3]
        api_key_id = self.path.split("/")[-1]
        self.appsync_backend.delete_api_key(api_id=api_id, api_key_id=api_key_id)
        return "{}"

    def list_api_keys(self) -> str:
        # /v1/apis/[api_id]/apikeys
        api_id = self.path.split("/")[-2]
        api_keys = self.appsync_backend.list_api_keys(api_id=api_id)
        return json.dumps(dict(apiKeys=[key.to_json() for key in api_keys]))

    def update_api_key(self) -> str:
        api_id = self.path.split("/")[-3]
        api_key_id = self.path.split("/")[-1]
        params = json.loads(self.body)
        description = params.get("description")
        expires = params.get("expires")
        api_key = self.appsync_backend.update_api_key(
            api_id=api_id,
            api_key_id=api_key_id,
            description=description,
            expires=expires,
        )
        return json.dumps(dict(apiKey=api_key.to_json()))

    def start_schema_creation(self) -> str:
        params = json.loads(self.body)
        api_id = self.path.split("/")[-2]
        definition = params.get("definition")
        status = self.appsync_backend.start_schema_creation(
            api_id=api_id, definition=definition
        )
        return json.dumps({"status": status})

    def get_schema_creation_status(self) -> str:
        api_id = self.path.split("/")[-2]
        status, details = self.appsync_backend.get_schema_creation_status(api_id=api_id)
        return json.dumps(dict(status=status, details=details))

    def tag_resource(self) -> str:
        resource_arn = self._extract_arn_from_path()
        params = json.loads(self.body)
        tags = params.get("tags")
        self.appsync_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return "{}"

    def untag_resource(self) -> str:
        resource_arn = self._extract_arn_from_path()
        tag_keys = self.querystring.get("tagKeys", [])
        self.appsync_backend.untag_resource(
            resource_arn=resource_arn, tag_keys=tag_keys
        )
        return "{}"

    def list_tags_for_resource(self) -> str:
        resource_arn = self._extract_arn_from_path()
        tags = self.appsync_backend.list_tags_for_resource(resource_arn=resource_arn)
        return json.dumps(dict(tags=tags))

    def _extract_arn_from_path(self) -> str:
        # /v1/tags/arn_that_may_contain_a_slash
        path = unquote(self.path)
        return "/".join(path.split("/")[3:])

    def get_type(self) -> str:
        api_id = unquote(self.path.split("/")[-3])
        type_name = self.path.split("/")[-1]
        type_format = self.querystring.get("format")[0]  # type: ignore[index]
        graphql_type = self.appsync_backend.get_type(
            api_id=api_id, type_name=type_name, type_format=type_format
        )
        return json.dumps(dict(type=graphql_type))

    def get_introspection_schema(self) -> str:
        api_id = self.path.split("/")[-2]
        format_ = self.querystring.get("format")[0]  # type: ignore[index]
        if self.querystring.get("includeDirectives"):
            include_directives = (
                self.querystring.get("includeDirectives")[0].lower() == "true"  # type: ignore[index]
            )
        else:
            include_directives = True
        graphql_schema = self.appsync_backend.get_graphql_schema(api_id=api_id)

        schema = graphql_schema.get_introspection_schema(
            format_=format_, include_directives=include_directives
        )
        return schema

    def get_api_cache(self) -> str:
        api_id = self.path.split("/")[-2]
        api_cache = self.appsync_backend.get_api_cache(
            api_id=api_id,
        )
        return json.dumps(dict(apiCache=api_cache.to_json()))

    def delete_api_cache(self) -> str:
        api_id = self.path.split("/")[-2]
        self.appsync_backend.delete_api_cache(
            api_id=api_id,
        )
        return "{}"

    def create_api_cache(self) -> str:
        params = json.loads(self.body)
        api_id = self.path.split("/")[-2]
        ttl = params.get("ttl")
        transit_encryption_enabled = params.get("transitEncryptionEnabled")
        at_rest_encryption_enabled = params.get("atRestEncryptionEnabled")
        api_caching_behavior = params.get("apiCachingBehavior")
        type = params.get("type")
        health_metrics_config = params.get("healthMetricsConfig")
        api_cache = self.appsync_backend.create_api_cache(
            api_id=api_id,
            ttl=ttl,
            transit_encryption_enabled=transit_encryption_enabled,
            at_rest_encryption_enabled=at_rest_encryption_enabled,
            api_caching_behavior=api_caching_behavior,
            type=type,
            health_metrics_config=health_metrics_config,
        )
        return json.dumps(dict(apiCache=api_cache.to_json()))

    def update_api_cache(self) -> str:
        api_id = self.path.split("/")[-3]
        params = json.loads(self.body)
        ttl = params.get("ttl")
        api_caching_behavior = params.get("apiCachingBehavior")
        type = params.get("type")
        health_metrics_config = params.get("healthMetricsConfig")
        api_cache = self.appsync_backend.update_api_cache(
            api_id=api_id,
            ttl=ttl,
            api_caching_behavior=api_caching_behavior,
            type=type,
            health_metrics_config=health_metrics_config,
        )
        return json.dumps(dict(apiCache=api_cache.to_json()))

    def flush_api_cache(self) -> str:
        api_id = self.path.split("/")[-2]
        self.appsync_backend.flush_api_cache(
            api_id=api_id,
        )
        return "{}"
