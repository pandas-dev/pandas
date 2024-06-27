import json
from typing import Any, Dict, Tuple
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import AppConfigBackend, appconfig_backends


class AppConfigResponse(BaseResponse):
    def tags(self, request: Any, full_url: str, headers: Any) -> str:  # type: ignore[return]
        self.setup_class(request, full_url, headers)
        if request.method == "GET":
            return self.list_tags_for_resource()
        if request.method == "POST":
            return self.tag_resource()
        if request.method == "DELETE":
            return self.untag_resource()

    def __init__(self) -> None:
        super().__init__(service_name="appconfig")

    @property
    def appconfig_backend(self) -> AppConfigBackend:
        return appconfig_backends[self.current_account][self.region]

    def create_application(self) -> str:
        name = self._get_param("Name")
        description = self._get_param("Description")
        tags = self._get_param("Tags")
        app = self.appconfig_backend.create_application(
            name=name,
            description=description,
            tags=tags,
        )
        return json.dumps(app.to_json())

    def delete_application(self) -> str:
        app_id = self._get_param("ApplicationId")
        self.appconfig_backend.delete_application(app_id)
        return "{}"

    def get_application(self) -> str:
        app_id = self._get_param("ApplicationId")
        app = self.appconfig_backend.get_application(app_id)
        return json.dumps(app.to_json())

    def update_application(self) -> str:
        app_id = self._get_param("ApplicationId")
        name = self._get_param("Name")
        description = self._get_param("Description")
        app = self.appconfig_backend.update_application(
            application_id=app_id,
            name=name,
            description=description,
        )
        return json.dumps(app.to_json())

    def create_configuration_profile(self) -> str:
        app_id = self._get_param("ApplicationId")
        name = self._get_param("Name")
        description = self._get_param("Description")
        location_uri = self._get_param("LocationUri")
        retrieval_role_arn = self._get_param("RetrievalRoleArn")
        validators = self._get_param("Validators")
        _type = self._get_param("Type")
        tags = self._get_param("Tags")
        config_profile = self.appconfig_backend.create_configuration_profile(
            application_id=app_id,
            name=name,
            description=description,
            location_uri=location_uri,
            retrieval_role_arn=retrieval_role_arn,
            validators=validators,
            _type=_type,
            tags=tags,
        )
        return json.dumps(config_profile.to_json())

    def delete_configuration_profile(self) -> str:
        app_id = self._get_param("ApplicationId")
        config_profile_id = self._get_param("ConfigurationProfileId")
        self.appconfig_backend.delete_configuration_profile(app_id, config_profile_id)
        return "{}"

    def get_configuration_profile(self) -> str:
        app_id = self._get_param("ApplicationId")
        config_profile_id = self._get_param("ConfigurationProfileId")
        config_profile = self.appconfig_backend.get_configuration_profile(
            app_id, config_profile_id
        )
        return json.dumps(config_profile.to_json())

    def update_configuration_profile(self) -> str:
        app_id = self._get_param("ApplicationId")
        config_profile_id = self._get_param("ConfigurationProfileId")
        name = self._get_param("Name")
        description = self._get_param("Description")
        retrieval_role_arn = self._get_param("RetrievalRoleArn")
        validators = self._get_param("Validators")
        config_profile = self.appconfig_backend.update_configuration_profile(
            application_id=app_id,
            config_profile_id=config_profile_id,
            name=name,
            description=description,
            retrieval_role_arn=retrieval_role_arn,
            validators=validators,
        )
        return json.dumps(config_profile.to_json())

    def list_configuration_profiles(self) -> str:
        app_id = self._get_param("ApplicationId")
        profiles = self.appconfig_backend.list_configuration_profiles(app_id)
        return json.dumps({"Items": [p.to_json() for p in profiles]})

    def list_tags_for_resource(self) -> str:
        arn = unquote(self.path.split("/tags/")[-1])
        tags = self.appconfig_backend.list_tags_for_resource(arn)
        return json.dumps({"Tags": tags})

    def tag_resource(self) -> str:
        arn = unquote(self.path.split("/tags/")[-1])
        tags = self._get_param("Tags")
        self.appconfig_backend.tag_resource(arn, tags)
        return "{}"

    def untag_resource(self) -> str:
        arn = unquote(self.path.split("/tags/")[-1])
        tag_keys = self.querystring.get("tagKeys")
        self.appconfig_backend.untag_resource(arn, tag_keys)  # type: ignore[arg-type]
        return "{}"

    def create_hosted_configuration_version(self) -> Tuple[str, Dict[str, Any]]:
        app_id = self._get_param("ApplicationId")
        config_profile_id = self._get_param("ConfigurationProfileId")
        description = self.headers.get("Description")
        content = self.body
        content_type = self.headers.get("Content-Type")
        version_label = self.headers.get("VersionLabel")
        version = self.appconfig_backend.create_hosted_configuration_version(
            app_id=app_id,
            config_profile_id=config_profile_id,
            description=description,
            content=content,
            content_type=content_type,
            version_label=version_label,
        )
        return version.content, version.get_headers()

    def get_hosted_configuration_version(self) -> Tuple[str, Dict[str, Any]]:
        app_id = self._get_param("ApplicationId")
        config_profile_id = self._get_param("ConfigurationProfileId")
        version_number = self._get_int_param("VersionNumber")
        version = self.appconfig_backend.get_hosted_configuration_version(
            app_id=app_id,
            config_profile_id=config_profile_id,
            version=version_number,
        )
        return version.content, version.get_headers()

    def delete_hosted_configuration_version(self) -> str:
        app_id = self._get_param("ApplicationId")
        config_profile_id = self._get_param("ConfigurationProfileId")
        version_number = self._get_int_param("VersionNumber")
        self.appconfig_backend.delete_hosted_configuration_version(
            app_id=app_id,
            config_profile_id=config_profile_id,
            version=version_number,
        )
        return "{}"
