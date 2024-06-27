from typing import Any, Dict, Iterable, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import (
    AppNotFoundException,
    ConfigurationProfileNotFound,
    ConfigurationVersionNotFound,
)


class HostedConfigurationVersion(BaseModel):
    def __init__(
        self,
        app_id: str,
        config_id: str,
        version: int,
        description: str,
        content: str,
        content_type: str,
        version_label: str,
    ):
        self.app_id = app_id
        self.config_id = config_id
        self.version = version
        self.description = description
        self.content = content
        self.content_type = content_type
        self.version_label = version_label

    def get_headers(self) -> Dict[str, Any]:
        return {
            "application-id": self.app_id,
            "configuration-profile-id": self.config_id,
            "version-number": self.version,
            "description": self.description,
            "content-type": self.content_type,
            "VersionLabel": self.version_label,
        }


class ConfigurationProfile(BaseModel):
    def __init__(
        self,
        application_id: str,
        name: str,
        region: str,
        account_id: str,
        description: str,
        location_uri: str,
        retrieval_role_arn: str,
        validators: List[Dict[str, str]],
        _type: str,
    ):
        self.id = mock_random.get_random_hex(7)
        self.arn = f"arn:{get_partition(region)}:appconfig:{region}:{account_id}:application/{application_id}/configurationprofile/{self.id}"
        self.application_id = application_id
        self.name = name
        self.description = description
        self.location_uri = location_uri
        self.retrieval_role_arn = retrieval_role_arn
        self.validators = validators
        self._type = _type
        self.config_versions: Dict[int, HostedConfigurationVersion] = dict()

    def create_version(
        self,
        app_id: str,
        config_id: str,
        description: str,
        content: str,
        content_type: str,
        version_label: str,
    ) -> HostedConfigurationVersion:
        if self.config_versions:
            version = sorted(self.config_versions.keys())[-1] + 1
        else:
            version = 1
        self.config_versions[version] = HostedConfigurationVersion(
            app_id=app_id,
            config_id=config_id,
            version=version,
            description=description,
            content=content,
            content_type=content_type,
            version_label=version_label,
        )
        return self.config_versions[version]

    def get_version(self, version: int) -> HostedConfigurationVersion:
        if version not in self.config_versions:
            raise ConfigurationVersionNotFound
        return self.config_versions[version]

    def delete_version(self, version: int) -> None:
        self.config_versions.pop(version)

    def to_json(self) -> Dict[str, Any]:
        return {
            "Id": self.id,
            "Name": self.name,
            "ApplicationId": self.application_id,
            "Description": self.description,
            "LocationUri": self.location_uri,
            "RetrievalRoleArn": self.retrieval_role_arn,
            "Validators": self.validators,
            "Type": self._type,
        }


class Application(BaseModel):
    def __init__(
        self, name: str, description: Optional[str], region: str, account_id: str
    ):
        self.id = mock_random.get_random_hex(7)
        self.arn = f"arn:{get_partition(region)}:appconfig:{region}:{account_id}:application/{self.id}"
        self.name = name
        self.description = description

        self.config_profiles: Dict[str, ConfigurationProfile] = dict()

    def to_json(self) -> Dict[str, Any]:
        return {
            "Id": self.id,
            "Name": self.name,
            "Description": self.description,
        }


class AppConfigBackend(BaseBackend):
    """Implementation of AppConfig APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.applications: Dict[str, Application] = dict()
        self.tagger = TaggingService()

    def create_application(
        self, name: str, description: Optional[str], tags: Dict[str, str]
    ) -> Application:
        app = Application(
            name, description, region=self.region_name, account_id=self.account_id
        )
        self.applications[app.id] = app
        self.tag_resource(app.arn, tags)
        return app

    def delete_application(self, app_id: str) -> None:
        self.applications.pop(app_id, None)

    def get_application(self, app_id: str) -> Application:
        if app_id not in self.applications:
            raise AppNotFoundException
        return self.applications[app_id]

    def update_application(
        self, application_id: str, name: str, description: str
    ) -> Application:
        app = self.get_application(application_id)
        if name is not None:
            app.name = name
        if description is not None:
            app.description = description
        return app

    def create_configuration_profile(
        self,
        application_id: str,
        name: str,
        description: str,
        location_uri: str,
        retrieval_role_arn: str,
        validators: List[Dict[str, str]],
        _type: str,
        tags: Dict[str, str],
    ) -> ConfigurationProfile:
        config_profile = ConfigurationProfile(
            application_id=application_id,
            name=name,
            region=self.region_name,
            account_id=self.account_id,
            description=description,
            location_uri=location_uri,
            retrieval_role_arn=retrieval_role_arn,
            validators=validators,
            _type=_type,
        )
        self.tag_resource(config_profile.arn, tags)
        self.get_application(application_id).config_profiles[config_profile.id] = (
            config_profile
        )
        return config_profile

    def delete_configuration_profile(self, app_id: str, config_profile_id: str) -> None:
        self.get_application(app_id).config_profiles.pop(config_profile_id)

    def get_configuration_profile(
        self, app_id: str, config_profile_id: str
    ) -> ConfigurationProfile:
        app = self.get_application(app_id)
        if config_profile_id not in app.config_profiles:
            raise ConfigurationProfileNotFound
        return app.config_profiles[config_profile_id]

    def update_configuration_profile(
        self,
        application_id: str,
        config_profile_id: str,
        name: str,
        description: str,
        retrieval_role_arn: str,
        validators: List[Dict[str, str]],
    ) -> ConfigurationProfile:
        config_profile = self.get_configuration_profile(
            application_id, config_profile_id
        )
        if name is not None:
            config_profile.name = name
        if description is not None:
            config_profile.description = description
        if retrieval_role_arn is not None:
            config_profile.retrieval_role_arn = retrieval_role_arn
        if validators is not None:
            config_profile.validators = validators
        return config_profile

    def list_configuration_profiles(
        self, app_id: str
    ) -> Iterable[ConfigurationProfile]:
        app = self.get_application(app_id)
        return app.config_profiles.values()

    def create_hosted_configuration_version(
        self,
        app_id: str,
        config_profile_id: str,
        description: str,
        content: str,
        content_type: str,
        version_label: str,
    ) -> HostedConfigurationVersion:
        """
        The LatestVersionNumber-parameter is not yet implemented
        """
        profile = self.get_configuration_profile(app_id, config_profile_id)
        return profile.create_version(
            app_id=app_id,
            config_id=config_profile_id,
            description=description,
            content=content,
            content_type=content_type,
            version_label=version_label,
        )

    def get_hosted_configuration_version(
        self, app_id: str, config_profile_id: str, version: int
    ) -> HostedConfigurationVersion:
        profile = self.get_configuration_profile(
            app_id=app_id, config_profile_id=config_profile_id
        )
        return profile.get_version(version)

    def delete_hosted_configuration_version(
        self, app_id: str, config_profile_id: str, version: int
    ) -> None:
        profile = self.get_configuration_profile(
            app_id=app_id, config_profile_id=config_profile_id
        )
        profile.delete_version(version=version)

    def list_tags_for_resource(self, arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(arn)

    def tag_resource(self, arn: str, tags: Dict[str, str]) -> None:
        self.tagger.tag_resource(arn, TaggingService.convert_dict_to_tags_input(tags))

    def untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


appconfig_backends = BackendDict(AppConfigBackend, "appconfig")
