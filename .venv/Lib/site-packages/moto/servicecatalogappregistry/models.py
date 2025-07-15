"""AppRegistryBackend class with methods for supported APIs."""

import datetime
import re
from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.resourcegroups.models import FakeResourceGroup
from moto.servicecatalogappregistry.exceptions import (
    ResourceNotFoundException,
    ValidationException,
)
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition


class Application(BaseModel):
    def __init__(
        self,
        name: str,
        description: str,
        region: str,
        account_id: str,
    ):
        self.id = mock_random.get_random_string(
            length=27, include_digits=True, lower_case=True
        )
        self.arn = f"arn:{get_partition(region)}:servicecatalog:{region}:{account_id}:applications/{self.id}"
        self.name = name
        self.description = description
        self.creationTime = datetime.datetime.now()
        self.lastUpdateTime = self.creationTime
        self.tags: Dict[str, str] = dict()
        self.applicationTag: Dict[str, str] = {"awsApplication": self.arn}

        self.associated_resources: Dict[str, AssociatedResource] = dict()

    def to_json(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "arn": self.arn,
            "name": self.name,
            "description": self.description,
            "creationTime": str(self.creationTime),
            "lastUpdateTime": str(self.lastUpdateTime),
            "tags": self.tags,
            "applicationTag": self.applicationTag,
        }


class AssociatedResource(BaseBackend):
    def __init__(
        self,
        resource_type: str,
        resource: str,
        options: List[str],
        application: Application,
        account_id: str,
        region_name: str,
    ):
        if resource_type == "CFN_STACK":
            from moto.cloudformation.exceptions import ValidationError
            from moto.cloudformation.models import cloudformation_backends

            self.resource = resource
            match = re.search(
                r"^arn:aws:cloudformation:(us(-gov)?|ap|ca|cn|eu|sa)-(central|(north|south)?(east|west)?)-\d:\d{12}:stack/\w[a-zA-Z0-9\-]{0,127}/[a-f0-9]{8}(-[a-f0-9]{4}){3}-[a-f0-9]{12}$",
                resource,
            )
            if match is not None:
                self.name = resource.split("/")[1]
            else:
                self.name = resource
            cf_backend = cloudformation_backends[account_id][region_name]
            try:
                cf_backend.get_stack(self.name)
            except ValidationError:
                raise ResourceNotFoundException(
                    f"No CloudFormation stack called '{self.name}' found"
                )
        elif resource_type == "RESOURCE_TAG_VALUE":
            tags = {
                "EnableAWSServiceCatalogAppRegistry": "true",
                "aws:servicecatalog:applicationName": application.name,
                "aws:servicecatalog:applicationId": application.id,
                "aws:servicecatalog:applicationArn": application.arn,
            }
            new_resource_group = FakeResourceGroup(
                account_id,
                region_name,
                f"AWS_AppRegistry_AppTag_{account_id}-{application.name}",
                {"Type": "TAG_FILTERS_1_0", "Query": resource},
                None,
                tags,
            )
            self.resource = new_resource_group.arn
            self.name = new_resource_group._name
            self.query = resource
        else:
            raise ValidationException
        self.resource_type = resource_type
        self.options = options

    def to_json(self) -> Dict[str, Any]:
        return_dict = {
            "name": self.name,
            "arn": self.resource,
            "resourceType": self.resource_type,
            "options": self.options,
        }
        if self.resource_type == "RESOURCE_TAG_VALUE":
            return_dict["resourceDetails"] = {"tagValue": self.query}  # type: ignore
        return return_dict


class AppRegistryBackend(BaseBackend):
    """Implementation of AppRegistry APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.applications: Dict[str, Application] = dict()
        self.tagger = TaggingService()

    def create_application(
        self, name: str, description: str, tags: Dict[str, str], client_token: str
    ) -> Application:
        app = Application(
            name,
            description,
            region=self.region_name,
            account_id=self.account_id,
        )
        self.applications[app.arn] = app
        self._tag_resource(app.arn, tags)
        return app

    def list_applications(self) -> List[Application]:
        return list(self.applications.values())

    def associate_resource(
        self, application: str, resource_type: str, resource: str, options: List[str]
    ) -> Dict[str, Any]:
        app = self.applications[application]
        new_resource = AssociatedResource(
            resource_type, resource, options, app, self.account_id, self.region_name
        )
        app.associated_resources[new_resource.resource] = new_resource
        return {"applicationArn": app.arn, "resourceArn": resource, "options": options}

    def _list_tags_for_resource(self, arn: str) -> Dict[str, str]:
        return self.tagger.get_tag_dict_for_resource(arn)

    def _tag_resource(self, arn: str, tags: Dict[str, str]) -> None:
        self.tagger.tag_resource(arn, TaggingService.convert_dict_to_tags_input(tags))

    def _untag_resource(self, arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


servicecatalogappregistry_backends = BackendDict(
    AppRegistryBackend, "servicecatalog-appregistry"
)
