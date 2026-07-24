"""CloudDirectoryBackend class with methods for supported APIs."""

import datetime
from collections.abc import Iterator

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.resource_tagging import TaggableResourcesMixin, TaggedResource
from moto.utilities.paginator import paginate
from moto.utilities.tagging_service import TaggingService

from .exceptions import (
    InvalidArnException,
    ResourceNotFoundException,
    SchemaAlreadyPublishedException,
)

PAGINATION_MODEL = {
    "list_directories": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "directory_arn",
    },
    "list_development_schema_arns": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "schema_arn",
    },
    "list_published_schema_arns": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "unique_attribute": "schema_arn",
    },
}


class Directory(BaseModel):
    def __init__(
        self, account_id: str, region: str, name: str, schema_arn: str
    ) -> None:
        self.name = name
        self.schema_arn = schema_arn
        self.directory_arn = (
            f"arn:aws:clouddirectory:{region}:{account_id}:directory/{name}"
        )
        self.state = "ENABLED"
        self.creation_date_time = datetime.datetime.now()
        self.object_identifier = f"directory-{name}"

    def to_dict(self) -> dict[str, str]:
        return {
            "Name": self.name,
            "SchemaArn": self.schema_arn,
            "DirectoryArn": self.directory_arn,
            "State": self.state,
            "CreationDateTime": str(self.creation_date_time),
            "ObjectIdentifier": self.object_identifier,
        }


class CloudDirectoryBackend(BaseBackend, TaggableResourcesMixin):
    """Implementation of CloudDirectory APIs."""

    SERVICE_NAMESPACE = "clouddirectory"

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.directories: dict[str, Directory] = {}
        self.schemas_states: dict[str, list[str]] = {
            "development": [],
            "published": [],
            "applied": [],
        }
        self.tagger = TaggingService()

    def apply_schema(self, directory_arn: str, published_schema_arn: str) -> None:
        directory = self.directories.get(directory_arn)
        if not directory:
            raise ResourceNotFoundException(directory_arn)
        if published_schema_arn not in self.schemas_states["published"]:
            raise ResourceNotFoundException(published_schema_arn)
        directory.schema_arn = published_schema_arn
        return

    def publish_schema(
        self, name: str, version: str, development_schema_arn: str, minor_version: str
    ) -> str:
        schema_arn = f"arn:aws:clouddirectory:{self.region_name}:{self.account_id}:schema/published/{name}/{version}/{minor_version}"
        if development_schema_arn in self.schemas_states["published"]:
            raise SchemaAlreadyPublishedException(development_schema_arn)
        if development_schema_arn in self.schemas_states["development"]:
            self.schemas_states["development"].remove(development_schema_arn)
            self.schemas_states["published"].append(schema_arn)
        else:
            raise ResourceNotFoundException(development_schema_arn)
        return schema_arn

    def create_directory(self, name: str, schema_arn: str) -> Directory:
        directory = Directory(self.account_id, self.region_name, name, schema_arn)
        self.directories[directory.directory_arn] = directory
        return directory

    def create_schema(self, name: str) -> str:
        self.schema_arn = f"arn:aws:clouddirectory:{self.region_name}:{self.account_id}:schema/development/{name}"
        self.schemas_states["development"].append(self.schema_arn)
        return self.schema_arn

    def delete_schema(self, schema_arn: str) -> None:
        if schema_arn in self.schemas_states["development"]:
            self.schemas_states["development"].remove(schema_arn)
        elif schema_arn in self.schemas_states["published"]:
            self.schemas_states["published"].remove(schema_arn)
        else:
            raise ResourceNotFoundException(schema_arn)
        return

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_development_schema_arns(self) -> list[str]:
        return self.schemas_states["development"]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_published_schema_arns(self) -> list[str]:
        return self.schemas_states["published"]

    @paginate(pagination_model=PAGINATION_MODEL)
    def list_directories(self, state: str) -> list[Directory]:
        directories = list(self.directories.values())
        if state:
            directories = [
                directory for directory in directories if directory.state == state
            ]
        return directories

    def delete_directory(self, directory_arn: str) -> str:
        directory = self.directories.pop(directory_arn)
        return directory.directory_arn

    def get_directory(self, directory_arn: str) -> Directory:
        directory = self.directories.get(directory_arn)
        if not directory:
            raise InvalidArnException(directory_arn)
        return directory

    def list_tags_for_resource(
        self, resource_arn: str, next_token: str, max_results: int
    ) -> list[dict[str, str]]:
        tags = self.tagger.list_tags_for_resource(resource_arn)["Tags"]
        return tags

    # Resource Groups Tagging API (TaggableResourcesMixin method overrides)
    def iter_tagged_resources(self) -> Iterator[TaggedResource]:
        for directory in self.directories.values():
            yield TaggedResource(
                arn=directory.directory_arn,
                tags=self.tagger.get_tag_dict_for_resource(directory.directory_arn),
                resource_type="clouddirectory:directory",
            )

    def tag_resource(self, arn: str, tags: dict[str, str]) -> None:
        self.tagger.tag_resource(arn, self.tagger.convert_dict_to_tags_input(tags))

    def untag_resource(self, arn: str, tag_keys: list[str]) -> None:
        self.tagger.untag_resource_using_names(arn, tag_keys)


clouddirectory_backends = BackendDict(CloudDirectoryBackend, "clouddirectory")
