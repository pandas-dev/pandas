"""Handles incoming clouddirectory requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import CloudDirectoryBackend, clouddirectory_backends


class CloudDirectoryResponse(BaseResponse):
    """Handler for CloudDirectory requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="clouddirectory")

    @property
    def clouddirectory_backend(self) -> "CloudDirectoryBackend":
        """Return backend instance specific for this region."""
        return clouddirectory_backends[self.current_account][self.region]

    def apply_schema(self) -> str:
        directory_arn = self.headers.get("x-amz-data-partition")
        published_schema_arn = self._get_param("PublishedSchemaArn")
        self.clouddirectory_backend.apply_schema(
            directory_arn=directory_arn,
            published_schema_arn=published_schema_arn,
        )
        return json.dumps(
            {
                "AppliedSchemaArn": published_schema_arn,
                "DirectoryArn": directory_arn,
            }
        )

    def publish_schema(self) -> str:
        development_schema_arn = self.headers.get("x-amz-data-partition")
        version = self._get_param("Version")
        minor_version = self._get_param("MinorVersion")
        name = self._get_param("Name")
        schema = self.clouddirectory_backend.publish_schema(
            name=name,
            version=version,
            minor_version=minor_version,
            development_schema_arn=development_schema_arn,
        )
        return json.dumps({"PublishedSchemaArn": schema})

    def create_directory(self) -> str:
        name = self._get_param("Name")
        schema_arn = self.headers.get("x-amz-data-partition")
        directory = self.clouddirectory_backend.create_directory(
            name=name,
            schema_arn=schema_arn,
        )

        return json.dumps(
            dict(
                DirectoryArn=directory.directory_arn,
                Name=name,
                ObjectIdentifier=directory.object_identifier,
                AppliedSchemaArn=directory.schema_arn,
            )
        )

    def create_schema(self) -> str:
        name = self._get_param("Name")
        schema = self.clouddirectory_backend.create_schema(
            name=name,
        )
        return json.dumps(dict(SchemaArn=schema))

    def list_directories(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        state = self._get_param("state")
        directories, next_token = self.clouddirectory_backend.list_directories(
            state=state, next_token=next_token, max_results=max_results
        )
        directory_list = [directory.to_dict() for directory in directories]
        return json.dumps(dict(Directories=directory_list, NextToken=next_token))

    def tag_resource(self) -> str:
        resource_arn = self._get_param("ResourceArn")
        tags = self._get_param("Tags")
        self.clouddirectory_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps(dict())

    def untag_resource(self) -> str:
        resource_arn = self._get_param("ResourceArn")
        tag_keys = self._get_param("TagKeys")
        self.clouddirectory_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return json.dumps(dict())

    def delete_directory(self) -> str:
        # Retrieve arn from headers
        # https://docs.aws.amazon.com/clouddirectory/latest/APIReference/API_DeleteDirectory.html
        arn = self.headers.get("x-amz-data-partition")
        directory_arn = self.clouddirectory_backend.delete_directory(
            directory_arn=arn,
        )
        return json.dumps(dict(DirectoryArn=directory_arn))

    def delete_schema(self) -> str:
        # Retrieve arn from headers
        # https://docs.aws.amazon.com/clouddirectory/latest/APIReference/API_DeleteSchema.html
        arn = self.headers.get("x-amz-data-partition")
        self.clouddirectory_backend.delete_schema(
            schema_arn=arn,
        )
        return json.dumps(dict(SchemaArn=arn))

    def list_development_schema_arns(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        schemas, next_token = self.clouddirectory_backend.list_development_schema_arns(
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(dict(SchemaArns=schemas, NextToken=next_token))

    def list_published_schema_arns(self) -> str:
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        schemas, next_token = self.clouddirectory_backend.list_published_schema_arns(
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(dict(SchemaArns=schemas, NextToken=next_token))

    def get_directory(self) -> str:
        # Retrieve arn from headers
        # https://docs.aws.amazon.com/clouddirectory/latest/APIReference/API_GetDirectory.html
        arn = self.headers.get("x-amz-data-partition")
        directory = self.clouddirectory_backend.get_directory(
            directory_arn=arn,
        )
        return json.dumps(dict(Directory=directory.to_dict()))

    def list_tags_for_resource(self) -> str:
        resource_arn = self._get_param("ResourceArn")
        next_token = self._get_param("NextToken")
        max_results = self._get_param("MaxResults")
        tags = self.clouddirectory_backend.list_tags_for_resource(
            resource_arn=resource_arn,
            next_token=next_token,
            max_results=max_results,
        )
        return json.dumps(dict(Tags=tags, NextToken=next_token))
