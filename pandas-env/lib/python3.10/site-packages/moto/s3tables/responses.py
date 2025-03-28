"""Handles incoming s3tables requests, invokes methods, returns responses."""

import json
from typing import Any, Dict
from urllib.parse import unquote

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import S3TablesBackend, s3tables_backends


class S3TablesResponse(BaseResponse):
    """Handler for S3Tables requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="s3tables")
        self.default_response_headers = {"Content-Type": "application/json"}

    @property
    def s3tables_backend(self) -> S3TablesBackend:
        """Return backend instance specific for this region."""
        return s3tables_backends[self.current_account][self.region]

    def create_table_bucket(self) -> TYPE_RESPONSE:
        name = json.loads(self.body)["name"]
        bucket = self.s3tables_backend.create_table_bucket(
            name=name,
        )
        return 200, self.default_response_headers, json.dumps(dict(arn=bucket.arn))

    def list_table_buckets(self) -> TYPE_RESPONSE:
        params = self._get_params()
        prefix = params.get("prefix")
        continuation_token = params.get("continuationToken")
        max_buckets = params.get("maxBuckets")
        table_buckets, continuation_token = self.s3tables_backend.list_table_buckets(
            prefix=prefix,
            continuation_token=continuation_token,
            max_buckets=int(max_buckets) if max_buckets else None,
        )

        body: Dict[str, Any] = {
            "tableBuckets": [
                {
                    "arn": b.arn,
                    "name": b.name,
                    "ownerAccountId": b.account_id,
                    "createdAt": b.creation_date.isoformat(),
                }
                for b in table_buckets
            ]
        }
        if continuation_token:
            body.update(continuationToken=continuation_token)

        return 200, self.default_response_headers, json.dumps(body)

    def get_table_bucket(self) -> TYPE_RESPONSE:
        _, table_bucket_arn = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        bucket = self.s3tables_backend.get_table_bucket(
            table_bucket_arn=table_bucket_arn,
        )

        return (
            200,
            self.default_response_headers,
            json.dumps(
                dict(
                    arn=bucket.arn,
                    name=bucket.name,
                    ownerAccountId=bucket.account_id,
                    createdAt=bucket.creation_date.isoformat(),
                )
            ),
        )

    def delete_table_bucket(self) -> TYPE_RESPONSE:
        _, table_bucket_arn = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        self.s3tables_backend.delete_table_bucket(
            table_bucket_arn=table_bucket_arn,
        )

        return 204, {}, ""

    def create_namespace(self) -> TYPE_RESPONSE:
        _, table_bucket_arn = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        name = json.loads(self.body)["namespace"][0]
        namespace = self.s3tables_backend.create_namespace(
            table_bucket_arn=table_bucket_arn,
            namespace=name,
        )
        return (
            200,
            self.default_response_headers,
            json.dumps(
                dict(tableBucketArn=table_bucket_arn, namespace=[namespace.name])
            ),
        )

    def list_namespaces(self) -> TYPE_RESPONSE:
        _, table_bucket_arn = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)

        params = self._get_params()
        continuation_token = params.get("continuationToken")
        max_namespaces = params.get("maxNamespaces")
        prefix = params.get("prefix")

        namespaces, continuation_token = self.s3tables_backend.list_namespaces(
            table_bucket_arn=table_bucket_arn,
            prefix=prefix,
            continuation_token=continuation_token,
            max_namespaces=int(max_namespaces) if max_namespaces else None,
        )

        body: Dict[str, Any] = {
            "namespaces": [
                {
                    "namespace": [ns.name],
                    "createdAt": ns.creation_date.isoformat(),
                    "createdBy": ns.created_by,
                    "ownerAccountId": ns.account_id,
                }
                for ns in namespaces
            ]
        }
        if continuation_token:
            body.update(continuationToken=continuation_token)

        return 200, self.default_response_headers, json.dumps(body)

    def get_namespace(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, name = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        namespace = self.s3tables_backend.get_namespace(
            table_bucket_arn=table_bucket_arn,
            namespace=name,
        )
        return (
            200,
            self.default_response_headers,
            json.dumps(
                dict(
                    namespace=[namespace.name],
                    createdAt=namespace.creation_date.isoformat(),
                    createdBy=namespace.created_by,
                    ownerAccountId=namespace.account_id,
                )
            ),
        )

    def delete_namespace(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        self.s3tables_backend.delete_namespace(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
        )
        return 204, self.default_response_headers, ""

    def create_table(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        body = json.loads(self.body)
        name = body["name"]
        format = body["format"]
        table = self.s3tables_backend.create_table(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            name=name,
            format=format,
        )
        return (
            200,
            self.default_response_headers,
            json.dumps(dict(tableARN=table.arn, versionToken=table.version_token)),
        )

    def get_table(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace, name = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        table = self.s3tables_backend.get_table(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            name=name,
        )
        return (
            200,
            self.default_response_headers,
            json.dumps(
                dict(
                    name=table.name,
                    type=table.type,
                    tableARN=table.arn,
                    namespace=[namespace],
                    versionToken=table.version_token,
                    metadataLocation=table.metadata_location,
                    warehouseLocation=table.warehouse_location,
                    createdAt=table.creation_date.isoformat(),
                    createdBy=table.account_id,
                    managedByService=table.managed_by_service,
                    modifiedAt=table.last_modified.isoformat(),
                    modifiedBy=table.modified_by,
                    ownerAccountId=table.account_id,
                    format=table.format,
                )
            ),
        )

    def list_tables(self) -> TYPE_RESPONSE:
        _, table_bucket_arn = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        params = self._get_params()
        namespace = params.get("namespace")
        prefix = params.get("prefix")
        continuation_token = params.get("continuationToken")
        max_tables = params.get("maxTables")
        tables, continuation_token = self.s3tables_backend.list_tables(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            prefix=prefix,
            continuation_token=continuation_token,
            max_tables=int(max_tables) if max_tables else None,
        )
        body: Dict[str, Any] = {
            "tables": [
                {
                    "namespace": [table.namespace],
                    "name": table.name,
                    "createdAt": table.creation_date.isoformat(),
                    "modifiedAt": table.last_modified.isoformat(),
                }
                for table in tables
            ]
        }

        if continuation_token:
            body.update(continuationToken=continuation_token)

        return 200, self.default_response_headers, json.dumps(body)

    def delete_table(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace, name = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        params = self._get_params()
        version_token = params.get("versionToken")
        self.s3tables_backend.delete_table(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            name=name,
            version_token=version_token,
        )
        return 204, {}, ""

    def get_table_metadata_location(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace, name, _ = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        table = self.s3tables_backend.get_table(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            name=name,
        )
        return (
            200,
            self.default_response_headers,
            json.dumps(
                dict(
                    versionToken=table.version_token,
                    metadataLocation=table.metadata_location,
                    warehouseLocation=table.warehouse_location,
                )
            ),
        )

    def update_table_metadata_location(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace, name, _ = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        body = json.loads(self.body)
        metadata_location = body["metadataLocation"]
        version_token = body["versionToken"]
        table = self.s3tables_backend.update_table_metadata_location(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            name=name,
            version_token=version_token,
            metadata_location=metadata_location,
        )
        return (
            200,
            self.default_response_headers,
            json.dumps(
                dict(
                    name=table.name,
                    tableArn=table.arn,
                    namespace=namespace,
                    versionToken=table.version_token,
                    metadataLocation=table.metadata_location,
                )
            ),
        )

    def rename_table(self) -> TYPE_RESPONSE:
        _, table_bucket_arn, namespace, name, _ = self.raw_path.lstrip("/").split("/")
        table_bucket_arn = unquote(table_bucket_arn)
        body = json.loads(self.body)
        version_token = body.get("versionToken")
        new_namespace_name = body.get("newNamespaceName")
        new_name = body.get("newName")
        self.s3tables_backend.rename_table(
            table_bucket_arn=table_bucket_arn,
            namespace=namespace,
            name=name,
            new_namespace_name=new_namespace_name,
            new_name=new_name,
            version_token=version_token,
        )
        return 200, {}, ""
