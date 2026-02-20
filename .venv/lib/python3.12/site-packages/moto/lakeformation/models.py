from __future__ import annotations

import json
from collections import defaultdict
from enum import Enum
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.tagging_service import TaggingService

from .exceptions import AlreadyExists, EntityNotFound, InvalidInput


class RessourceType(Enum):
    catalog = "CATALOG"
    database = "DATABASE"
    table = "TABLE"
    data_location = "DATA_LOCATION"


class Resource(BaseModel):
    def __init__(self, arn: str, role_arn: str):
        self.arn = arn
        self.role_arn = role_arn

    def to_dict(self) -> dict[str, Any]:
        return {
            "ResourceArn": self.arn,
            "RoleArn": self.role_arn,
        }


class Permission:
    def __init__(
        self,
        principal: dict[str, str],
        resource: dict[str, Any],
        permissions: list[str],
        permissions_with_grant_options: list[str],
    ):
        self.principal = principal
        self.resource = resource
        self.permissions = permissions
        self.permissions_with_grant_options = permissions_with_grant_options

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, Permission):
            return (
                (self.principal == other.principal)
                and (self.resource == other.resource)
                and (self.permissions == other.permissions)
                and (
                    self.permissions_with_grant_options
                    == other.permissions_with_grant_options
                )
            )
        return False

    def __hash__(self) -> int:
        return hash(
            (
                json.dumps(self.principal),
                json.dumps(self.resource),
                json.dumps(self.permissions),
                json.dumps(self.permissions_with_grant_options),
            )
        )

    def equal_principal_and_resouce(self, other: Permission) -> bool:
        return (self.principal == other.principal) and (self.resource == other.resource)

    def merge(self, other: Permission) -> None:
        self.permissions = list(set(self.permissions).union(other.permissions))
        self.permissions_with_grant_options = list(
            set(self.permissions_with_grant_options).union(
                other.permissions_with_grant_options
            )
        )

    def diff(self, other: Permission) -> None:
        if self.permissions is not None:
            self.permissions = list(set(self.permissions).difference(other.permissions))
        if self.permissions_with_grant_options is not None:
            self.permissions_with_grant_options = list(
                set(self.permissions_with_grant_options).difference(
                    other.permissions_with_grant_options
                )
            )

    def is_empty(self) -> bool:
        return (
            len(self.permissions) == 0 and len(self.permissions_with_grant_options) == 0
        )

    def to_external_form(self) -> dict[str, Any]:
        return {
            "Permissions": self.permissions,
            "PermissionsWithGrantOption": self.permissions_with_grant_options,
            "Resource": self.resource,
            "Principal": self.principal,
        }


class PermissionCatalog:
    def __init__(self) -> None:
        self.permissions: set[Permission] = set()

    def add_permission(self, permission: Permission) -> None:
        for existing_permission in self.permissions:
            if permission.equal_principal_and_resouce(existing_permission):
                # Permission with same principal and resouce, only once of these can exist
                existing_permission.merge(permission)
                return
        # found no match
        self.permissions.add(permission)

    def remove_permission(self, permission: Permission) -> None:
        for existing_permission in self.permissions:
            if permission.equal_principal_and_resouce(existing_permission):
                # Permission with same principal and resouce, only once of these can exist
                # remove and readd to recalculate the hash value after the diff
                self.permissions.remove(existing_permission)
                existing_permission.diff(permission)
                self.permissions.add(existing_permission)
                if existing_permission.is_empty():
                    self.permissions.remove(existing_permission)
                return


class ListPermissionsResourceDatabase:
    def __init__(self, catalog_id: Optional[str], name: str):
        self.name = name
        self.catalog_id = catalog_id


class ListPermissionsResourceTable:
    def __init__(
        self,
        catalog_id: Optional[str],
        database_name: str,
        name: Optional[str],
        table_wildcard: Optional[
            dict[str, str]
        ],  # Placeholder type, table_wildcard is an empty dict in docs
    ):
        if name is None and table_wildcard is None:
            raise InvalidInput("Table name and table wildcard cannot both be empty.")
        if name is not None and table_wildcard is not None:
            raise InvalidInput("Table name and table wildcard cannot both be present.")
        self.database_name = database_name
        self.name = name
        self.catalog_id = catalog_id
        self.table_wildcard = table_wildcard


class ExcludedColumnNames:
    def __init__(self, excluded_column_names: list[str]):
        self.excluded_column_names = excluded_column_names


class ListPermissionsResourceTableWithColumns:
    def __init__(
        self,
        catalog_id: Optional[str],
        database_name: str,
        name: str,
        column_names: list[str],
        column_wildcard: ExcludedColumnNames,
    ):
        self.database_name = database_name
        self.name = name
        self.catalog_id = catalog_id
        self.column_names = column_names
        self.column_wildcard = column_wildcard


class ListPermissionsResourceDataLocation:
    def __init__(self, catalog_id: Optional[str], resource_arn: str):
        self.catalog_id = catalog_id
        self.resource_arn = resource_arn


class ListPermissionsResourceDataCellsFilter:
    def __init__(
        self, table_catalog_id: str, database_name: str, table_name: str, name: str
    ):
        self.table_catalog_id = table_catalog_id
        self.database_name = database_name
        self.table_name = table_name
        self.name = name


class ListPermissionsResourceLFTag:
    def __init__(self, catalog_id: str, tag_key: str, tag_values: list[str]):
        self.catalog_id = catalog_id
        self.tag_key = tag_key
        self.tag_values = tag_values


class LFTag:
    def __init__(self, tag_key: str, tag_values: list[str]):
        self.tag_key = tag_key
        self.tag_values = tag_values


class ListPermissionsResourceLFTagPolicy:
    def __init__(self, catalog_id: str, resource_type: str, expression: list[LFTag]):
        self.catalog_id = catalog_id
        self.resource_type = resource_type
        self.expression = expression


class ListPermissionsResource:
    def __init__(
        self,
        catalog: Optional[
            dict[str, str]
        ],  # Placeholder type, catalog is an empty dict in docs
        database: Optional[ListPermissionsResourceDatabase],
        table: Optional[ListPermissionsResourceTable],
        table_with_columns: Optional[ListPermissionsResourceTableWithColumns],
        data_location: Optional[ListPermissionsResourceDataLocation],
        data_cells_filter: Optional[ListPermissionsResourceDataCellsFilter],
        lf_tag: Optional[ListPermissionsResourceLFTag],
        lf_tag_policy: Optional[ListPermissionsResourceLFTagPolicy],
    ):
        if (
            catalog is None
            and database is None
            and table is None
            and data_location is None
        ):
            # Error message is the exact string returned by the AWS-CLI eventhough it is valid
            # to not populate the respective fields as long as data_location is given.
            raise InvalidInput(
                "Resource must have either the catalog, table or database field populated."
            )
        self.catalog = catalog
        self.database = database
        self.table = table
        self.table_with_columns = table_with_columns
        self.data_location = data_location
        self.data_cells_filter = data_cells_filter
        self.lf_tag = lf_tag
        self.lf_tag_policy = lf_tag_policy


def default_settings() -> dict[str, Any]:
    return {
        "DataLakeAdmins": [],
        "CreateDatabaseDefaultPermissions": [
            {
                "Principal": {"DataLakePrincipalIdentifier": "IAM_ALLOWED_PRINCIPALS"},
                "Permissions": ["ALL"],
            }
        ],
        "CreateTableDefaultPermissions": [
            {
                "Principal": {"DataLakePrincipalIdentifier": "IAM_ALLOWED_PRINCIPALS"},
                "Permissions": ["ALL"],
            }
        ],
        "TrustedResourceOwners": [],
        "AllowExternalDataFiltering": False,
        "ExternalDataFilteringAllowList": [],
    }


class LakeFormationBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.resources: dict[str, Resource] = {}
        self.settings: dict[str, dict[str, Any]] = defaultdict(default_settings)
        self.grants: dict[str, PermissionCatalog] = {}
        self.tagger = TaggingService()
        self.lf_database_tags: dict[tuple[str, str], list[dict[str, str]]] = {}
        self.lf_table_tags: dict[tuple[str, str, str], list[dict[str, str]]] = {}
        self.lf_columns_tags: dict[tuple[str, ...], list[dict[str, str]]] = {}

    def describe_resource(self, resource_arn: str) -> Resource:
        if resource_arn not in self.resources:
            raise EntityNotFound
        return self.resources[resource_arn]

    def deregister_resource(self, resource_arn: str) -> None:
        if resource_arn not in self.resources:
            raise EntityNotFound
        del self.resources[resource_arn]

    def register_resource(self, resource_arn: str, role_arn: str) -> None:
        if resource_arn in self.resources:
            raise AlreadyExists(
                "An error occurred (AlreadyExistsException) when calling the RegisterResource operation: Resource is already registered"
            )
        self.resources[resource_arn] = Resource(resource_arn, role_arn)

    def list_resources(self) -> list[Resource]:
        return list(self.resources.values())

    def get_data_lake_settings(self, catalog_id: str) -> dict[str, Any]:
        return self.settings[catalog_id]

    def put_data_lake_settings(self, catalog_id: str, settings: dict[str, Any]) -> None:
        self.settings[catalog_id] = settings

    def grant_permissions(
        self,
        catalog_id: str,
        principal: dict[str, str],
        resource: dict[str, Any],
        permissions: list[str],
        permissions_with_grant_options: list[str],
    ) -> None:
        if catalog_id not in self.grants:
            self.grants[catalog_id] = PermissionCatalog()

        self.grants[catalog_id].add_permission(
            Permission(
                principal=principal,
                resource=resource,
                permissions=permissions or [],
                permissions_with_grant_options=permissions_with_grant_options or [],
            )
        )

    def revoke_permissions(
        self,
        catalog_id: str,
        principal: dict[str, str],
        resource: dict[str, Any],
        permissions_to_revoke: list[str],
        permissions_with_grant_options_to_revoke: list[str],
    ) -> None:
        if catalog_id not in self.grants:
            return

        catalog = self.grants[catalog_id]
        catalog.remove_permission(
            Permission(
                principal=principal,
                resource=resource,
                permissions=permissions_to_revoke or [],
                permissions_with_grant_options=permissions_with_grant_options_to_revoke
                or [],
            )
        )

    def list_permissions(
        self,
        catalog_id: str,
        principal: Optional[dict[str, str]] = None,
        resource: Optional[ListPermissionsResource] = None,
        resource_type: Optional[RessourceType] = None,
    ) -> list[dict[str, Any]]:
        """
        No pagination has been implemented yet.
        """
        if catalog_id not in self.grants:
            return []

        permissions = list(self.grants[catalog_id].permissions)

        def filter_for_principal(permission: Permission) -> bool:
            return permission.principal == principal

        if principal is not None:
            permissions = list(filter(filter_for_principal, permissions))

        def filter_for_resource_type(permission: Permission) -> bool:
            if resource_type is None:  # Check for mypy
                return False
            resource = permission.resource
            if resource_type == RessourceType.catalog:
                return "Catalog" in resource
            elif resource_type == RessourceType.database:
                return "Database" in resource
            elif resource_type == RessourceType.data_location:
                return "DataLocation" in resource
            elif resource_type == RessourceType.table:
                return "Table" in resource or "TableWithColumns" in resource
            return False

        if resource_type is not None:
            permissions = list(filter(filter_for_resource_type, permissions))

        def filter_for_resource(permission: Permission) -> bool:
            """
            If catalog is provided:
                only matching permissions with resource-type "Catalog" are returned;
            if catalog is not provided and database is provided:
                only matching permissions with resource-type "Database" are returned;
            if catalog and database are not provided and table is provided:
                only matching permissions with resource-type "Table" are returned;
            if catalog and database and table are not provided and data location is provided:
                only matching permissions with resource-type "DataLocation" are returned;
            """
            if resource is None:  # Check for linter
                return False
            permission_resource = permission.resource
            catalog = resource.catalog
            if catalog is not None and "Catalog" in permission_resource:
                return catalog == permission_resource["Catalog"]

            database = resource.database
            if database is not None and "Database" in permission_resource:
                equals = database.name == permission_resource["Database"]["Name"]
                if database.catalog_id is not None:
                    equals = equals and (
                        database.catalog_id
                        == permission_resource["Database"].get("CatalogId")
                    )
                return equals

            table = resource.table
            if table is not None and "Table" in permission_resource:
                equals = (
                    table.database_name == permission_resource["Table"]["DatabaseName"]
                )
                if table.catalog_id is not None:
                    equals = equals and (
                        table.catalog_id
                        == permission_resource["Table"].get("CatalogId")
                    )
                if table.name is not None and table.table_wildcard is None:
                    equals = equals and (
                        table.name == permission_resource["Table"]["Name"]
                    )
                if table.name is None and table.table_wildcard is not None:
                    equals = equals and (
                        table.table_wildcard
                        == permission_resource["Table"]["TableWildcard"]
                    )
                return equals

            data_location = resource.data_location
            if data_location is not None and "DataLocation" in permission_resource:
                equals = (
                    data_location.resource_arn
                    == permission_resource["DataLocation"]["ResourceArn"]
                )
                if data_location.catalog_id is not None:
                    equals = equals and (
                        data_location.catalog_id
                        == permission_resource["DataLocation"].get("CatalogId")
                    )
                return equals

            return False

        if resource is not None:
            permissions = list(filter(filter_for_resource, permissions))

        return [permission.to_external_form() for permission in permissions]

    def create_lf_tag(self, catalog_id: str, key: str, values: list[str]) -> None:
        # There is no ARN that we can use, so just create another  unique identifier that's easy to recognize and reproduce
        arn = f"arn:lakeformation:{catalog_id}"
        tag_list = TaggingService.convert_dict_to_tags_input({key: values})  # type: ignore
        self.tagger.tag_resource(arn=arn, tags=tag_list)

    def get_lf_tag(self, catalog_id: str, key: str) -> list[str]:
        # There is no ARN that we can use, so just create another  unique identifier that's easy to recognize and reproduce
        arn = f"arn:lakeformation:{catalog_id}"
        all_tags = self.tagger.get_tag_dict_for_resource(arn=arn)
        return all_tags.get(key, [])  # type: ignore

    def delete_lf_tag(self, catalog_id: str, key: str) -> None:
        # There is no ARN that we can use, so just create another  unique identifier that's easy to recognize and reproduce
        arn = f"arn:lakeformation:{catalog_id}"
        self.tagger.untag_resource_using_names(arn, tag_names=[key])

        # Also remove any LF resource tags that used this tag-key
        for db_name in self.lf_database_tags:
            self.lf_database_tags[db_name] = [
                tag for tag in self.lf_database_tags[db_name] if tag["TagKey"] != key
            ]
        for table in self.lf_table_tags:
            self.lf_table_tags[table] = [
                tag for tag in self.lf_table_tags[table] if tag["TagKey"] != key
            ]
        for column in self.lf_columns_tags:
            self.lf_columns_tags[column] = [
                tag for tag in self.lf_columns_tags[column] if tag["TagKey"] != key
            ]

    def list_lf_tags(self, catalog_id: str) -> dict[str, str]:
        # There is no ARN that we can use, so just create another  unique identifier that's easy to recognize and reproduce
        arn = f"arn:lakeformation:{catalog_id}"
        return self.tagger.get_tag_dict_for_resource(arn=arn)

    def update_lf_tag(
        self, catalog_id: str, tag_key: str, to_delete: list[str], to_add: list[str]
    ) -> None:
        arn = f"arn:lakeformation:{catalog_id}"
        existing_tags = self.list_lf_tags(catalog_id)
        existing_tags[tag_key].extend(to_add or [])  # type: ignore
        for tag in to_delete or []:
            existing_tags[tag_key].remove(tag)  # type: ignore
        self.tagger.tag_resource(
            arn, TaggingService.convert_dict_to_tags_input(existing_tags)
        )

    def list_data_cells_filter(self) -> list[dict[str, Any]]:
        """
        This currently just returns an empty list, as the corresponding Create is not yet implemented
        """
        return []

    def batch_grant_permissions(
        self, catalog_id: str, entries: list[dict[str, Any]]
    ) -> None:
        for entry in entries:
            self.grant_permissions(
                catalog_id=catalog_id,
                principal=entry.get("Principal"),  # type: ignore[arg-type]
                resource=entry.get("Resource"),  # type: ignore[arg-type]
                permissions=entry.get("Permissions"),  # type: ignore[arg-type]
                permissions_with_grant_options=entry.get("PermissionsWithGrantOptions"),  # type: ignore[arg-type]
            )

    def batch_revoke_permissions(
        self, catalog_id: str, entries: list[dict[str, Any]]
    ) -> None:
        for entry in entries:
            self.revoke_permissions(
                catalog_id=catalog_id,
                principal=entry.get("Principal"),  # type: ignore[arg-type]
                resource=entry.get("Resource"),  # type: ignore[arg-type]
                permissions_to_revoke=entry.get("Permissions"),  # type: ignore[arg-type]
                permissions_with_grant_options_to_revoke=entry.get(  # type: ignore[arg-type]
                    "PermissionsWithGrantOptions"
                ),
            )

    def add_lf_tags_to_resource(
        self, catalog_id: str, resource: dict[str, Any], tags: list[dict[str, str]]
    ) -> list[dict[str, Any]]:
        existing_lf_tags = self.list_lf_tags(catalog_id)
        failures = []

        for tag in tags:
            if "CatalogId" not in tag:
                tag["CatalogId"] = catalog_id
                if tag["TagKey"] not in existing_lf_tags:
                    failures.append(
                        {
                            "LFTag": tag,
                            "Error": {
                                "ErrorCode": "EntityNotFoundException",
                                "ErrorMessage": "Tag or tag value does not exist.",
                            },
                        }
                    )

        if failures:
            return failures

        if "Database" in resource:
            db_catalog_id = resource["Database"].get("CatalogId", self.account_id)
            db_name = resource["Database"]["Name"]
            self.lf_database_tags[(db_catalog_id, db_name)] = tags
        if "Table" in resource:
            db_catalog_id = resource["Table"].get("CatalogId", self.account_id)
            db_name = resource["Table"]["DatabaseName"]
            name = resource["Table"]["Name"]
            self.lf_table_tags[(db_catalog_id, db_name, name)] = tags
        if "TableWithColumns" in resource:
            db_catalog_id = resource["TableWithColumns"].get(
                "CatalogId", self.account_id
            )
            db_name = resource["TableWithColumns"]["DatabaseName"]
            name = resource["TableWithColumns"]["Name"]
            for column in resource["TableWithColumns"]["ColumnNames"]:
                self.lf_columns_tags[(db_catalog_id, db_name, name, column)] = tags
        return failures

    def get_resource_lf_tags(
        self,
        catalog_id: str,
        resource: dict[str, Any],
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
        database_tags = []
        table_tags = []
        column_tags = []
        if "Database" in resource:
            database_catalog_id = resource["Database"].get("CatalogId", self.account_id)
            database_name = resource["Database"]["Name"]
            database_tags = self.lf_database_tags[(database_catalog_id, database_name)]
        if "Table" in resource:
            db_catalog_id = resource["Table"].get("CatalogId", self.account_id)
            db_name = resource["Table"]["DatabaseName"]
            name = resource["Table"]["Name"]
            table_tags = self.lf_table_tags[(db_catalog_id, db_name, name)]
        if "TableWithColumns" in resource:
            for column in resource["TableWithColumns"]["ColumnNames"]:
                db_catalog_id = resource["TableWithColumns"].get(
                    "CatalogId", self.account_id
                )
                db_name = resource["TableWithColumns"]["DatabaseName"]
                name = resource["TableWithColumns"]["Name"]
                dct_key = (db_catalog_id, db_name, name, column)
                if self.lf_columns_tags.get(dct_key):
                    column_tags.append(
                        {"Name": column, "LFTags": self.lf_columns_tags[dct_key]}
                    )
        return database_tags, table_tags, column_tags

    def remove_lf_tags_from_resource(
        self, catalog_id: str, resource: dict[str, Any], tags: list[dict[str, str]]
    ) -> None:
        for tag in tags:
            if "CatalogId" not in tag:
                tag["CatalogId"] = catalog_id
        if "Database" in resource:
            database_catalog_id = resource["Database"].get("CatalogId", self.account_id)
            database_name = resource["Database"]["Name"]
            existing_tags = self.lf_database_tags[(database_catalog_id, database_name)]
            for tag in tags:
                existing_tags.remove(tag)
        if "Table" in resource:
            db_catalog_id = resource["Table"].get("CatalogId", self.account_id)
            db_name = resource["Table"]["DatabaseName"]
            name = resource["Table"]["Name"]
            existing_tags = self.lf_table_tags[(db_catalog_id, db_name, name)]
            for tag in tags:
                existing_tags.remove(tag)
        if "TableWithColumns" in resource:
            for column in resource["TableWithColumns"]["ColumnNames"]:
                db_catalog_id = resource["TableWithColumns"].get(
                    "CatalogId", self.account_id
                )
                db_name = resource["TableWithColumns"]["DatabaseName"]
                name = resource["TableWithColumns"]["Name"]
                dct_key = (db_catalog_id, db_name, name, column)
                existing_tags = self.lf_columns_tags[dct_key]
                for tag in tags:
                    existing_tags.remove(tag)


lakeformation_backends = BackendDict(LakeFormationBackend, "lakeformation")
