import json
import re
from typing import Any, Dict, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from .exceptions import BadRequestException


class FakeResourceGroup(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        resource_query: Dict[str, str],
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        configuration: Optional[List[Dict[str, Any]]] = None,
    ):
        self.errors: List[str] = []
        description = description or ""
        tags = tags or {}
        if self._validate_description(value=description):
            self._description = description
        if self._validate_name(value=name):
            self._name = name
        if self._validate_resource_query(value=resource_query):
            self._resource_query = resource_query
        if self._validate_tags(value=tags):
            self._tags = tags
        self._raise_errors()
        self.arn = f"arn:{get_partition(region_name)}:resource-groups:us-west-1:{account_id}:{name}"
        self.configuration = configuration

    @staticmethod
    def _format_error(key: str, value: Any, constraint: str) -> str:  # type: ignore[misc]
        return f"Value '{value}' at '{key}' failed to satisfy constraint: {constraint}"

    def _raise_errors(self) -> None:
        if self.errors:
            errors_len = len(self.errors)
            plural = "s" if len(self.errors) > 1 else ""
            errors = "; ".join(self.errors)
            raise BadRequestException(
                f"{errors_len} validation error{plural} detected: {errors}"
            )

    def _validate_description(self, value: str) -> bool:
        errors = []
        if len(value) > 511:
            errors.append(
                self._format_error(
                    key="description",
                    value=value,
                    constraint="Member must have length less than or equal to 512",
                )
            )
        if not re.match(r"^[\sa-zA-Z0-9_.-]*$", value):
            errors.append(
                self._format_error(
                    key="name",
                    value=value,
                    constraint=r"Member must satisfy regular expression pattern: [\sa-zA-Z0-9_\.-]*",
                )
            )
        if errors:
            self.errors += errors
            return False
        return True

    def _validate_name(self, value: str) -> bool:
        errors = []
        if len(value) > 128:
            errors.append(
                self._format_error(
                    key="name",
                    value=value,
                    constraint="Member must have length less than or equal to 128",
                )
            )
        # Note \ is a character to match not an escape.
        if not re.match(r"^[a-zA-Z0-9_\\.-]+$", value):
            errors.append(
                self._format_error(
                    key="name",
                    value=value,
                    constraint=r"Member must satisfy regular expression pattern: [a-zA-Z0-9_\.-]+",
                )
            )
        if errors:
            self.errors += errors
            return False
        return True

    def _validate_resource_query(self, value: Dict[str, str]) -> bool:
        if not value:
            return True
        errors = []
        if value["Type"] not in {"CLOUDFORMATION_STACK_1_0", "TAG_FILTERS_1_0"}:
            errors.append(
                self._format_error(
                    key="resourceQuery.type",
                    value=value,
                    constraint="Member must satisfy enum value set: [CLOUDFORMATION_STACK_1_0, TAG_FILTERS_1_0]",
                )
            )
        if len(value["Query"]) > 2048:
            errors.append(
                self._format_error(
                    key="resourceQuery.query",
                    value=value,
                    constraint="Member must have length less than or equal to 2048",
                )
            )
        if errors:
            self.errors += errors
            return False
        return True

    def _validate_tags(self, value: Dict[str, str]) -> bool:
        errors = []
        # AWS only outputs one error for all keys and one for all values.
        error_keys = None
        error_values = None
        regex = re.compile(r"^([\\p{L}\\p{Z}\\p{N}_.:/=+\-@]*)$")
        for tag_key, tag_value in value.items():
            # Validation for len(tag_key) >= 1 is done by botocore.
            if len(tag_key) > 128 or re.match(regex, tag_key):
                error_keys = self._format_error(
                    key="tags",
                    value=value,
                    constraint=(
                        "Map value must satisfy constraint: ["
                        "Member must have length less than or equal to 128, "
                        "Member must have length greater than or equal to 1, "
                        r"Member must satisfy regular expression pattern: ^([\p{L}\p{Z}\p{N}_.:/=+\-@]*)$"
                        "]"
                    ),
                )
            # Validation for len(tag_value) >= 0 is nonsensical.
            if len(tag_value) > 256 or re.match(regex, tag_key):
                error_values = self._format_error(
                    key="tags",
                    value=value,
                    constraint=(
                        "Map value must satisfy constraint: ["
                        "Member must have length less than or equal to 256, "
                        "Member must have length greater than or equal to 0, "
                        r"Member must satisfy regular expression pattern: ^([\p{L}\p{Z}\p{N}_.:/=+\-@]*)$"
                        "]"
                    ),
                )
        if error_keys:
            errors.append(error_keys)
        if error_values:
            errors.append(error_values)
        if errors:
            self.errors += errors
            return False
        return True

    @property
    def description(self) -> str:
        return self._description

    @description.setter
    def description(self, value: str) -> None:
        if not self._validate_description(value=value):
            self._raise_errors()
        self._description = value

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        if not self._validate_name(value=value):
            self._raise_errors()
        self._name = value

    @property
    def resource_query(self) -> Dict[str, str]:
        return self._resource_query

    @resource_query.setter
    def resource_query(self, value: Dict[str, str]) -> None:
        if not self._validate_resource_query(value=value):
            self._raise_errors()
        self._resource_query = value

    @property
    def tags(self) -> Dict[str, str]:
        return self._tags

    @tags.setter
    def tags(self, value: Dict[str, str]) -> None:
        if not self._validate_tags(value=value):
            self._raise_errors()
        self._tags = value


class ResourceGroups:
    def __init__(self) -> None:
        self.by_name: Dict[str, FakeResourceGroup] = {}
        self.by_arn: Dict[str, FakeResourceGroup] = {}

    def __contains__(self, item: str) -> bool:
        return item in self.by_name

    def append(self, resource_group: FakeResourceGroup) -> None:
        self.by_name[resource_group.name] = resource_group
        self.by_arn[resource_group.arn] = resource_group

    def delete(self, name: str) -> FakeResourceGroup:
        group = self.by_name[name]
        del self.by_name[name]
        del self.by_arn[group.arn]
        return group


class ResourceGroupsBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.groups = ResourceGroups()

    @staticmethod
    def _validate_resource_query(resource_query: Dict[str, str]) -> None:
        if not resource_query:
            return
        query_type = resource_query["Type"]
        query = json.loads(resource_query["Query"])
        query_keys = set(query.keys())
        invalid_json_exception = BadRequestException(
            "Invalid query: Invalid query format: check JSON syntax"
        )
        if not isinstance(query["ResourceTypeFilters"], list):
            raise invalid_json_exception
        if query_type == "CLOUDFORMATION_STACK_1_0":
            if query_keys != {"ResourceTypeFilters", "StackIdentifier"}:
                raise invalid_json_exception
            stack_identifier = query["StackIdentifier"]
            if not isinstance(stack_identifier, str):
                raise invalid_json_exception
            if not re.match(
                r"^arn:aws:cloudformation:[a-z]{2}-[a-z]+-[0-9]+:[0-9]+:stack/[-0-9A-z]+/[-0-9a-f]+$",
                stack_identifier,
            ):
                raise BadRequestException(
                    "Invalid query: Verify that the specified ARN is formatted correctly."
                )
            # Once checking other resources is implemented.
            # if stack_identifier not in self.cloudformation_backend.stacks:
            #   raise BadRequestException("Invalid query: The specified CloudFormation stack doesn't exist.")
        if query_type == "TAG_FILTERS_1_0":
            if query_keys != {"ResourceTypeFilters", "TagFilters"}:
                raise invalid_json_exception
            tag_filters = query["TagFilters"]
            if not isinstance(tag_filters, list):
                raise invalid_json_exception
            if not tag_filters or len(tag_filters) > 50:
                raise BadRequestException(
                    "Invalid query: The TagFilters list must contain between 1 and 50 elements"
                )
            for tag_filter in tag_filters:
                if not isinstance(tag_filter, dict):
                    raise invalid_json_exception
                if set(tag_filter.keys()) != {"Key", "Values"}:
                    raise invalid_json_exception
                key = tag_filter["Key"]
                if not isinstance(key, str):
                    raise invalid_json_exception
                if not key:
                    raise BadRequestException(
                        "Invalid query: The TagFilter element cannot have empty or null Key field"
                    )
                if len(key) > 128:
                    raise BadRequestException(
                        "Invalid query: The maximum length for a tag Key is 128"
                    )
                values = tag_filter["Values"]
                if not isinstance(values, list):
                    raise invalid_json_exception
                if len(values) > 20:
                    raise BadRequestException(
                        "Invalid query: The TagFilter Values list must contain between 0 and 20 elements"
                    )
                for value in values:
                    if not isinstance(value, str):
                        raise invalid_json_exception
                    if len(value) > 256:
                        raise BadRequestException(
                            "Invalid query: The maximum length for a tag Value is 256"
                        )

    @staticmethod
    def _validate_tags(tags: Dict[str, str]) -> None:
        for tag in tags:
            if tag.lower().startswith("aws:"):
                raise BadRequestException("Tag keys must not start with 'aws:'")

    def create_group(
        self,
        name: str,
        resource_query: Dict[str, str],
        description: Optional[str] = None,
        tags: Optional[Dict[str, str]] = None,
        configuration: Optional[List[Dict[str, Any]]] = None,
    ) -> FakeResourceGroup:
        tags = tags or {}
        group = FakeResourceGroup(
            account_id=self.account_id,
            region_name=self.region_name,
            name=name,
            resource_query=resource_query,
            description=description,
            tags=tags,
            configuration=configuration,
        )
        if name in self.groups:
            raise BadRequestException("Cannot create group: group already exists")
        if name.upper().startswith("AWS"):
            raise BadRequestException("Group name must not start with 'AWS'")
        self._validate_tags(tags)
        self._validate_resource_query(resource_query)
        self.groups.append(group)
        return group

    def delete_group(self, group_name: str) -> FakeResourceGroup:
        return self.groups.delete(name=group_name)

    def get_group(self, group_name: str) -> FakeResourceGroup:
        return self.groups.by_name[group_name]

    def get_tags(self, arn: str) -> Dict[str, str]:
        return self.groups.by_arn[arn].tags

    def list_groups(self) -> Dict[str, FakeResourceGroup]:
        """
        Pagination or the Filters-parameter is not yet implemented
        """
        return self.groups.by_name

    def tag(self, arn: str, tags: Dict[str, str]) -> None:
        all_tags = self.groups.by_arn[arn].tags
        all_tags.update(tags)
        self._validate_tags(all_tags)
        self.groups.by_arn[arn].tags = all_tags

    def untag(self, arn: str, keys: List[str]) -> None:
        group = self.groups.by_arn[arn]
        for key in keys:
            del group.tags[key]

    def update_group(
        self, group_name: str, description: Optional[str] = None
    ) -> FakeResourceGroup:
        if description:
            self.groups.by_name[group_name].description = description
        return self.groups.by_name[group_name]

    def update_group_query(
        self, group_name: str, resource_query: Dict[str, str]
    ) -> FakeResourceGroup:
        self._validate_resource_query(resource_query)
        self.groups.by_name[group_name].resource_query = resource_query
        return self.groups.by_name[group_name]

    def get_group_configuration(
        self, group_name: str
    ) -> Optional[List[Dict[str, Any]]]:
        group = self.groups.by_name[group_name]
        return group.configuration

    def put_group_configuration(
        self, group_name: str, configuration: List[Dict[str, Any]]
    ) -> FakeResourceGroup:
        self.groups.by_name[group_name].configuration = configuration
        return self.groups.by_name[group_name]


resourcegroups_backends = BackendDict(ResourceGroupsBackend, "resource-groups")
