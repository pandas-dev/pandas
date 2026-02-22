import json
import re
from datetime import datetime
from typing import Any, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.id_generator import generate_str_id
from moto.utilities.utils import ARN_PARTITION_REGEX, get_partition

from .exceptions import BadRequestException, NotFoundException


class FakeResourceGroup(BaseModel):
    def __init__(
        self,
        account_id: str,
        region_name: str,
        name: str,
        resource_query: dict[str, str],
        description: Optional[str] = None,
        tags: Optional[dict[str, str]] = None,
        configuration: Optional[list[dict[str, Any]]] = None,
    ):
        self.errors: list[str] = []
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
        self.arn = f"arn:{get_partition(region_name)}:resource-groups:us-west-1:{account_id}:group/{name}"
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

    def _validate_resource_query(self, value: dict[str, str]) -> bool:
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

    def _validate_tags(self, value: dict[str, str]) -> bool:
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
    def resource_query(self) -> dict[str, str]:
        return self._resource_query

    @resource_query.setter
    def resource_query(self, value: dict[str, str]) -> None:
        if not self._validate_resource_query(value=value):
            self._raise_errors()
        self._resource_query = value

    @property
    def tags(self) -> dict[str, str]:
        return self._tags

    @tags.setter
    def tags(self, value: dict[str, str]) -> None:
        if not self._validate_tags(value=value):
            self._raise_errors()
        self._tags = value


class FakeTagSyncTask(BaseModel):
    def __init__(
        self,
        group_arn: str,
        group_name: str,
        role_arn: str,
        tag_key: Optional[str] = None,
        tag_value: Optional[str] = None,
        resource_query: Optional[dict[str, str]] = None,
    ):
        self.group_arn = group_arn
        self.group_name = group_name
        if tag_key is not None and tag_value is not None and resource_query is not None:
            raise BadRequestException(
                "To define a task, you can use TagKey and TagValue with non-null values, or use a ResourceQuery"
            )
        self.task_id = generate_str_id(
            resource_identifier=None,
            existing_ids=None,
            tags=None,
            length=27,
            include_digits=True,
            lower_case=True,
        )
        self.task_arn = f"{group_arn}/tag-sync-task/{self.task_id}"
        self.tag_key = tag_key
        self.tag_value = tag_value
        self.resource_query = resource_query
        self.role_arn = role_arn
        self.status = "ACTIVE"
        self.error_message = None
        self.created_at = datetime.now()

    def as_dict(self) -> dict[str, Any]:
        to_return: dict[str, Any] = {}
        to_return["GroupArn"] = self.group_arn
        to_return["GroupName"] = self.group_name
        to_return["TaskArn"] = self.task_arn
        if (
            self.resource_query is None
            and self.tag_key is not None
            and self.tag_value is not None
        ):
            to_return["TagKey"] = self.tag_key
            to_return["TagValue"] = self.tag_value
        elif (
            self.resource_query is not None
            and self.tag_key is None
            and self.tag_value is None
        ):
            to_return["ResourceQuery"] = self.resource_query
        to_return["RoleArn"] = self.role_arn
        to_return["Status"] = self.status
        to_return["CreatedAt"] = self.created_at

        return to_return


class ResourceGroups:
    def __init__(self) -> None:
        self.by_name: dict[str, FakeResourceGroup] = {}
        self.by_arn: dict[str, FakeResourceGroup] = {}

    def __contains__(self, item: str) -> bool:
        return item in self.by_name or item in self.by_arn

    def __getitem__(self, item: str) -> FakeResourceGroup:
        if item in self.by_name:
            return self.by_name[item]
        if item in self.by_arn:
            return self.by_arn[item]
        raise KeyError(f"Resource group '{item}' not found")

    def append(self, resource_group: FakeResourceGroup) -> None:
        self.by_name[resource_group.name] = resource_group
        self.by_arn[resource_group.arn] = resource_group

    def delete(self, name: str) -> FakeResourceGroup:
        group = self[name]
        del self.by_name[group.name]
        del self.by_arn[group.arn]
        return group


class ResourceGroupsBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.groups = ResourceGroups()
        self.tag_sync_tasks: dict[str, FakeTagSyncTask] = {}

    @staticmethod
    def _validate_resource_query(resource_query: dict[str, str]) -> None:
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
                rf"{ARN_PARTITION_REGEX}:cloudformation:[a-z]{{2}}-[a-z]+-[0-9]+:[0-9]+:stack/[-0-9A-z]+/[-0-9a-f]+$",
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
    def _validate_tags(tags: dict[str, str]) -> None:
        for tag in tags:
            if tag.lower().startswith("aws:"):
                raise BadRequestException("Tag keys must not start with 'aws:'")

    def create_group(
        self,
        name: str,
        resource_query: dict[str, str],
        description: Optional[str] = None,
        tags: Optional[dict[str, str]] = None,
        configuration: Optional[list[dict[str, Any]]] = None,
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
        group = self.get_group(group_name)
        return self.groups.delete(name=group.name)

    def get_group(self, group_name: str) -> FakeResourceGroup:
        try:
            group = self.groups[group_name]
        except KeyError:
            raise NotFoundException()
        return group

    def get_tag_sync_task(self, task_arn: str) -> dict[str, Any]:
        tag_sync_task = self.tag_sync_tasks[task_arn]
        return tag_sync_task.as_dict()

    def cancel_tag_sync_task(self, task_arn: str) -> None:
        del self.tag_sync_tasks[task_arn]
        return

    def get_tags(self, arn: str) -> dict[str, str]:
        return self.groups.by_arn[arn].tags

    def list_groups(self) -> dict[str, FakeResourceGroup]:
        """
        Pagination or the Filters-parameter is not yet implemented
        """
        return self.groups.by_name

    def tag(self, arn: str, tags: dict[str, str]) -> None:
        all_tags = self.groups.by_arn[arn].tags
        all_tags.update(tags)
        self._validate_tags(all_tags)
        self.groups.by_arn[arn].tags = all_tags

    def untag(self, arn: str, keys: list[str]) -> None:
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
        self, group_name: str, resource_query: dict[str, str]
    ) -> FakeResourceGroup:
        self._validate_resource_query(resource_query)
        self.groups.by_name[group_name].resource_query = resource_query
        return self.groups.by_name[group_name]

    def get_group_configuration(
        self, group_name: str
    ) -> Optional[list[dict[str, Any]]]:
        group = self.groups.by_name[group_name]
        return group.configuration

    def put_group_configuration(
        self, group_name: str, configuration: list[dict[str, Any]]
    ) -> FakeResourceGroup:
        self.groups.by_name[group_name].configuration = configuration
        return self.groups.by_name[group_name]

    def list_tag_sync_tasks(
        self,
        filters: Optional[list[dict[str, str]]] = None,
        max_results: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        tag_sync_tasks = []
        for task in self.tag_sync_tasks.values():
            tag_sync_tasks.append(task.as_dict())
        return tag_sync_tasks

    def start_tag_sync_task(
        self,
        group: str,
        role_arn: str,
        tag_key: Optional[str] = None,
        tag_value: Optional[str] = None,
        resource_query: Optional[dict[str, str]] = None,
    ) -> dict[str, Any]:
        group_arn, group_name = None, None
        group_arn_regex = r"arn:aws(-[a-z]+)*:resource-groups:[a-z]{2}(-[a-z]+)+-\d{1}:[0-9]{12}:group/([a-zA-Z0-9_\.-]{1,300}|[a-zA-Z0-9_\.-]{1,150}/[a-z0-9]{26})"
        match = re.search(group_arn_regex, group)
        if match is not None:
            group_arn = group
            split_group = re.split(r":[0-9]{12}:group/", group_arn)
            group_name = split_group[-1]
        else:
            group_name = group
            group_arn = "".join(
                [
                    "arn:",
                    get_partition(self.region_name),
                    ":resource-groups:",
                    self.region_name,
                    ":",
                    self.account_id,
                    ":",
                    group_name,
                ]
            )
        task = FakeTagSyncTask(
            group_arn, group_name, role_arn, tag_key, tag_value, resource_query
        )
        self.tag_sync_tasks[task.task_arn] = task
        task_dict = task.as_dict()
        del task_dict["Status"]
        del task_dict["CreatedAt"]
        return task_dict


resourcegroups_backends = BackendDict(ResourceGroupsBackend, "resource-groups")
