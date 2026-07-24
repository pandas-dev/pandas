from moto.core.responses import ActionResult, BaseResponse

from .exceptions import ResourceGroupsTaggingAPIError
from .models import ResourceGroupsTaggingAPIBackend, resourcegroupstaggingapi_backends


class ResourceGroupsTaggingAPIResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="resourcegroupstaggingapi")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> ResourceGroupsTaggingAPIBackend:
        return resourcegroupstaggingapi_backends[self.current_account][self.region]

    def get_resources(self) -> ActionResult:
        pagination_token = self._get_param("PaginationToken")
        tag_filters = self._get_param("TagFilters", [])
        resources_per_page = self._get_int_param("ResourcesPerPage", 50)
        tags_per_page = self._get_int_param("TagsPerPage", 100)
        resource_type_filters = self._get_param("ResourceTypeFilters", [])
        resource_arn_list = self._get_param("ResourceARNList", [])

        # Simple range checking
        if tags_per_page not in range(100, 501):
            raise ResourceGroupsTaggingAPIError(
                "InvalidParameterException", "TagsPerPage must be between 100 and 500"
            )
        if resources_per_page not in range(1, 101):
            raise ResourceGroupsTaggingAPIError(
                "InvalidParameterException",
                "ResourcesPerPage must be between 1 and 100",
            )
        pagination_token, resource_tag_mapping_list = self.backend.get_resources(
            pagination_token=pagination_token,
            tag_filters=tag_filters,
            resources_per_page=resources_per_page,
            tags_per_page=tags_per_page,
            resource_type_filters=resource_type_filters,
            resource_arn_list=resource_arn_list,
        )

        # Format tag response
        response = {
            "ResourceTagMappingList": [
                {
                    "ResourceARN": resource.arn,
                    "Tags": [
                        {"Key": key, "Value": value}
                        for key, value in resource.tags.items()
                    ],
                }
                for resource in resource_tag_mapping_list
            ],
            "PaginationToken": pagination_token,
        }
        return ActionResult(response)

    def get_tag_keys(self) -> ActionResult:
        pagination_token = self._get_param("PaginationToken")
        pagination_token, tag_keys = self.backend.get_tag_keys(
            pagination_token=pagination_token
        )

        response = {"TagKeys": tag_keys, "PaginationToken": pagination_token}
        return ActionResult(response)

    def get_tag_values(self) -> ActionResult:
        pagination_token = self._get_param("PaginationToken")
        key = self._get_param("Key")

        pagination_token, tag_values = self.backend.get_tag_values(
            pagination_token=pagination_token, key=key
        )

        response = {"TagValues": tag_values, "PaginationToken": pagination_token}
        return ActionResult(response)

    def tag_resources(self) -> ActionResult:
        resource_arns = self._get_param("ResourceARNList", [])
        tags = self._get_param("Tags", {})
        failed_resources = self.backend.tag_resources(
            resource_arns=resource_arns, tags=tags
        )

        return ActionResult({"FailedResourcesMap": failed_resources})

    def untag_resources(self) -> ActionResult:
        resource_arn_list = self._get_param("ResourceARNList", [])
        tag_keys = self._get_param("TagKeys", [])

        failed_resources = self.backend.untag_resources(
            resource_arn_list=resource_arn_list,
            tag_keys=tag_keys,
        )

        return ActionResult({"FailedResourcesMap": failed_resources})
