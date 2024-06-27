import json

from moto.core.responses import BaseResponse

from .models import ResourceGroupsTaggingAPIBackend, resourcegroupstaggingapi_backends


class ResourceGroupsTaggingAPIResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="resourcegroupstaggingapi")

    @property
    def backend(self) -> ResourceGroupsTaggingAPIBackend:
        return resourcegroupstaggingapi_backends[self.current_account][self.region]

    def get_resources(self) -> str:
        pagination_token = self._get_param("PaginationToken")
        tag_filters = self._get_param("TagFilters", [])
        resources_per_page = self._get_int_param("ResourcesPerPage", 50)
        tags_per_page = self._get_int_param("TagsPerPage", 100)
        resource_type_filters = self._get_param("ResourceTypeFilters", [])

        pagination_token, resource_tag_mapping_list = self.backend.get_resources(
            pagination_token=pagination_token,
            tag_filters=tag_filters,
            resources_per_page=resources_per_page,
            tags_per_page=tags_per_page,
            resource_type_filters=resource_type_filters,
        )

        # Format tag response
        response = {"ResourceTagMappingList": resource_tag_mapping_list}
        if pagination_token:
            response["PaginationToken"] = pagination_token

        return json.dumps(response)

    def get_tag_keys(self) -> str:
        pagination_token = self._get_param("PaginationToken")
        pagination_token, tag_keys = self.backend.get_tag_keys(
            pagination_token=pagination_token
        )

        response = {"TagKeys": tag_keys}
        if pagination_token:
            response["PaginationToken"] = pagination_token

        return json.dumps(response)

    def get_tag_values(self) -> str:
        pagination_token = self._get_param("PaginationToken")
        key = self._get_param("Key")
        pagination_token, tag_values = self.backend.get_tag_values(
            pagination_token=pagination_token, key=key
        )

        response = {"TagValues": tag_values}
        if pagination_token:
            response["PaginationToken"] = pagination_token

        return json.dumps(response)

    def tag_resources(self) -> str:
        resource_arns = self._get_param("ResourceARNList")
        tags = self._get_param("Tags")
        failed_resources = self.backend.tag_resources(
            resource_arns=resource_arns, tags=tags
        )

        return json.dumps({"FailedResourcesMap": failed_resources})

    # def untag_resources(self):
    #     resource_arn_list = self._get_list_prefix("ResourceARNList.member")
    #     tag_keys = self._get_list_prefix("TagKeys.member")
    #     failed_resources_map = self.backend.untag_resources(
    #         resource_arn_list=resource_arn_list,
    #         tag_keys=tag_keys,
    #     )
    #
    #     # failed_resources_map should be {'resource': {'ErrorCode': str, 'ErrorMessage': str, 'StatusCode': int}}
    #     return json.dumps({'FailedResourcesMap': failed_resources_map})
