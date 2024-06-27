import json
import urllib

from moto.core.responses import BaseResponse

from .models import PanoramaBackend, panorama_backends


class PanoramaResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="panorama")

    @property
    def panorama_backend(self) -> PanoramaBackend:
        return panorama_backends[self.current_account][self.region]

    def provision_device(self) -> str:
        description = self._get_param("Description")
        name = self._get_param("Name")
        networking_configuration = self._get_param("NetworkingConfiguration")
        tags = self._get_param("Tags")
        device = self.panorama_backend.provision_device(
            description=description,
            name=name,
            networking_configuration=networking_configuration,
            tags=tags,
        )
        return json.dumps(device.response_provision)

    def describe_device(self) -> str:
        device_id = urllib.parse.unquote(self._get_param("DeviceId"))
        device = self.panorama_backend.describe_device(device_id=device_id)
        return json.dumps(device.response_object())

    def list_devices(
        self,
    ) -> str:
        device_aggregated_status_filter = self._get_param(
            "DeviceAggregatedStatusFilter"
        )
        max_results = self._get_int_param("MaxResults")
        name_filter = self._get_param("NameFilter")
        next_token = self._get_param("NextToken")
        sort_by = self._get_param("SortBy")
        sort_order = self._get_param("SortOrder")
        list_devices, next_token = self.panorama_backend.list_devices(
            device_aggregated_status_filter=device_aggregated_status_filter,
            max_results=max_results,
            name_filter=name_filter,
            next_token=next_token,
            sort_by=sort_by,
            sort_order=sort_order,
        )
        return json.dumps(
            {
                "Devices": [device.response_listed() for device in list_devices],
                "NextToken": next_token,
            }
        )

    def update_device_metadata(self) -> str:
        device_id = urllib.parse.unquote(self._get_param("DeviceId"))
        description = self._get_param("Description")
        device = self.panorama_backend.update_device_metadata(
            device_id=device_id, description=description
        )
        return json.dumps(device.response_updated)

    def delete_device(self) -> str:
        device_id = urllib.parse.unquote(self._get_param("DeviceId"))
        device = self.panorama_backend.delete_device(device_id=device_id)
        return json.dumps(device.response_deleted)

    def create_node_from_template_job(self) -> str:
        job_tags = self._get_param("JobTags")
        node_description = self._get_param("NodeDescription")
        node_name = self._get_param("NodeName")
        output_package_name = self._get_param("OutputPackageName")
        output_package_version = self._get_param("OutputPackageVersion")
        template_parameters = self._get_param("TemplateParameters")
        template_type = self._get_param("TemplateType")
        node = self.panorama_backend.create_node_from_template_job(
            job_tags=job_tags,
            node_description=node_description,
            node_name=node_name,
            output_package_name=output_package_name,
            output_package_version=output_package_version,
            template_parameters=template_parameters,
            template_type=template_type,
        )
        return json.dumps(node.response_created())

    def describe_node_from_template_job(self) -> str:
        job_id = self._get_param("JobId")
        node = self.panorama_backend.describe_node_from_template_job(job_id=job_id)
        return json.dumps(node.response_described())

    def list_nodes(self) -> str:
        category = self._get_param("category")
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")
        list_nodes, next_token = self.panorama_backend.list_nodes(
            category=category, max_results=max_results, next_token=next_token
        )
        return json.dumps(
            {
                "Nodes": [node.response_listed() for node in list_nodes],
                "NextToken": next_token,
            }
        )
