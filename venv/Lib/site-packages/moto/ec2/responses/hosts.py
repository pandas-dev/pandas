from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class HostsResponse(EC2BaseResponse):
    def allocate_hosts(self) -> ActionResult:
        quantity = self._get_int_param("Quantity")
        host_recovery = self._get_param("HostRecovery")
        zone = self._get_param("AvailabilityZone")
        instance_type = self._get_param("InstanceType")
        instance_family = self._get_param("InstanceFamily")
        auto_placement = self._get_param("AutoPlacement")
        tags = self._parse_tag_specification()
        host_tags = tags.get("dedicated-host", {})
        host_ids = self.ec2_backend.allocate_hosts(
            quantity,
            host_recovery,
            zone,
            instance_type,
            instance_family,
            auto_placement,
            host_tags,
        )
        return ActionResult({"HostIds": host_ids})

    def describe_hosts(self) -> ActionResult:
        host_ids = self._get_param("HostIds", [])
        filters = self._filters_from_querystring()
        hosts = self.ec2_backend.describe_hosts(host_ids, filters)
        return ActionResult({"Hosts": hosts})

    def modify_hosts(self) -> ActionResult:
        host_ids = self._get_param("HostIds", [])
        auto_placement = self._get_param("AutoPlacement")
        host_recovery = self._get_param("HostRecovery")
        instance_type = self._get_param("InstanceType")
        instance_family = self._get_param("InstanceFamily")
        self.ec2_backend.modify_hosts(
            host_ids, auto_placement, host_recovery, instance_type, instance_family
        )
        return ActionResult({"Successful": host_ids})

    def release_hosts(self) -> ActionResult:
        host_ids = self._get_param("HostIds", [])
        self.ec2_backend.release_hosts(host_ids)
        return ActionResult({"Successful": host_ids})
