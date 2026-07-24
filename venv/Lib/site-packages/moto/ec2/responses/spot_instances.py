from moto.core.responses import ActionResult

from ..utils import parse_user_data
from ._base_response import EC2BaseResponse


class SpotInstances(EC2BaseResponse):
    def cancel_spot_instance_requests(self) -> ActionResult:
        request_ids = self._get_param("SpotInstanceRequestIds", [])

        self.error_on_dryrun()

        requests = self.ec2_backend.cancel_spot_instance_requests(request_ids)
        cancelled = [
            {"SpotInstanceRequestId": r.id, "State": "cancelled"} for r in requests
        ]
        return ActionResult({"CancelledSpotInstanceRequests": cancelled})

    def create_spot_datafeed_subscription(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "SpotInstances.create_spot_datafeed_subscription is not yet implemented"
        )

    def delete_spot_datafeed_subscription(self) -> None:
        self.error_on_dryrun()

        raise NotImplementedError(
            "SpotInstances.delete_spot_datafeed_subscription is not yet implemented"
        )

    def describe_spot_datafeed_subscription(self) -> None:
        raise NotImplementedError(
            "SpotInstances.describe_spot_datafeed_subscription is not yet implemented"
        )

    def describe_spot_instance_requests(self) -> ActionResult:
        spot_instance_ids = self._get_param("SpotInstanceRequestIds", [])
        filters = self._filters_from_querystring()
        requests = self.ec2_backend.describe_spot_instance_requests(
            filters=filters, spot_instance_ids=spot_instance_ids
        )
        return ActionResult({"SpotInstanceRequests": requests})

    def describe_spot_price_history(self) -> ActionResult:
        instance_types_filters = self._get_param("InstanceTypes", [])
        filter_dict = self._filters_from_querystring()
        prices = self.ec2_backend.describe_spot_price_history(
            instance_types_filters, filter_dict
        )
        history = [
            {
                "InstanceType": p["InstanceType"],
                "ProductDescription": "Linux/UNIX (Amazon VPC)",
                "SpotPrice": "0.00001",
                "AvailabilityZone": p["Location"],
                "Timestamp": "2006-01-02T15:04:05.999999999Z",
            }
            for p in prices
        ]
        return ActionResult({"SpotPriceHistory": history})

    def request_spot_instances(self) -> ActionResult:
        price = self._get_param("SpotPrice")
        image_id = self._get_param("LaunchSpecification.ImageId")
        count = self._get_int_param("InstanceCount", 1)
        spot_instance_type = self._get_param("Type", "one-time")
        valid_from = self._get_param("ValidFrom")
        valid_until = self._get_param("ValidUntil")
        launch_group = self._get_param("LaunchGroup")
        availability_zone_group = self._get_param("AvailabilityZoneGroup")
        key_name = self._get_param("LaunchSpecification.KeyName")
        security_groups = self._get_param("LaunchSpecification.SecurityGroups", [])
        user_data = parse_user_data(self._get_param("LaunchSpecification.UserData"))
        instance_type = self._get_param("LaunchSpecification.InstanceType", "m1.small")
        placement = self._get_param("LaunchSpecification.Placement.AvailabilityZone")
        kernel_id = self._get_param("LaunchSpecification.KernelId")
        ramdisk_id = self._get_param("LaunchSpecification.RamdiskId")
        monitoring_enabled = self._get_param("LaunchSpecification.Monitoring.Enabled")
        subnet_id = self._get_param("LaunchSpecification.SubnetId")
        instance_interruption_behaviour = self._get_param(
            "InstanceInterruptionBehavior"
        )
        tags = self._parse_tag_specification()

        self.error_on_dryrun()

        requests = self.ec2_backend.request_spot_instances(
            price=price,
            image_id=image_id,
            count=count,
            spot_instance_type=spot_instance_type,
            valid_from=valid_from,
            valid_until=valid_until,
            launch_group=launch_group,
            availability_zone_group=availability_zone_group,
            key_name=key_name,
            security_groups=security_groups,
            user_data=user_data,
            instance_type=instance_type,
            placement=placement,
            kernel_id=kernel_id,
            ramdisk_id=ramdisk_id,
            monitoring_enabled=monitoring_enabled,
            subnet_id=subnet_id,
            instance_interruption_behaviour=instance_interruption_behaviour,
            tags=tags,
        )

        return ActionResult({"SpotInstanceRequests": requests})
