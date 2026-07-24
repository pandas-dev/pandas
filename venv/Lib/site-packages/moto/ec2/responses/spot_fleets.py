from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class SpotFleets(EC2BaseResponse):
    def cancel_spot_fleet_requests(self) -> ActionResult:
        spot_fleet_request_ids = self._get_param("SpotFleetRequestIds", [])
        terminate_instances = self._get_bool_param("TerminateInstances")
        spot_fleets = self.ec2_backend.cancel_spot_fleet_requests(
            spot_fleet_request_ids, terminate_instances
        )
        successful = [
            {
                "SpotFleetRequestId": sf.id,
                "CurrentSpotFleetRequestState": "cancelled_terminating",
                "PreviousSpotFleetRequestState": "active",
            }
            for sf in spot_fleets
        ]
        return ActionResult({"SuccessfulFleetRequests": successful})

    def describe_spot_fleet_instances(self) -> ActionResult:
        spot_fleet_request_id = self._get_param("SpotFleetRequestId")
        spot_requests = self.ec2_backend.describe_spot_fleet_instances(
            spot_fleet_request_id
        )
        active_instances = [
            {
                "InstanceId": sr.instance.id,
                "SpotInstanceRequestId": sr.id,
                "InstanceType": sr.instance.instance_type,
            }
            for sr in spot_requests
        ]
        return ActionResult(
            {
                "SpotFleetRequestId": spot_fleet_request_id,
                "ActiveInstances": active_instances,
            }
        )

    def describe_spot_fleet_requests(self) -> ActionResult:
        spot_fleet_request_ids = self._get_param("SpotFleetRequestIds", [])
        requests = self.ec2_backend.describe_spot_fleet_requests(spot_fleet_request_ids)
        return ActionResult({"SpotFleetRequestConfigs": requests})

    def modify_spot_fleet_request(self) -> ActionResult:
        spot_fleet_request_id = self._get_param("SpotFleetRequestId")
        target_capacity = self._get_int_param("TargetCapacity")
        terminate_instances = self._get_param(
            "ExcessCapacityTerminationPolicy", if_none="Default"
        )
        self.ec2_backend.modify_spot_fleet_request(
            spot_fleet_request_id, target_capacity, terminate_instances
        )
        return ActionResult({"Return": True})

    def request_spot_fleet(self) -> ActionResult:
        spot_config = self._get_param("SpotFleetRequestConfig", {})
        spot_price = spot_config.get("SpotPrice")
        target_capacity = spot_config["TargetCapacity"]
        iam_fleet_role = spot_config["IamFleetRole"]
        allocation_strategy = spot_config["AllocationStrategy"]
        instance_interruption_behaviour = spot_config.get(
            "InstanceInterruptionBehavior"
        )

        launch_specs = spot_config.get("LaunchSpecifications")
        launch_template_config = self._get_param(
            "SpotFleetRequestConfig.LaunchTemplateConfigs", []
        )
        tag_specifications = spot_config.get("TagSpecifications", [])

        request = self.ec2_backend.request_spot_fleet(
            spot_price=spot_price,
            target_capacity=target_capacity,
            iam_fleet_role=iam_fleet_role,
            allocation_strategy=allocation_strategy,
            launch_specs=launch_specs,
            launch_template_config=launch_template_config,
            instance_interruption_behaviour=instance_interruption_behaviour,
            tag_specifications=tag_specifications,
        )
        return ActionResult({"SpotFleetRequestId": request.id})
