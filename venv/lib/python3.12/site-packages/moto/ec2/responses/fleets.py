from moto.core.responses import ActionResult
from moto.core.utils import get_value

from ._base_response import EC2BaseResponse


class Fleets(EC2BaseResponse):
    def delete_fleets(self) -> ActionResult:
        fleet_ids = self._get_param("FleetIds", [])
        terminate_instances = self._get_bool_param("TerminateInstances")
        fleets = self.ec2_backend.delete_fleets(fleet_ids, terminate_instances)
        result = {
            "SuccessfulFleetDeletions": [
                {
                    "CurrentFleetState": fleet.state,
                    "PreviousFleetState": "active",
                    "FleetId": fleet.id,
                }
                for fleet in fleets
            ],
            "UnsuccessfulFleetDeletions": [],
        }
        return ActionResult(result)

    def describe_fleet_instances(self) -> ActionResult:
        fleet_id = self._get_param("FleetId")
        fleet_instances = self.ec2_backend.describe_fleet_instances(fleet_id)
        # The ec2_backend.describe_fleet_instances method returns a list of incompatible types,
        # objects and dictionaries, so we have to use `get_value` to pull from either.
        # TODO: This should be addressed on the backend.
        result = {
            "ActiveInstances": [
                {
                    "InstanceId": get_value(fi, "instance.id", None),
                    "InstanceType": get_value(fi, "instance.instance_type", None),
                    "SpotInstanceRequestId": get_value(fi, "id", None),
                    "InstanceHealth": "healthy",
                }
                for fi in fleet_instances
            ],
            "FleetId": fleet_id,
        }
        return ActionResult(result)

    def describe_fleets(self) -> ActionResult:
        fleet_ids = self._get_param("FleetIds", [])
        fleets = self.ec2_backend.describe_fleets(fleet_ids)
        result = {"Fleets": fleets}
        return ActionResult(result)

    def create_fleet(self) -> ActionResult:
        on_demand_options = self._get_param("OnDemandOptions", {})
        spot_options = self._get_param("SpotOptions", {})
        target_capacity_specification = self._get_param(
            "TargetCapacitySpecification", {}
        )
        launch_template_configs = self._get_param("LaunchTemplateConfigs", [])
        excess_capacity_termination_policy = self._get_param(
            "ExcessCapacityTerminationPolicy"
        )
        replace_unhealthy_instances = self._get_param("ReplaceUnhealthyInstances")
        terminate_instances_with_expiration = self._get_param(
            "TerminateInstancesWithExpiration", if_none=True
        )
        fleet_type = self._get_param("Type", if_none="maintain")
        valid_from = self._get_param("ValidFrom")
        valid_until = self._get_param("ValidUntil")

        tag_specifications = self._get_param("TagSpecifications", [])

        request = self.ec2_backend.create_fleet(
            on_demand_options=on_demand_options,
            spot_options=spot_options,
            target_capacity_specification=target_capacity_specification,
            launch_template_configs=launch_template_configs,
            excess_capacity_termination_policy=excess_capacity_termination_policy,
            replace_unhealthy_instances=replace_unhealthy_instances,
            terminate_instances_with_expiration=terminate_instances_with_expiration,
            fleet_type=fleet_type,
            valid_from=valid_from,
            valid_until=valid_until,
            tag_specifications=tag_specifications,
        )
        result = {"FleetId": request.id}
        if request.fleet_type == "instant":
            # On Demand and Spot Instances are stored as incompatible types on the backend,
            # so we have to do some extra work here.
            # TODO: This should be addressed on the backend.
            on_demand_instances = [
                {
                    "Lifecycle": "on-demand",
                    "InstanceIds": [instance["instance"].id],
                    "InstanceType": instance["instance"].instance_type,
                }
                for instance in request.on_demand_instances
            ]
            spot_requests = [
                {
                    "Lifecycle": "spot",
                    "InstanceIds": [instance.instance.id],
                    "InstanceType": instance.instance.instance_type,
                }
                for instance in request.spot_requests
            ]
            result["Instances"] = on_demand_instances + spot_requests
        return ActionResult(result)
