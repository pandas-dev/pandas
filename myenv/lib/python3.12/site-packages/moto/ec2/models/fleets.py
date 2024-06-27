import datetime
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from moto.ec2.models.spot_requests import SpotFleetLaunchSpec, SpotInstanceRequest

from ..utils import (
    convert_tag_spec,
    random_fleet_id,
)
from .core import TaggedEC2Resource


class Fleet(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        fleet_id: str,
        on_demand_options: Dict[str, Any],
        spot_options: Dict[str, Any],
        target_capacity_specification: Dict[str, Any],
        launch_template_configs: List[Dict[str, Any]],
        excess_capacity_termination_policy: str,
        replace_unhealthy_instances: bool,
        terminate_instances_with_expiration: bool,
        fleet_type: str,
        valid_from: str,
        valid_until: str,
        tag_specifications: List[Dict[str, Any]],
    ):
        self.ec2_backend = ec2_backend
        self.id = fleet_id
        self.spot_options = spot_options
        self.on_demand_options = on_demand_options
        self.target_capacity_specification = target_capacity_specification
        self.launch_template_configs = launch_template_configs
        self.excess_capacity_termination_policy = (
            excess_capacity_termination_policy or "termination"
        )
        self.replace_unhealthy_instances = replace_unhealthy_instances
        self.terminate_instances_with_expiration = terminate_instances_with_expiration
        self.fleet_type = fleet_type
        self.valid_from = valid_from or datetime.datetime.now(tz=datetime.timezone.utc)
        self.valid_until = valid_until
        tag_map = convert_tag_spec(tag_specifications).get("fleet", {})
        self.add_tags(tag_map)
        self.tags = self.get_tags()

        self.state = "active"
        self.fulfilled_capacity = 0.0
        self.fulfilled_on_demand_capacity = 0.0
        self.fulfilled_spot_capacity = 0.0

        self.launch_specs: List[SpotFleetLaunchSpec] = []

        launch_specs_from_config: List[Dict[str, Any]] = []
        for config in launch_template_configs or []:
            launch_spec = config["LaunchTemplateSpecification"]
            if "LaunchTemplateId" in launch_spec:
                launch_template = self.ec2_backend.get_launch_template(
                    template_id=launch_spec["LaunchTemplateId"]
                )
            elif "LaunchTemplateName" in launch_spec:
                launch_template = self.ec2_backend.get_launch_template_by_name(
                    name=launch_spec["LaunchTemplateName"]
                )
            else:
                continue
            launch_template_data = launch_template.latest_version().data
            new_launch_template = launch_template_data.copy()
            if config.get("Overrides"):
                for override in config["Overrides"]:
                    new_launch_template.update(override)
            launch_specs_from_config.append(new_launch_template)

        for spec in launch_specs_from_config:
            tag_spec_set = spec.get("TagSpecification", [])
            tags = convert_tag_spec(tag_spec_set)
            self.launch_specs.append(
                SpotFleetLaunchSpec(
                    ebs_optimized=spec.get("EbsOptimized"),
                    group_set=spec.get("GroupSet", []),
                    iam_instance_profile=spec.get("IamInstanceProfile"),
                    image_id=spec["ImageId"],
                    instance_type=spec["InstanceType"],
                    key_name=spec.get("KeyName"),
                    monitoring=spec.get("Monitoring"),
                    spot_price=spec.get("SpotPrice"),
                    subnet_id=spec.get("SubnetId"),
                    tag_specifications=tags,
                    user_data=spec.get("UserData"),
                    weighted_capacity=spec.get("WeightedCapacity", 1),
                )
            )

        self.spot_requests: List[SpotInstanceRequest] = []
        self.on_demand_instances: List[Dict[str, Any]] = []
        default_capacity = (
            target_capacity_specification.get("DefaultTargetCapacityType")
            or "on-demand"
        )
        self.target_capacity = int(
            target_capacity_specification.get("TotalTargetCapacity")  # type: ignore[arg-type]
        )
        self.spot_target_capacity = int(
            target_capacity_specification.get("SpotTargetCapacity", 0)
        )
        if self.spot_target_capacity > 0:
            self.create_spot_requests(self.spot_target_capacity)
        self.on_demand_target_capacity = int(
            target_capacity_specification.get(
                "OnDemandTargetCapacity", self.target_capacity
            )
        )
        if self.on_demand_target_capacity > 0:
            self.create_on_demand_requests(self.on_demand_target_capacity)

        remaining_capacity = self.target_capacity - self.fulfilled_capacity
        if remaining_capacity > 0:
            if default_capacity == "on-demand":
                self.create_on_demand_requests(remaining_capacity)
            elif default_capacity == "spot":
                self.create_spot_requests(remaining_capacity)

    @property
    def physical_resource_id(self) -> str:
        return self.id

    def create_spot_requests(self, weight_to_add: float) -> List[SpotInstanceRequest]:
        weight_map, added_weight = self.get_launch_spec_counts(weight_to_add)
        for launch_spec, count in weight_map.items():
            requests = self.ec2_backend.request_spot_instances(
                price=launch_spec.spot_price,
                image_id=launch_spec.image_id,
                count=count,
                spot_instance_type="persistent",
                valid_from=None,
                valid_until=None,
                launch_group=None,
                availability_zone_group=None,
                key_name=launch_spec.key_name,
                security_groups=launch_spec.group_set,
                user_data=launch_spec.user_data,
                instance_type=launch_spec.instance_type,
                placement=None,
                kernel_id=None,
                ramdisk_id=None,
                monitoring_enabled=launch_spec.monitoring,
                subnet_id=launch_spec.subnet_id,
                spot_fleet_id=self.id,
                tags=launch_spec.tag_specifications,
            )
            self.spot_requests.extend(requests)
        self.fulfilled_capacity += added_weight
        return self.spot_requests

    def create_on_demand_requests(self, weight_to_add: float) -> None:
        weight_map, added_weight = self.get_launch_spec_counts(weight_to_add)
        for launch_spec, count in weight_map.items():
            reservation = self.ec2_backend.run_instances(
                image_id=launch_spec.image_id,
                count=count,
                instance_type=launch_spec.instance_type,
                is_instance_type_default=False,
                key_name=launch_spec.key_name,
                security_group_names=launch_spec.group_set,
                user_data=launch_spec.user_data,
                placement=None,
                kernel_id=None,
                ramdisk_id=None,
                monitoring_enabled=launch_spec.monitoring,
                subnet_id=launch_spec.subnet_id,
                fleet_id=self.id,
                tags=launch_spec.tag_specifications,
            )

            # get the instance from the reservation
            instance = reservation.instances[0]
            self.on_demand_instances.append(
                {
                    "id": reservation.id,
                    "instance": instance,
                }
            )
        self.fulfilled_capacity += added_weight

    def get_launch_spec_counts(
        self, weight_to_add: float
    ) -> Tuple[Dict[SpotFleetLaunchSpec, int], float]:
        weight_map: Dict[SpotFleetLaunchSpec, int] = defaultdict(int)

        weight_so_far = 0.0
        if (
            self.spot_options
            and self.spot_options["AllocationStrategy"] == "diversified"
        ):
            launch_spec_index = 0
            while True:
                launch_spec = self.launch_specs[
                    launch_spec_index % len(self.launch_specs)
                ]
                weight_map[launch_spec] += 1
                weight_so_far += launch_spec.weighted_capacity
                if weight_so_far >= weight_to_add:
                    break
                launch_spec_index += 1
        else:  # lowestPrice
            cheapest_spec = sorted(
                # FIXME: change `+inf` to the on demand price scaled to weighted capacity when it's not present
                self.launch_specs,
                key=lambda spec: float(spec.spot_price or "+inf"),
            )[0]
            weight_so_far = weight_to_add + (
                weight_to_add % cheapest_spec.weighted_capacity
            )
            weight_map[cheapest_spec] = int(
                weight_so_far // cheapest_spec.weighted_capacity
            )

        return weight_map, weight_so_far

    def terminate_instances(self) -> None:
        instance_ids = []
        new_fulfilled_capacity = self.fulfilled_capacity
        for req in self.spot_requests + self.on_demand_instances:
            instance = None
            try:
                instance = req.instance  # type: ignore
            except AttributeError:
                instance = req["instance"]  # type: ignore[index]

            if instance.state == "terminated":
                continue

            # stop when we hit the target capacity
            if new_fulfilled_capacity <= self.target_capacity:
                break

            instance_ids.append(instance.id)
            new_fulfilled_capacity -= 1

        self.spot_requests = [
            req for req in self.spot_requests if req.instance.id not in instance_ids
        ]
        self.on_demand_instances = [
            req
            for req in self.on_demand_instances
            if req["instance"].id not in instance_ids
        ]
        self.ec2_backend.terminate_instances(instance_ids)


class FleetsBackend:
    def __init__(self) -> None:
        self.fleets: Dict[str, Fleet] = {}

    def create_fleet(
        self,
        on_demand_options: Dict[str, Any],
        spot_options: Dict[str, Any],
        target_capacity_specification: Dict[str, Any],
        launch_template_configs: List[Dict[str, Any]],
        excess_capacity_termination_policy: str,
        replace_unhealthy_instances: bool,
        terminate_instances_with_expiration: bool,
        fleet_type: str,
        valid_from: str,
        valid_until: str,
        tag_specifications: List[Dict[str, Any]],
    ) -> Fleet:
        fleet_id = random_fleet_id()
        fleet = Fleet(
            self,
            fleet_id,
            on_demand_options,
            spot_options,
            target_capacity_specification,
            launch_template_configs,
            excess_capacity_termination_policy,
            replace_unhealthy_instances,
            terminate_instances_with_expiration,
            fleet_type,
            valid_from,
            valid_until,
            tag_specifications,
        )
        self.fleets[fleet_id] = fleet
        return fleet

    def get_fleet(self, fleet_id: str) -> Optional[Fleet]:
        return self.fleets.get(fleet_id)

    def describe_fleet_instances(self, fleet_id: str) -> List[Any]:
        fleet = self.get_fleet(fleet_id)
        if not fleet:
            return []
        return fleet.spot_requests + fleet.on_demand_instances

    def describe_fleets(self, fleet_ids: Optional[List[str]]) -> List[Fleet]:
        fleets = list(self.fleets.values())

        if fleet_ids:
            fleets = [fleet for fleet in fleets if fleet.id in fleet_ids]

        return fleets

    def delete_fleets(
        self, fleet_ids: List[str], terminate_instances: bool
    ) -> List[Fleet]:
        fleets = []
        for fleet_id in fleet_ids:
            fleet = self.fleets[fleet_id]
            if terminate_instances:
                # State indicates the fleet is in the process of being terminated
                # AWS will change the state to `deleted` after a few seconds/minutes
                # Note that it will stay in the `deleted`-state for at least a few hours
                fleet.state = "deleted_terminating"
                fleet.target_capacity = 0
                fleet.terminate_instances()
            else:
                # State is different, and indicates that instances are still running
                fleet.state = "deleted_running"
            fleets.append(fleet)
        return fleets
