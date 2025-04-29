from collections import defaultdict
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from moto.core.common_models import BaseModel, CloudFormationModel
from moto.ec2.exceptions import InvalidParameterValueErrorTagSpotFleetRequest

if TYPE_CHECKING:
    from moto.ec2.models.instances import Instance
    from moto.ec2.models.security_groups import SecurityGroup
from ..utils import (
    convert_tag_spec,
    generic_filter,
    random_spot_fleet_request_id,
    random_spot_request_id,
)
from .core import TaggedEC2Resource
from .instance_types import INSTANCE_TYPE_OFFERINGS


class LaunchSpecification(BaseModel):
    def __init__(
        self,
        kernel_id: Optional[str],
        ramdisk_id: Optional[str],
        image_id: Optional[str],
        key_name: Optional[str],
        instance_type: str,
        placement: Optional[str],
        monitored: bool,
        subnet_id: str,
    ):
        self.key_name = key_name
        self.instance_type = instance_type
        self.image_id = image_id
        self.groups: List[SecurityGroup] = []
        self.placement = placement
        self.kernel = kernel_id
        self.ramdisk = ramdisk_id
        self.monitored = monitored
        self.subnet_id = subnet_id
        self.ebs_optimized = False


class SpotInstanceRequest(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        spot_request_id: str,
        price: str,
        image_id: str,
        spot_instance_type: str,
        valid_from: Optional[str],
        valid_until: Optional[str],
        launch_group: Optional[str],
        availability_zone_group: Optional[str],
        key_name: str,
        security_groups: List[str],
        user_data: Dict[str, Any],
        instance_type: str,
        placement: Optional[str],
        kernel_id: Optional[str],
        ramdisk_id: Optional[str],
        monitoring_enabled: bool,
        subnet_id: str,
        tags: Dict[str, Dict[str, str]],
        spot_fleet_id: Optional[str],
        instance_interruption_behaviour: Optional[str],
    ):
        super().__init__()
        self.ec2_backend = ec2_backend
        self.launch_specification = LaunchSpecification(
            kernel_id=kernel_id,
            ramdisk_id=ramdisk_id,
            image_id=image_id,
            key_name=key_name,
            instance_type=instance_type,
            placement=placement,
            monitored=monitoring_enabled,
            subnet_id=subnet_id,
        )
        self.id = spot_request_id
        self.state = "open"
        self.status = "pending-evaluation"
        self.status_message = "Your Spot request has been submitted for review, and is pending evaluation."
        if price:
            price = f"{float(price):.6f}"  # round up/down to 6 decimals
        self.price = price
        self.type = spot_instance_type
        self.valid_from = valid_from
        self.valid_until = valid_until
        self.launch_group = launch_group
        self.availability_zone_group = availability_zone_group
        self.instance_interruption_behaviour = (
            instance_interruption_behaviour or "terminate"
        )
        self.user_data = user_data  # NOT
        self.spot_fleet_id = spot_fleet_id
        tag_map = tags.get("spot-instances-request", {})
        self.add_tags(tag_map)
        self.all_tags = tags

        if security_groups:
            for group_name in security_groups:
                group = self.ec2_backend.get_security_group_by_name_or_id(group_name)
                if group:
                    self.launch_specification.groups.append(group)
        else:
            # If not security groups, add the default
            default_group = self.ec2_backend.get_security_group_by_name_or_id("default")
            self.launch_specification.groups.append(default_group)

        self.instance = self.launch_instance()
        self.state = "active"
        self.status = "fulfilled"
        self.status_message = ""

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "state":
            return self.state
        elif filter_name == "spot-instance-request-id":
            return self.id
        else:
            return super().get_filter_value(filter_name, "DescribeSpotInstanceRequests")

    def launch_instance(self) -> "Instance":
        reservation = self.ec2_backend.run_instances(
            image_id=self.launch_specification.image_id,
            count=1,
            user_data=self.user_data,
            instance_type=self.launch_specification.instance_type,
            is_instance_type_default=not self.launch_specification.instance_type,
            subnet_id=self.launch_specification.subnet_id,
            key_name=self.launch_specification.key_name,
            security_group_names=[],
            security_group_ids=self.launch_specification.groups,
            spot_fleet_id=self.spot_fleet_id,
            tags=self.all_tags,
            lifecycle="spot",
        )
        instance = reservation.instances[0]
        return instance


class SpotFleetLaunchSpec:
    def __init__(
        self,
        ebs_optimized: Any,
        group_set: List[str],
        iam_instance_profile: Any,
        image_id: str,
        instance_type: str,
        key_name: Any,
        monitoring: Any,
        spot_price: Any,
        subnet_id: Any,
        tag_specifications: Dict[str, Dict[str, str]],
        user_data: Any,
        weighted_capacity: float,
    ):
        self.ebs_optimized = ebs_optimized
        self.group_set = group_set
        self.iam_instance_profile = iam_instance_profile
        self.image_id = image_id
        self.instance_type = instance_type
        self.key_name = key_name
        self.monitoring = monitoring
        self.spot_price = spot_price
        self.subnet_id = subnet_id
        self.tag_specifications = tag_specifications
        self.user_data = user_data
        self.weighted_capacity = float(weighted_capacity)


class SpotFleetRequest(TaggedEC2Resource, CloudFormationModel):
    def __init__(
        self,
        ec2_backend: Any,
        spot_backend: "SpotRequestBackend",
        spot_fleet_request_id: str,
        spot_price: str,
        target_capacity: str,
        iam_fleet_role: str,
        allocation_strategy: str,
        launch_specs: List[Dict[str, Any]],
        launch_template_config: Optional[List[Dict[str, Any]]],
        instance_interruption_behaviour: Optional[str],
        tag_specifications: Optional[List[Dict[str, Any]]],
    ):
        self.ec2_backend = ec2_backend
        self.spot_backend = spot_backend
        self.id = spot_fleet_request_id
        self.spot_price = spot_price
        self.target_capacity = int(target_capacity)
        self.iam_fleet_role = iam_fleet_role
        self.allocation_strategy = allocation_strategy
        self.instance_interruption_behaviour = (
            instance_interruption_behaviour or "terminate"
        )
        self.state = "active"
        self.fulfilled_capacity = 0.0

        self.launch_specs = []

        self.tags = {}
        if tag_specifications is not None:
            tags = convert_tag_spec(tag_specifications)
            for resource_type in tags:
                if resource_type != "spot-fleet-request":
                    raise InvalidParameterValueErrorTagSpotFleetRequest(resource_type)
            self.tags.update(tags)

        launch_specs_from_config = []
        for config in launch_template_config or []:
            spec = config["LaunchTemplateSpecification"]
            if "LaunchTemplateId" in spec:
                launch_template = self.ec2_backend.get_launch_template(
                    template_id=spec["LaunchTemplateId"]
                )
            elif "LaunchTemplateName" in spec:
                launch_template = self.ec2_backend.get_launch_template_by_name(
                    name=spec["LaunchTemplateName"]
                )
            else:
                continue
            launch_template_data = launch_template.latest_version().data
            new_launch_template = launch_template_data.copy()
            if config.get("Overrides"):
                overrides = list(config["Overrides"].values())[0]
                new_launch_template.update(overrides)
            launch_specs_from_config.append(new_launch_template)

        for spec in (launch_specs or []) + launch_specs_from_config:
            tag_spec_set = spec.get("TagSpecificationSet", [])
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
                    spot_price=spec.get("SpotPrice", self.spot_price),
                    subnet_id=spec.get("SubnetId"),
                    tag_specifications=tags,
                    user_data=spec.get("UserData"),
                    weighted_capacity=spec.get("WeightedCapacity", 1),
                )
            )

        self.spot_requests: List[SpotInstanceRequest] = []
        self.create_spot_requests(self.target_capacity)

    @property
    def physical_resource_id(self) -> str:
        return self.id

    @staticmethod
    def cloudformation_name_type() -> str:
        return ""

    @staticmethod
    def cloudformation_type() -> str:
        # https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-ec2-spotfleet.html
        return "AWS::EC2::SpotFleet"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Any,
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "SpotFleetRequest":
        from ..models import ec2_backends

        properties = cloudformation_json["Properties"]["SpotFleetRequestConfigData"]
        ec2_backend = ec2_backends[account_id][region_name]

        spot_price = properties.get("SpotPrice")
        target_capacity = properties["TargetCapacity"]
        iam_fleet_role = properties["IamFleetRole"]
        allocation_strategy = properties["AllocationStrategy"]
        launch_specs = properties["LaunchSpecifications"]

        spot_fleet_request = ec2_backend.request_spot_fleet(
            spot_price,
            target_capacity,
            iam_fleet_role,
            allocation_strategy,
            launch_specs,
        )

        return spot_fleet_request

    def get_launch_spec_counts(
        self, weight_to_add: float
    ) -> Tuple[Dict[Any, int], float]:
        weight_map: Dict[Any, int] = defaultdict(int)

        weight_so_far = 0.0
        if self.allocation_strategy == "diversified":
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

    def create_spot_requests(self, weight_to_add: float) -> None:
        weight_map, added_weight = self.get_launch_spec_counts(weight_to_add)
        for launch_spec, count in weight_map.items():
            requests = self.spot_backend.request_spot_instances(
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

    def terminate_instances(self) -> None:
        instance_ids = []
        new_fulfilled_capacity = self.fulfilled_capacity
        for req in self.spot_requests:
            instance = req.instance
            for spec in self.launch_specs:
                if (
                    spec.instance_type == instance.instance_type
                    and spec.subnet_id == instance.subnet_id
                ):
                    break

            if new_fulfilled_capacity - spec.weighted_capacity < self.target_capacity:
                continue
            new_fulfilled_capacity -= spec.weighted_capacity  # pylint: disable=W0631
            instance_ids.append(instance.id)

        self.spot_requests = [
            req for req in self.spot_requests if req.instance.id not in instance_ids
        ]
        self.ec2_backend.terminate_instances(instance_ids)


class SpotRequestBackend:
    def __init__(self) -> None:
        self.spot_instance_requests: Dict[str, SpotInstanceRequest] = {}
        self.spot_fleet_requests: Dict[str, SpotFleetRequest] = {}

    def request_spot_instances(
        self,
        price: str,
        image_id: str,
        count: int,
        spot_instance_type: str,
        valid_from: Optional[str],
        valid_until: Optional[str],
        launch_group: Optional[str],
        availability_zone_group: Optional[str],
        key_name: str,
        security_groups: List[str],
        user_data: Dict[str, Any],
        instance_type: str,
        placement: Optional[str],
        kernel_id: Optional[str],
        ramdisk_id: Optional[str],
        monitoring_enabled: bool,
        subnet_id: str,
        tags: Optional[Dict[str, Dict[str, str]]] = None,
        spot_fleet_id: Optional[str] = None,
        instance_interruption_behaviour: Optional[str] = None,
    ) -> List[SpotInstanceRequest]:
        requests = []
        tags = tags or {}
        for _ in range(count):
            spot_request_id = random_spot_request_id()
            request = SpotInstanceRequest(
                self,
                spot_request_id,
                price,
                image_id,
                spot_instance_type,
                valid_from,
                valid_until,
                launch_group,
                availability_zone_group,
                key_name,
                security_groups,
                user_data,
                instance_type,
                placement,
                kernel_id,
                ramdisk_id,
                monitoring_enabled,
                subnet_id,
                tags,
                spot_fleet_id,
                instance_interruption_behaviour,
            )
            self.spot_instance_requests[spot_request_id] = request
            requests.append(request)
        return requests

    def describe_spot_instance_requests(
        self, filters: Any = None, spot_instance_ids: Optional[List[str]] = None
    ) -> List[SpotInstanceRequest]:
        requests = list(self.spot_instance_requests.values())

        if spot_instance_ids:
            requests = [i for i in requests if i.id in spot_instance_ids]

        return generic_filter(filters, requests)

    def cancel_spot_instance_requests(
        self, request_ids: List[str]
    ) -> List[SpotInstanceRequest]:
        requests = []
        for request_id in request_ids:
            requests.append(self.spot_instance_requests.pop(request_id))
        return requests

    def request_spot_fleet(
        self,
        spot_price: str,
        target_capacity: str,
        iam_fleet_role: str,
        allocation_strategy: str,
        launch_specs: List[Dict[str, Any]],
        launch_template_config: Optional[List[Dict[str, Any]]] = None,
        instance_interruption_behaviour: Optional[str] = None,
        tag_specifications: Optional[List[Dict[str, Any]]] = None,
    ) -> SpotFleetRequest:
        spot_fleet_request_id = random_spot_fleet_request_id()
        request = SpotFleetRequest(
            ec2_backend=self,
            spot_backend=self,
            spot_fleet_request_id=spot_fleet_request_id,
            spot_price=spot_price,
            target_capacity=target_capacity,
            iam_fleet_role=iam_fleet_role,
            allocation_strategy=allocation_strategy,
            launch_specs=launch_specs,
            launch_template_config=launch_template_config,
            instance_interruption_behaviour=instance_interruption_behaviour,
            tag_specifications=tag_specifications,
        )
        self.spot_fleet_requests[spot_fleet_request_id] = request
        return request

    def get_spot_fleet_request(
        self, spot_fleet_request_id: str
    ) -> Optional[SpotFleetRequest]:
        return self.spot_fleet_requests.get(spot_fleet_request_id)

    def describe_spot_fleet_instances(
        self, spot_fleet_request_id: str
    ) -> List[SpotInstanceRequest]:
        spot_fleet = self.get_spot_fleet_request(spot_fleet_request_id)
        if not spot_fleet:
            return []
        return spot_fleet.spot_requests

    def describe_spot_fleet_requests(
        self, spot_fleet_request_ids: List[str]
    ) -> List[SpotFleetRequest]:
        requests = list(self.spot_fleet_requests.values())

        if spot_fleet_request_ids:
            requests = [
                request for request in requests if request.id in spot_fleet_request_ids
            ]

        return requests

    def cancel_spot_fleet_requests(
        self, spot_fleet_request_ids: List[str], terminate_instances: bool
    ) -> List[SpotFleetRequest]:
        spot_requests = []
        for spot_fleet_request_id in spot_fleet_request_ids:
            spot_fleet = self.spot_fleet_requests[spot_fleet_request_id]
            if terminate_instances:
                spot_fleet.target_capacity = 0
                spot_fleet.terminate_instances()
                del self.spot_fleet_requests[spot_fleet_request_id]
            else:
                spot_fleet.state = "cancelled_running"
            spot_requests.append(spot_fleet)
        return spot_requests

    def modify_spot_fleet_request(
        self, spot_fleet_request_id: str, target_capacity: int, terminate_instances: str
    ) -> None:
        if target_capacity < 0:
            raise ValueError("Cannot reduce spot fleet capacity below 0")
        spot_fleet_request = self.spot_fleet_requests[spot_fleet_request_id]
        delta = target_capacity - spot_fleet_request.fulfilled_capacity
        spot_fleet_request.target_capacity = target_capacity
        if delta > 0:
            spot_fleet_request.create_spot_requests(delta)
        elif delta < 0 and terminate_instances == "Default":
            spot_fleet_request.terminate_instances()

    def describe_spot_price_history(
        self, instance_types: Optional[List[str]] = None, filters: Any = None
    ) -> List[Dict[str, str]]:
        matches = INSTANCE_TYPE_OFFERINGS["availability-zone"]
        matches = matches.get(self.region_name, [])  # type: ignore[attr-defined]

        def matches_filters(offering: Dict[str, Any], filters: Any) -> bool:
            def matches_filter(key: str, values: List[str]) -> bool:
                if key == "availability-zone":
                    return offering.get("Location") in values
                elif key == "instance-type":
                    return offering.get("InstanceType") in values
                else:
                    return False

            return all([matches_filter(key, values) for key, values in filters.items()])

        matches = [o for o in matches if matches_filters(o, filters)]

        if instance_types:
            matches = [t for t in matches if t.get("InstanceType") in instance_types]

        return matches
