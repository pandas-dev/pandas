from collections.abc import Iterator
from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.resource_tagging import TaggableResourcesMixin, TaggedResource

from ..exceptions import (
    EC2ClientError,
    InvalidID,
    MissingParameterError,
    MotoNotImplementedError,
)
from ..utils import (
    EC2_RESOURCE_TO_PREFIX,
    get_prefix,
    is_valid_resource_id,
)
from .amis import AmiBackend
from .availability_zones_and_regions import RegionsAndZonesBackend
from .carrier_gateways import CarrierGatewayBackend
from .customer_gateways import CustomerGatewayBackend
from .dhcp_options import DHCPOptionsSetBackend
from .elastic_block_store import EBSBackend
from .elastic_ip_addresses import ElasticAddressBackend
from .elastic_network_interfaces import NetworkInterfaceBackend
from .fleets import FleetsBackend
from .flow_logs import FlowLogsBackend
from .hosts import HostsBackend
from .iam_instance_profile import IamInstanceProfileAssociationBackend
from .instance_types import InstanceTypeBackend, InstanceTypeOfferingBackend
from .instances import InstanceBackend
from .internet_gateways import (
    EgressOnlyInternetGatewayBackend,
    InternetGatewayBackend,
)
from .key_pairs import KeyPairBackend
from .launch_templates import LaunchTemplateBackend
from .managed_prefixes import ManagedPrefixListBackend
from .nat_gateways import NatGatewayBackend
from .network_acls import NetworkAclBackend
from .reserved_instances import ReservedInstancesBackend
from .route_tables import RouteBackend
from .security_groups import SecurityGroupBackend
from .spot_requests import SpotRequestBackend
from .subnets import SubnetBackend
from .tags import TagBackend
from .traffic_mirrors import TrafficMirrorsBackend
from .transit_gateway import TransitGatewayBackend
from .transit_gateway_attachments import TransitGatewayAttachmentBackend
from .transit_gateway_route_tables import TransitGatewayRouteTableBackend
from .vpc_peering_connections import VPCPeeringConnectionBackend
from .vpc_service_configuration import VPCServiceConfigurationBackend
from .vpcs import VPCBackend
from .vpn_connections import VPNConnectionBackend
from .vpn_gateway import VpnGatewayBackend
from .windows import WindowsBackend


def validate_resource_ids(resource_ids: list[str]) -> bool:
    if not resource_ids:
        raise MissingParameterError(parameter="resourceIdSet")
    for resource_id in resource_ids:
        if not is_valid_resource_id(resource_id):
            raise InvalidID(resource_id=resource_id)
    return True


class SettingsBackend:
    def __init__(self) -> None:
        self.ebs_encryption_by_default = False

    def disable_ebs_encryption_by_default(self) -> None:
        ec2_backend = ec2_backends[self.account_id][self.region_name]  # type: ignore[attr-defined]
        ec2_backend.ebs_encryption_by_default = False

    def enable_ebs_encryption_by_default(self) -> None:
        ec2_backend = ec2_backends[self.account_id][self.region_name]  # type: ignore[attr-defined]
        ec2_backend.ebs_encryption_by_default = True

    def get_ebs_encryption_by_default(self) -> bool:
        ec2_backend = ec2_backends[self.account_id][self.region_name]  # type: ignore[attr-defined]
        return ec2_backend.ebs_encryption_by_default


class EC2Backend(
    BaseBackend,
    TaggableResourcesMixin,
    InstanceBackend,
    InstanceTypeBackend,
    InstanceTypeOfferingBackend,
    TagBackend,
    EBSBackend,
    RegionsAndZonesBackend,
    AmiBackend,
    SecurityGroupBackend,
    VPCBackend,
    ManagedPrefixListBackend,
    SubnetBackend,
    FlowLogsBackend,
    NetworkInterfaceBackend,
    VPNConnectionBackend,
    VPCServiceConfigurationBackend,
    VPCPeeringConnectionBackend,
    RouteBackend,
    InternetGatewayBackend,
    EgressOnlyInternetGatewayBackend,
    SpotRequestBackend,
    ElasticAddressBackend,
    KeyPairBackend,
    SettingsBackend,
    DHCPOptionsSetBackend,
    NetworkAclBackend,
    VpnGatewayBackend,
    CustomerGatewayBackend,
    NatGatewayBackend,
    TrafficMirrorsBackend,
    TransitGatewayBackend,
    TransitGatewayRouteTableBackend,
    TransitGatewayAttachmentBackend,
    LaunchTemplateBackend,
    IamInstanceProfileAssociationBackend,
    CarrierGatewayBackend,
    FleetsBackend,
    WindowsBackend,
    HostsBackend,
    ReservedInstancesBackend,
):
    """
    moto includes a limited set of AMIs in `moto/ec2/resources/amis.json`.
    Additionally, the default AMI's specified by SSM will be provided.

    If you require specific AMIs to be available during your tests, you can provide your own AMI definitions by setting the
    environment variable `MOTO_AMIS_PATH` to point to a JSON file containing definitions of the required AMIs.
    No other AMI's will be loaded if this environment variable is set.

    To create such a file, refer to `scripts/get_amis.py`

    .. note:: You must set `MOTO_AMIS_PATH` before importing moto.

    """

    SERVICE_NAMESPACE = "ec2"

    def __init__(self, region_name: str, account_id: str):
        BaseBackend.__init__(self, region_name, account_id)
        for backend in EC2Backend.__mro__:
            if backend not in [EC2Backend, BaseBackend, object]:
                backend.__init__(self)  # type: ignore

        # Default VPC exists by default, which is the current behavior
        # of EC2-VPC. See for detail:
        #
        #   docs.aws.amazon.com/AmazonVPC/latest/UserGuide/default-vpc.html
        #
        if not self.vpcs:
            vpc = self.create_default_vpc()
        else:
            # For now this is included for potential
            # backward-compatibility issues
            vpc = list(self.vpcs.values())[0]

        self.default_vpc = vpc

    # Use this to generate a proper error template response when in a response
    # handler.
    def raise_error(self, code: str, message: str) -> None:
        raise EC2ClientError(code, message)

    def raise_not_implemented_error(self, blurb: str) -> None:
        raise MotoNotImplementedError(blurb)

    def do_resources_exist(self, resource_ids: list[str]) -> bool:
        for resource_id in resource_ids:
            resource_prefix = get_prefix(resource_id)
            if resource_prefix == EC2_RESOURCE_TO_PREFIX["customer-gateway"]:
                self.get_customer_gateway(customer_gateway_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["dhcp-options"]:
                self.describe_dhcp_options(dhcp_options_ids=[resource_id])
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["image"]:
                self.describe_images(ami_ids=[resource_id])
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["instance"]:
                self.get_instance_by_id(instance_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["internet-gateway"]:
                self.describe_internet_gateways(internet_gateway_ids=[resource_id])
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["launch-template"]:
                self.get_launch_template(resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["network-acl"]:
                self.describe_network_acls()
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["network-interface"]:
                self.describe_network_interfaces(
                    filters={"network-interface-id": resource_id}
                )
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["reserved-instance"]:
                self.raise_not_implemented_error("DescribeReservedInstances")
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["route-table"]:
                self.get_route_table(route_table_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["security-group"]:
                self.describe_security_groups(group_ids=[resource_id])
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["snapshot"]:
                self.get_snapshot(snapshot_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["spot-instance-request"]:
                self.describe_spot_instance_requests(
                    filters={"spot-instance-request-id": resource_id}
                )
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["subnet"]:
                self.get_subnet(subnet_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["volume"]:
                self.get_volume(volume_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["vpc"]:
                self.get_vpc(vpc_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["vpc-endpoint-service"]:
                self.get_vpc_endpoint_service(resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["vpc-peering-connection"]:
                self.get_vpc_peering_connection(vpc_pcx_id=resource_id)
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["vpn-connection"]:
                self.describe_vpn_connections(vpn_connection_ids=[resource_id])
            elif resource_prefix == EC2_RESOURCE_TO_PREFIX["vpn-gateway"]:
                self.get_vpn_gateway(vpn_gateway_id=resource_id)
            elif (
                resource_prefix
                == EC2_RESOURCE_TO_PREFIX["iam-instance-profile-association"]
            ):
                self.describe_iam_instance_profile_associations(
                    association_ids=[resource_id]
                )
        return True

    # Resource Groups Tagging API (TaggableResourcesMixin method overrides)
    def iter_tagged_resources(self) -> Iterator[TaggedResource]:
        resource_collections: dict[str, Any] = {
            "ec2:customer-gateway": self.customer_gateways.values(),
            "ec2:flow-logs": self.flow_logs.values(),
            "ec2:image": self.amis.values(),
            "ec2:instance": (
                instance
                for reservation in self.reservations.values()
                for instance in reservation.instances
            ),
            "ec2:internet-gateway": self.internet_gateways.values(),
            "ec2:managed-prefix-lists": self.managed_prefix_lists.values(),
            "ec2:natgateway": self.nat_gateways.values(),
            "ec2:network-interface": self.enis.values(),
            "ec2:route-table": self.route_tables.values(),
            "ec2:security-group": (
                sg for vpc in self.groups.values() for sg in vpc.values()
            ),
            "ec2:snapshot": self.snapshots.values(),
            "ec2:spot-instance-request": self.spot_instance_requests.values(),
            "ec2:subnet": (
                subnet
                for subnets_by_az in self.subnets.values()
                for subnet in subnets_by_az.values()
            ),
            "ec2:transit-gateway": self.transit_gateways.values(),
            "ec2:transit-gateway-attachment": self.transit_gateway_attachments.values(),
            "ec2:volume": self.volumes.values(),
            "ec2:vpc": self.vpcs.values(),
            "ec2:vpc-peering-connection": self.vpc_pcxs.values(),
            "ec2:vpn-connection": self.vpn_connections.values(),
        }
        for resource_type, resources in resource_collections.items():
            type_part = resource_type.split(":", 1)[1]
            for resource in resources:
                yield TaggedResource(
                    arn=f"arn:{self.partition}:ec2:{self.region_name}:{self.account_id}:{type_part}/{resource.id}",
                    tags=dict(self.tags.get(resource.id, {})),
                    resource_type=resource_type,
                )

    def tag_resource(self, arn: str, tags: dict[str, str]) -> None:
        resource_id = arn.rsplit("/", 1)[-1]
        self.create_tags([resource_id], tags)

    def untag_resource(self, arn: str, tag_keys: list[str]) -> None:
        resource_id = arn.rsplit("/", 1)[-1]
        self.delete_tag_keys([resource_id], tag_keys)


ec2_backends = BackendDict(EC2Backend, "ec2")
