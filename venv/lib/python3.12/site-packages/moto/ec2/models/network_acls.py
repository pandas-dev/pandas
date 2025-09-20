from typing import Any, Dict, List, Optional

from ..exceptions import (
    DependencyViolationError,
    InvalidNetworkAclIdError,
    InvalidRouteTableIdError,
    NetworkAclEntryAlreadyExistsError,
)
from ..utils import (
    generic_filter,
    random_network_acl_id,
    random_network_acl_subnet_association_id,
)
from .core import TaggedEC2Resource


class NetworkAclBackend:
    def __init__(self) -> None:
        self.network_acls: Dict[str, "NetworkAcl"] = {}

    def get_network_acl(self, network_acl_id: str) -> "NetworkAcl":
        network_acl = self.network_acls.get(network_acl_id, None)
        if not network_acl:
            raise InvalidNetworkAclIdError(network_acl_id)
        return network_acl

    def create_network_acl(
        self,
        vpc_id: str,
        tags: Optional[List[Dict[str, str]]] = None,
        default: bool = False,
    ) -> "NetworkAcl":
        network_acl_id = random_network_acl_id()
        self.get_vpc(vpc_id)  # type: ignore[attr-defined]
        network_acl = NetworkAcl(self, network_acl_id, vpc_id, default)
        for tag in tags or []:
            network_acl.add_tag(tag["Key"], tag["Value"])
        self.network_acls[network_acl_id] = network_acl
        if default:
            self.add_default_entries(network_acl_id)
        return network_acl

    def add_default_entries(self, network_acl_id: str) -> None:
        default_acl_entries = [
            {"rule_number": "100", "rule_action": "allow", "egress": "true"},
            {"rule_number": "32767", "rule_action": "deny", "egress": "true"},
            {"rule_number": "100", "rule_action": "allow", "egress": "false"},
            {"rule_number": "32767", "rule_action": "deny", "egress": "false"},
        ]
        for entry in default_acl_entries:
            self.create_network_acl_entry(
                network_acl_id=network_acl_id,
                rule_number=entry["rule_number"],
                protocol="-1",
                rule_action=entry["rule_action"],
                egress=entry["egress"],
                cidr_block="0.0.0.0/0",
                icmp_code=None,
                icmp_type=None,
                port_range_from=None,
                port_range_to=None,
                ipv6_cidr_block=None,
            )

    def delete_network_acl(self, network_acl_id: str) -> "NetworkAcl":
        if any(
            network_acl.id == network_acl_id and len(network_acl.associations) > 0
            for network_acl in self.network_acls.values()
        ):
            raise DependencyViolationError(
                f"The network ACL '{network_acl_id}' has dependencies and cannot be deleted."
            )

        deleted = self.network_acls.pop(network_acl_id, None)
        if not deleted:
            raise InvalidNetworkAclIdError(network_acl_id)
        return deleted

    def create_network_acl_entry(
        self,
        network_acl_id: str,
        rule_number: str,
        protocol: str,
        rule_action: str,
        egress: str,
        cidr_block: str,
        icmp_code: Optional[int],
        icmp_type: Optional[int],
        port_range_from: Optional[int],
        port_range_to: Optional[int],
        ipv6_cidr_block: Optional[str],
    ) -> "NetworkAclEntry":
        network_acl = self.get_network_acl(network_acl_id)
        if any(
            entry.egress == egress and entry.rule_number == rule_number
            for entry in network_acl.network_acl_entries
        ):
            raise NetworkAclEntryAlreadyExistsError(rule_number)
        network_acl_entry = NetworkAclEntry(
            self,
            network_acl_id=network_acl_id,
            rule_number=rule_number,
            protocol=protocol,
            rule_action=rule_action,
            egress=egress,
            cidr_block=cidr_block,
            icmp_code=icmp_code,
            icmp_type=icmp_type,
            port_range_from=port_range_from,
            port_range_to=port_range_to,
            ipv6_cidr_block=ipv6_cidr_block,
        )

        network_acl.network_acl_entries.append(network_acl_entry)
        return network_acl_entry

    def delete_network_acl_entry(
        self, network_acl_id: str, rule_number: str, egress: str
    ) -> "NetworkAclEntry":
        network_acl = self.get_network_acl(network_acl_id)
        entry = next(
            entry
            for entry in network_acl.network_acl_entries
            if entry.egress == egress and entry.rule_number == rule_number
        )
        if entry is not None:
            network_acl.network_acl_entries.remove(entry)
        return entry

    def replace_network_acl_entry(
        self,
        network_acl_id: str,
        rule_number: str,
        protocol: str,
        rule_action: str,
        egress: str,
        cidr_block: str,
        icmp_code: int,
        icmp_type: int,
        port_range_from: int,
        port_range_to: int,
        ipv6_cidr_block: Optional[str],
    ) -> "NetworkAclEntry":
        self.delete_network_acl_entry(network_acl_id, rule_number, egress)
        network_acl_entry = self.create_network_acl_entry(
            network_acl_id=network_acl_id,
            rule_number=rule_number,
            protocol=protocol,
            rule_action=rule_action,
            egress=egress,
            cidr_block=cidr_block,
            icmp_code=icmp_code,
            icmp_type=icmp_type,
            port_range_from=port_range_from,
            port_range_to=port_range_to,
            ipv6_cidr_block=ipv6_cidr_block,
        )
        return network_acl_entry

    def replace_network_acl_association(
        self, association_id: str, network_acl_id: str
    ) -> "NetworkAclAssociation":
        # lookup existing association for subnet and delete it
        default_acl = next(
            value
            for key, value in self.network_acls.items()
            if association_id in value.associations.keys()
        )

        subnet_id = None
        for key in default_acl.associations:
            if key == association_id:
                subnet_id = default_acl.associations[key].subnet_id
                del default_acl.associations[key]
                break

        new_assoc_id = random_network_acl_subnet_association_id()
        association = NetworkAclAssociation(
            self, new_assoc_id, subnet_id, network_acl_id
        )
        new_acl = self.get_network_acl(network_acl_id)
        new_acl.associations[new_assoc_id] = association
        return association

    def associate_default_network_acl_with_subnet(
        self, subnet_id: str, vpc_id: str
    ) -> None:
        association_id = random_network_acl_subnet_association_id()
        acl = next(
            acl
            for acl in self.network_acls.values()
            if acl.default and acl.vpc_id == vpc_id
        )
        acl.associations[association_id] = NetworkAclAssociation(
            self, association_id, subnet_id, acl.id
        )

    def describe_network_acls(
        self, network_acl_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List["NetworkAcl"]:
        network_acls = list(self.network_acls.values())

        if network_acl_ids:
            network_acls = [
                network_acl
                for network_acl in network_acls
                if network_acl.id in network_acl_ids
            ]
            if len(network_acls) != len(network_acl_ids):
                invalid_id = list(
                    set(network_acl_ids).difference(
                        set([network_acl.id for network_acl in network_acls])
                    )
                )[0]
                raise InvalidRouteTableIdError(invalid_id)

        return generic_filter(filters, network_acls)


class NetworkAclAssociation:
    def __init__(
        self,
        ec2_backend: Any,
        new_association_id: str,
        subnet_id: Optional[str],
        network_acl_id: str,
    ):
        self.ec2_backend = ec2_backend
        self.id = new_association_id
        self.new_association_id = new_association_id
        self.subnet_id = subnet_id
        self.network_acl_id = network_acl_id


class NetworkAcl(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        network_acl_id: str,
        vpc_id: str,
        default: bool = False,
        owner_id: Optional[str] = None,
    ):
        self.ec2_backend = ec2_backend
        self.id = network_acl_id
        self.vpc_id = vpc_id
        self.owner_id = owner_id or ec2_backend.account_id
        self.network_acl_entries: List[NetworkAclEntry] = []
        self.associations: Dict[str, NetworkAclAssociation] = {}
        self.default = "true" if default is True else "false"

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        if filter_name == "default":
            return self.default
        elif filter_name == "vpc-id":
            return self.vpc_id
        elif filter_name == "association.network-acl-id":
            return self.id
        elif filter_name == "association.subnet-id":
            return [assoc.subnet_id for assoc in self.associations.values()]
        elif filter_name == "entry.cidr":
            return [entry.cidr_block for entry in self.network_acl_entries]
        elif filter_name == "entry.protocol":
            return [entry.protocol for entry in self.network_acl_entries]
        elif filter_name == "entry.rule-number":
            return [entry.rule_number for entry in self.network_acl_entries]
        elif filter_name == "entry.rule-action":
            return [entry.rule_action for entry in self.network_acl_entries]
        elif filter_name == "entry.egress":
            return [entry.egress for entry in self.network_acl_entries]
        elif filter_name == "owner-id":
            return self.owner_id
        else:
            return super().get_filter_value(filter_name, "DescribeNetworkAcls")


class NetworkAclEntry(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        network_acl_id: str,
        rule_number: str,
        protocol: str,
        rule_action: str,
        egress: str,
        cidr_block: str,
        icmp_code: Optional[int],
        icmp_type: Optional[int],
        port_range_from: Optional[int],
        port_range_to: Optional[int],
        ipv6_cidr_block: Optional[str],
    ):
        self.ec2_backend = ec2_backend
        self.network_acl_id = network_acl_id
        self.rule_number = rule_number
        self.protocol = protocol
        self.rule_action = rule_action
        self.egress = egress
        self.cidr_block = cidr_block
        self.ipv6_cidr_block = ipv6_cidr_block
        self.icmp_code = icmp_code
        self.icmp_type = icmp_type
        self.port_range_from = port_range_from
        self.port_range_to = port_range_to
