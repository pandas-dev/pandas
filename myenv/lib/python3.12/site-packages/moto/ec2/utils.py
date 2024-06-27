import base64
import fnmatch
import ipaddress
import re
from typing import Any, Dict, List, Optional, Set, Tuple, TypeVar, Union

from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives.asymmetric.ed25519 import (
    Ed25519PrivateKey,
    Ed25519PublicKey,
)
from cryptography.hazmat.primitives.asymmetric.rsa import RSAPublicKey

from moto.core.utils import utcnow
from moto.iam import iam_backends
from moto.moto_api._internal import mock_random as random
from moto.utilities.utils import md5_hash

EC2_RESOURCE_TO_PREFIX = {
    "customer-gateway": "cgw",
    "transit-gateway": "tgw",
    "transit-gateway-route-table": "tgw-rtb",
    "transit-gateway-attachment": "tgw-attach",
    "dedicated_host": "h",
    "dhcp-options": "dopt",
    "fleet": "fleet",
    "flow-logs": "fl",
    "image": "ami",
    "instance": "i",
    "internet-gateway": "igw",
    "egress-only-internet-gateway": "eigw",
    "launch-template": "lt",
    "nat-gateway": "nat",
    "network-acl": "acl",
    "network-acl-subnet-assoc": "aclassoc",
    "network-interface": "eni",
    "network-interface-attachment": "eni-attach",
    "reserved-instance": "uuid4",
    "route-table": "rtb",
    "route-table-association": "rtbassoc",
    "security-group": "sg",
    "security-group-rule": "sgr",
    "snapshot": "snap",
    "spot-instance-request": "sir",
    "spot-fleet-request": "sfr",
    "subnet": "subnet",
    "subnet-ipv6-cidr-block-association": "subnet-cidr-assoc",
    "reservation": "r",
    "volume": "vol",
    "vpc": "vpc",
    "vpc-endpoint": "vpce",
    "vpc-endpoint-service": "vpce-svc",
    "managed-prefix-list": "pl",
    "vpc-cidr-association-id": "vpc-cidr-assoc",
    "vpc-elastic-ip": "eipalloc",
    "vpc-elastic-ip-association": "eipassoc",
    "vpc-peering-connection": "pcx",
    "vpn-connection": "vpn",
    "vpn-gateway": "vgw",
    "iam-instance-profile-association": "iip-assoc",
    "carrier-gateway": "cagw",
    "key-pair": "key",
}


EC2_PREFIX_TO_RESOURCE = dict((v, k) for (k, v) in EC2_RESOURCE_TO_PREFIX.items())
HEX_CHARS = list(str(x) for x in range(10)) + ["a", "b", "c", "d", "e", "f"]


def random_resource_id(size: int = 8) -> str:
    return "".join(random.choice(HEX_CHARS) for _ in range(size))


def random_id(prefix: str = "", size: int = 8) -> str:
    return f"{prefix}-{random_resource_id(size)}"


def random_ami_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["image"])


def random_instance_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["instance"], size=17)


def random_reservation_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["reservation"])


def random_security_group_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["security-group"], size=17)


def random_security_group_rule_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["security-group-rule"], size=17)


def random_fleet_id() -> str:
    return f"fleet-{random_resource_id(size=8)}-{random_resource_id(size=4)}-{random_resource_id(size=4)}-{random_resource_id(size=4)}-{random_resource_id(size=12)}"


def random_flow_log_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["flow-logs"])


def random_snapshot_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["snapshot"])


def random_spot_request_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["spot-instance-request"])


def random_spot_fleet_request_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["spot-fleet-request"])


def random_subnet_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["subnet"])


def random_subnet_ipv6_cidr_block_association_id() -> str:
    return random_id(
        prefix=EC2_RESOURCE_TO_PREFIX["subnet-ipv6-cidr-block-association"]
    )


def random_subnet_association_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["route-table-association"])


def random_network_acl_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["network-acl"])


def random_network_acl_subnet_association_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["network-acl-subnet-assoc"])


def random_vpn_gateway_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpn-gateway"])


def random_vpn_connection_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpn-connection"])


def random_customer_gateway_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["customer-gateway"])


def random_volume_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["volume"])


def random_key_pair_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["key-pair"])


def random_vpc_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpc"])


def random_vpc_ep_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpc-endpoint"], size=8)


def random_vpc_cidr_association_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpc-cidr-association-id"])


def random_vpc_peering_connection_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpc-peering-connection"])


def random_eip_association_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpc-elastic-ip-association"])


def random_internet_gateway_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["internet-gateway"])


def random_egress_only_internet_gateway_id() -> str:
    return random_id(
        prefix=EC2_RESOURCE_TO_PREFIX["egress-only-internet-gateway"], size=17
    )


def random_route_table_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["route-table"])


def random_eip_allocation_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["vpc-elastic-ip"])


def random_dhcp_option_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["dhcp-options"])


def random_eni_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["network-interface"])


def random_eni_attach_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["network-interface-attachment"])


def random_nat_gateway_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["nat-gateway"], size=17)


def random_transit_gateway_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["transit-gateway"], size=17)


def random_transit_gateway_route_table_id() -> str:
    return random_id(
        prefix=EC2_RESOURCE_TO_PREFIX["transit-gateway-route-table"], size=17
    )


def random_transit_gateway_attachment_id() -> str:
    return random_id(
        prefix=EC2_RESOURCE_TO_PREFIX["transit-gateway-attachment"], size=17
    )


def random_launch_template_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["launch-template"], size=17)


def random_launch_template_name() -> str:
    return f"LaunchTemplate_{random_resource_id(size=12)}"


def random_iam_instance_profile_association_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["iam-instance-profile-association"])


def random_carrier_gateway_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["carrier-gateway"], size=17)


def random_public_ip() -> str:
    return f"54.214.{random.choice(range(255))}.{random.choice(range(255))}"


def random_dedicated_host_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["dedicated_host"], size=17)


def random_private_ip(cidr: Optional[str] = None, ipv6: bool = False) -> str:
    # prefix - ula.prefixlen : get number of remaing length for the IP.
    #                          prefix will be 32 for IPv4 and 128 for IPv6.
    #  random.getrandbits() will generate remaining bits for IPv6 or Ipv4 in decimal format
    if cidr:
        if ipv6:
            ula = ipaddress.IPv6Network(cidr)
            return str(ula.network_address + (random.getrandbits(128 - ula.prefixlen)))
        ula = ipaddress.IPv4Network(cidr)  # type: ignore[assignment]
        return str(ula.network_address + (random.getrandbits(32 - ula.prefixlen)))
    if ipv6:
        return f"2001::cafe:{random.getrandbits(16)}x/64"
    return f"10.{random.choice(range(255))}.{random.choice(range(255))}.{random.choice(range(255))}"


def random_ip() -> str:
    return f"127.{random.randint(0, 255)}.{random.randint(0, 255)}.{random.randint(0, 255)}"


def generate_dns_from_ip(ip: Any, dns_type: str = "internal") -> str:
    splits = ip.split("/")[0].split(".") if "/" in ip else ip.split(".")
    return f"ip-{splits[0]}-{splits[1]}-{splits[2]}-{splits[3]}.ec2.{dns_type}"


def random_mac_address() -> str:
    return f"02:00:00:{random.randint(0, 255)}02x:{random.randint(0, 255)}02x:{random.randint(0, 255)}02x"


def randor_ipv4_cidr() -> str:
    return f"10.0.{random.randint(0, 255)}.{random.randint(0, 255)}/16"


def random_ipv6_cidr() -> str:
    return f"2400:6500:{random_resource_id(4)}:{random_resource_id(2)}00::/56"


def generate_route_id(
    route_table_id: str,
    cidr_block: Optional[str],
    ipv6_cidr_block: Optional[str] = None,
    prefix_list: Optional[str] = None,
) -> str:
    if ipv6_cidr_block and not cidr_block:
        cidr_block = ipv6_cidr_block
    if prefix_list and not cidr_block:
        cidr_block = prefix_list
    return f"{route_table_id}~{cidr_block}"


def random_managed_prefix_list_id() -> str:
    return random_id(prefix=EC2_RESOURCE_TO_PREFIX["managed-prefix-list"], size=8)


def create_dns_entries(service_name: str, vpc_endpoint_id: str) -> Dict[str, str]:
    return {
        "dns_name": f"{vpc_endpoint_id}-{random_resource_id(8)}.{service_name}",
        "hosted_zone_id": random_resource_id(13).upper(),
    }


def utc_date_and_time() -> str:
    x = utcnow()
    # Better performing alternative to x.strftime("%Y-%m-%dT%H:%M:%S.000Z")
    return f"{x.year}-{x.month:02d}-{x.day:02d}T{x.hour:02d}:{x.minute:02d}:{x.second:02d}.000Z"


def split_route_id(route_id: str) -> Tuple[str, str]:
    values = route_id.split("~")
    return values[0], values[1]


def get_attribute_value(
    parameter: str, querystring_dict: Dict[str, List[str]]
) -> Union[None, bool, str]:
    for key, value in querystring_dict.items():
        match = re.search(rf"{parameter}.Value", key)
        if match:
            if value[0].lower() in ["true", "false"]:
                return True if value[0].lower() in ["true"] else False
            return value[0]
    return None


def get_object_value(obj: Any, attr: str) -> Any:
    keys = attr.split(".")
    val = obj
    for key in keys:
        if hasattr(val, key):
            val = getattr(val, key)
        elif isinstance(val, dict):
            val = val[key]
        elif isinstance(val, list):
            for item in val:
                item_val = get_object_value(item, key)
                if item_val:
                    return item_val
        elif key == "owner_id" and hasattr(val, "account_id"):
            val = getattr(val, "account_id")
        else:
            return None
    return val


def is_tag_filter(filter_name: str) -> bool:
    return (
        filter_name.startswith("tag:")
        or filter_name.startswith("tag-value")
        or filter_name.startswith("tag-key")
    )


def get_obj_tag(obj: Any, filter_name: str) -> Optional[str]:
    tag_name = filter_name.replace("tag:", "", 1)
    tags = dict((tag["key"], tag["value"]) for tag in obj.get_tags())
    return tags.get(tag_name)


def get_obj_tag_names(obj: Any) -> Set[str]:
    tags = set((tag["key"] for tag in obj.get_tags()))
    return tags


def get_obj_tag_values(obj: Any, key: Optional[str] = None) -> Set[str]:
    tags = set(
        (tag["value"] for tag in obj.get_tags() if tag["key"] == key or key is None)
    )
    return tags


def add_tag_specification(tags: Any) -> Dict[str, str]:
    tags = tags[0] if isinstance(tags, list) and len(tags) == 1 else tags
    tags = (tags or {}).get("Tag", [])
    tags = {t["Key"]: t["Value"] for t in tags}
    return tags


def tag_filter_matches(obj: Any, filter_name: str, filter_values: List[str]) -> bool:
    regex_filters = [re.compile(simple_aws_filter_to_re(f)) for f in filter_values]
    if filter_name == "tag-key":
        tag_values = get_obj_tag_names(obj)
    elif filter_name == "tag-value":
        tag_values = get_obj_tag_values(obj)
    elif filter_name.startswith("tag:"):
        key = filter_name[4:]
        tag_values = get_obj_tag_values(obj, key=key)
    else:
        tag_values = [get_obj_tag(obj, filter_name) or ""]  # type: ignore[assignment]

    for tag_value in tag_values:
        if any(regex.match(tag_value) for regex in regex_filters):
            return True
        if tag_value in filter_values:
            return True

    return False


filter_dict_attribute_mapping = {
    "instance-state-name": "state",
    "instance-id": "id",
    "state-reason-code": "_state_reason.code",
    "source-dest-check": "source_dest_check",
    "vpc-id": "vpc_id",
    "group-id": "security_groups.id",
    "instance.group-id": "security_groups.id",
    "instance.group-name": "security_groups.name",
    "instance-type": "instance_type",
    "private-ip-address": "private_ip",
    "ip-address": "public_ip",
    "availability-zone": "placement",
    "architecture": "architecture",
    "image-id": "image_id",
    "network-interface.private-dns-name": "private_dns",
    "private-dns-name": "private_dns",
    "owner-id": "owner_id",
    "subnet-id": "subnet_id",
    "dns-name": "public_dns",
    "key-name": "key_name",
    "product-code": "product_codes",
}


def passes_filter_dict(instance: Any, filter_dict: Dict[str, Any]) -> bool:
    for filter_name, filter_values in filter_dict.items():
        if filter_name in filter_dict_attribute_mapping:
            instance_attr = filter_dict_attribute_mapping[filter_name]
            instance_value = get_object_value(instance, instance_attr)
            if not instance_value_in_filter_values(instance_value, filter_values):
                return False

        elif is_tag_filter(filter_name):
            if not tag_filter_matches(instance, filter_name, filter_values):
                return False
        else:
            raise NotImplementedError(
                "Filter dicts have not been implemented in Moto for '%s' yet. Feel free to open an issue at https://github.com/getmoto/moto/issues"
                % filter_name
            )
    return True


def instance_value_in_filter_values(instance_value: Any, filter_values: Any) -> bool:
    if isinstance(instance_value, list):
        if not set(filter_values).intersection(set(instance_value)):
            return False
    elif instance_value not in filter_values:
        return False
    return True


FILTER_TYPE = TypeVar("FILTER_TYPE")


def filter_reservations(
    reservations: List[FILTER_TYPE], filter_dict: Any
) -> List[FILTER_TYPE]:
    result = []
    for reservation in reservations:
        new_instances = []
        for instance in reservation.instances:  # type: ignore[attr-defined]
            if passes_filter_dict(instance, filter_dict):
                new_instances.append(instance)
        if new_instances:
            reservation.instances = new_instances  # type: ignore[attr-defined]
            result.append(reservation)
    return result


filter_dict_igw_mapping = {
    "attachment.vpc-id": "vpc.id",
    "attachment.state": "attachment_state",
    "internet-gateway-id": "id",
}


def passes_igw_filter_dict(igw: Any, filter_dict: Dict[str, Any]) -> bool:
    for filter_name, filter_values in filter_dict.items():
        if filter_name in filter_dict_igw_mapping:
            igw_attr = filter_dict_igw_mapping[filter_name]
            if get_object_value(igw, igw_attr) not in filter_values:
                return False
        elif is_tag_filter(filter_name):
            if not tag_filter_matches(igw, filter_name, filter_values):
                return False
        else:
            raise NotImplementedError(
                "Internet Gateway filter dicts have not been implemented in Moto for '%s' yet. Feel free to open an issue at https://github.com/getmoto/moto/issues",
                filter_name,
            )
    return True


def filter_internet_gateways(
    igws: List[FILTER_TYPE], filter_dict: Any
) -> List[FILTER_TYPE]:
    result = []
    for igw in igws:
        if passes_igw_filter_dict(igw, filter_dict):
            result.append(igw)
    return result


def is_filter_matching(obj: Any, _filter: str, filter_value: Any) -> bool:
    value = obj.get_filter_value(_filter)

    if filter_value is None:
        return False

    if isinstance(value, str):
        if not isinstance(filter_value, list):
            filter_value = [filter_value]
        if any(fnmatch.fnmatch(value, pattern) for pattern in filter_value):
            return True
        return False

    if isinstance(value, type({}.keys())):
        if isinstance(filter_value, str) and filter_value in value:
            return True

    try:
        value = set(value)
        return (value and value.issubset(filter_value)) or value.issuperset(
            filter_value
        )
    except TypeError:
        return value in filter_value


def generic_filter(
    filters: Dict[str, Any], objects: List[FILTER_TYPE]
) -> List[FILTER_TYPE]:
    if filters:
        for _filter, _filter_value in filters.items():
            objects = [
                obj
                for obj in objects
                if is_filter_matching(obj, _filter, _filter_value)
            ]

    return objects


def simple_aws_filter_to_re(filter_string: str) -> str:
    tmp_filter = filter_string.replace(r"\?", "[?]")
    tmp_filter = tmp_filter.replace(r"\*", "[*]")
    tmp_filter = fnmatch.translate(tmp_filter)
    return tmp_filter


def random_ed25519_key_pair() -> Dict[str, str]:
    private_key = Ed25519PrivateKey.generate()
    private_key_material = private_key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.OpenSSH,
        encryption_algorithm=serialization.NoEncryption(),
    )
    public_key = private_key.public_key()
    public_key_material = public_key.public_bytes(
        encoding=serialization.Encoding.OpenSSH,
        format=serialization.PublicFormat.OpenSSH,
    )
    fingerprint = public_key_fingerprint(public_key)

    return {
        "fingerprint": fingerprint,
        "material": private_key_material.decode("ascii"),
        "material_public": public_key_material.decode("ascii"),
    }


def random_rsa_key_pair() -> Dict[str, str]:
    private_key = rsa.generate_private_key(
        public_exponent=65537, key_size=2048, backend=default_backend()
    )
    private_key_material = private_key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.TraditionalOpenSSL,
        encryption_algorithm=serialization.NoEncryption(),
    )
    public_key = private_key.public_key()
    public_key_material = public_key.public_bytes(
        encoding=serialization.Encoding.OpenSSH,
        format=serialization.PublicFormat.OpenSSH,
    )
    fingerprint = public_key_fingerprint(public_key)

    return {
        "fingerprint": fingerprint,
        "material": private_key_material.decode("ascii"),
        "material_public": public_key_material.decode("ascii"),
    }


def get_prefix(resource_id: str) -> str:
    resource_id_prefix, _, after = resource_id.partition("-")
    if resource_id_prefix == EC2_RESOURCE_TO_PREFIX["transit-gateway"]:
        if after.startswith("rtb"):
            resource_id_prefix = EC2_RESOURCE_TO_PREFIX["transit-gateway-route-table"]
        if after.startswith("attach"):
            resource_id_prefix = EC2_RESOURCE_TO_PREFIX["transit-gateway-attachment"]
    if resource_id_prefix == EC2_RESOURCE_TO_PREFIX["network-interface"]:
        if after.startswith("attach"):
            resource_id_prefix = EC2_RESOURCE_TO_PREFIX["network-interface-attachment"]
    if resource_id.startswith(EC2_RESOURCE_TO_PREFIX["vpc-endpoint-service"]):
        resource_id_prefix = EC2_RESOURCE_TO_PREFIX["vpc-endpoint-service"]
    if resource_id_prefix not in EC2_RESOURCE_TO_PREFIX.values():
        uuid4hex = re.compile(r"[0-9a-f]{12}4[0-9a-f]{3}[89ab][0-9a-f]{15}\Z", re.I)
        if uuid4hex.match(resource_id) is not None:
            resource_id_prefix = EC2_RESOURCE_TO_PREFIX["reserved-instance"]
        else:
            # We should probably raise an error here, to make it more obvious this is not yet supported
            return None  # type: ignore[return-value]
    return resource_id_prefix


def is_valid_resource_id(resource_id: str) -> bool:
    valid_prefixes = EC2_RESOURCE_TO_PREFIX.values()
    resource_id_prefix = get_prefix(resource_id)
    if resource_id_prefix not in valid_prefixes:
        return False
    resource_id_pattern = resource_id_prefix + "-[0-9a-f]{8}"
    resource_pattern_re = re.compile(resource_id_pattern)
    return resource_pattern_re.match(resource_id) is not None


def is_valid_cidr(cird: str) -> bool:
    cidr_pattern = r"^(([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])\.){3}([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])(\/(\d|[1-2]\d|3[0-2]))$"
    cidr_pattern_re = re.compile(cidr_pattern)
    return cidr_pattern_re.match(cird) is not None


def is_valid_ipv6_cidr(cird: str) -> bool:
    cidr_pattern = r"^s*((([0-9A-Fa-f]{1,4}:){7}([0-9A-Fa-f]{1,4}|:))|(([0-9A-Fa-f]{1,4}:){6}(:[0-9A-Fa-f]{1,4}|((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3})|:))|(([0-9A-Fa-f]{1,4}:){5}(((:[0-9A-Fa-f]{1,4}){1,2})|:((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3})|:))|(([0-9A-Fa-f]{1,4}:){4}(((:[0-9A-Fa-f]{1,4}){1,3})|((:[0-9A-Fa-f]{1,4})?:((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3}))|:))|(([0-9A-Fa-f]{1,4}:){3}(((:[0-9A-Fa-f]{1,4}){1,4})|((:[0-9A-Fa-f]{1,4}){0,2}:((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3}))|:))|(([0-9A-Fa-f]{1,4}:){2}(((:[0-9A-Fa-f]{1,4}){1,5})|((:[0-9A-Fa-f]{1,4}){0,3}:((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3}))|:))|(([0-9A-Fa-f]{1,4}:){1}(((:[0-9A-Fa-f]{1,4}){1,6})|((:[0-9A-Fa-f]{1,4}){0,4}:((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3}))|:))|(:(((:[0-9A-Fa-f]{1,4}){1,7})|((:[0-9A-Fa-f]{1,4}){0,5}:((25[0-5]|2[0-4]d|1dd|[1-9]?d)(.(25[0-5]|2[0-4]d|1dd|[1-9]?d)){3}))|:)))(%.+)?s*(\/([0-9]|[1-9][0-9]|1[0-1][0-9]|12[0-8]))?$"
    cidr_pattern_re = re.compile(cidr_pattern)
    return cidr_pattern_re.match(cird) is not None


def is_valid_security_group_id(sg_id: str) -> bool:
    security_group_id_pattern = r"^sg-[a-f0-9]{8,17}$"
    compiled_re = re.compile(security_group_id_pattern)
    return compiled_re.match(sg_id) is not None


def generate_instance_identity_document(instance: Any) -> Dict[str, Any]:
    """
    http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/instance-identity-documents.html

    A JSON file that describes an instance. Usually retrieved by URL:
    http://169.254.169.254/latest/dynamic/instance-identity/document
    Here we just fill a dictionary that represents the document

    Typically, this document is used by the amazon-ecs-agent when registering a
    new ContainerInstance
    """

    document = {
        "devPayProductCodes": None,
        "availabilityZone": instance.placement["AvailabilityZone"],
        "privateIp": instance.private_ip_address,
        "version": "2010-8-31",
        "region": instance.placement["AvailabilityZone"][:-1],
        "instanceId": instance.id,
        "billingProducts": None,
        "instanceType": instance.instance_type,
        "accountId": "012345678910",
        "pendingTime": "2015-11-19T16:32:11Z",
        "imageId": instance.image_id,
        "kernelId": instance.kernel_id,
        "ramdiskId": instance.ramdisk_id,
        "architecture": instance.architecture,
    }

    return document


def _convert_rfc4716(data: bytes) -> bytes:
    """Convert an RFC 4716 public key to OpenSSH authorized_keys format"""

    # Normalize line endings and join continuation lines
    data_normalized = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    data_joined = data_normalized.replace(b"\\\n", b"")
    lines = data_joined.splitlines()

    # Trim header and footer
    if lines[0] != b"---- BEGIN SSH2 PUBLIC KEY ----":
        raise ValueError("Invalid RFC4716 header line")
    if lines[-1] != b"---- END SSH2 PUBLIC KEY ----":
        raise ValueError("Invalid RFC4716 footer line")
    lines = lines[1:-1]

    # Leading lines containing a colon are headers
    headers = {}
    num_header_lines = 0
    for line in lines:
        if b":" not in line:
            break
        num_header_lines += 1
        header_name, header_value = line.split(b": ")
        headers[header_name.lower()] = header_value

    # Remaining lines are key data
    data_lines = lines[num_header_lines:]
    b64_key = b"".join(data_lines)

    # Extract the algo name from the binary packet
    packet = base64.b64decode(b64_key)
    alg_len = int.from_bytes(packet[:4], "big")
    alg = packet[4 : 4 + alg_len]

    result_parts = [alg, b64_key]
    if b"comment" in headers:
        result_parts.append(headers[b"comment"])
    return b" ".join(result_parts)


def public_key_parse(
    key_material: Union[str, bytes],
) -> Union[RSAPublicKey, Ed25519PublicKey]:
    try:
        if isinstance(key_material, str):
            key_material = key_material.encode("ascii")
        key_material = base64.b64decode(key_material)

        if key_material.startswith(b"---- BEGIN SSH2 PUBLIC KEY ----"):
            # cryptography doesn't parse RFC4716 key format, so we have to convert it first
            key_material = _convert_rfc4716(key_material)

        public_key = serialization.load_ssh_public_key(key_material)

        if not isinstance(public_key, (RSAPublicKey, Ed25519PublicKey)):
            raise ValueError("bad key")
    except UnicodeDecodeError:
        raise ValueError("bad key")

    return public_key


def public_key_fingerprint(public_key: Union[RSAPublicKey, Ed25519PublicKey]) -> str:
    # TODO: Use different fingerprint calculation methods based on key type and source
    # see https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/verify-keys.html#how-ec2-key-fingerprints-are-calculated
    key_data = public_key.public_bytes(
        encoding=serialization.Encoding.DER,
        format=serialization.PublicFormat.SubjectPublicKeyInfo,
    )
    fingerprint_hex = md5_hash(key_data).hexdigest()
    fingerprint = re.sub(r"([a-f0-9]{2})(?!$)", r"\1:", fingerprint_hex)
    return fingerprint


def filter_iam_instance_profile_associations(
    iam_instance_associations: List[FILTER_TYPE], filter_dict: Any
) -> List[FILTER_TYPE]:
    if not filter_dict:
        return iam_instance_associations
    result = []
    for iam_instance_association in iam_instance_associations:
        filter_passed = True
        if filter_dict.get("instance-id"):
            if (
                iam_instance_association.instance.id  # type: ignore[attr-defined]
                not in filter_dict.get("instance-id").values()
            ):
                filter_passed = False
        if filter_dict.get("state"):
            if iam_instance_association.state not in filter_dict.get("state").values():  # type: ignore[attr-defined]
                filter_passed = False
        if filter_passed:
            result.append(iam_instance_association)
    return result


def filter_iam_instance_profiles(
    account_id: str,
    partition: str,
    iam_instance_profile_arn: Optional[str],
    iam_instance_profile_name: Optional[str],
) -> Any:
    instance_profile = None
    instance_profile_by_name = None
    instance_profile_by_arn = None
    backend = iam_backends[account_id][partition]
    if iam_instance_profile_name:
        instance_profile_by_name = backend.get_instance_profile(
            iam_instance_profile_name
        )
        instance_profile = instance_profile_by_name
    if iam_instance_profile_arn:
        instance_profile_by_arn = backend.get_instance_profile_by_arn(
            iam_instance_profile_arn
        )
        instance_profile = instance_profile_by_arn
    # We would prefer instance profile that we found by arn
    if iam_instance_profile_arn and iam_instance_profile_name:
        if instance_profile_by_name == instance_profile_by_arn:
            instance_profile = instance_profile_by_arn
        else:
            instance_profile = None

    return instance_profile


def describe_tag_filter(
    filters: Any, instances: List[FILTER_TYPE]
) -> List[FILTER_TYPE]:
    result = instances.copy()
    for instance in instances:
        for key in filters:
            if key.startswith("tag:"):
                match = re.match(r"tag:(.*)", key)
                if match:
                    tag_key_name = match.group(1)
                    need_delete = True
                    for tag in instance.get_tags():  # type: ignore[attr-defined]
                        if tag.get("key") == tag_key_name and tag.get(
                            "value"
                        ) in filters.get(key):
                            need_delete = False
                        elif tag.get("key") == tag_key_name and tag.get(
                            "value"
                        ) not in filters.get(key):
                            need_delete = True
                    if need_delete:
                        result.remove(instance)
    return result


def gen_moto_amis(
    described_images: List[Dict[str, Any]], drop_images_missing_keys: bool = True
) -> List[Dict[str, Any]]:
    """Convert `boto3.EC2.Client.describe_images` output to form acceptable to `MOTO_AMIS_PATH`

    Parameters
    ==========
    described_images : list of dicts
        as returned by :ref:`boto3:EC2.Client.describe_images` in "Images" key
    drop_images_missing_keys : bool, default=True
        When `True` any entry in `images` that is missing a required key will silently
        be excluded from the returned list

    Throws
    ======
    `KeyError` when `drop_images_missing_keys` is `False` and a required key is missing
    from an element of `images`

    Returns
    =======
    list of dicts suitable to be serialized into JSON as a target for `MOTO_AMIS_PATH` environment
    variable.

    See Also
    ========
    * :ref:`moto.ec2.models.EC2Backend`
    """
    result = []
    for image in described_images:
        try:
            tmp = {
                "ami_id": image["ImageId"],
                "name": image["Name"],
                "description": image["Description"],
                "owner_id": image["OwnerId"],
                "public": image["Public"],
                "virtualization_type": image["VirtualizationType"],
                "architecture": image["Architecture"],
                "state": image["State"],
                "platform": image.get("Platform"),
                "image_type": image["ImageType"],
                "hypervisor": image["Hypervisor"],
                "root_device_name": image["RootDeviceName"],
                "root_device_type": image["RootDeviceType"],
                "sriov": image.get("SriovNetSupport", "simple"),
            }
            result.append(tmp)
        except Exception as err:
            if not drop_images_missing_keys:
                raise err

    return result


def convert_tag_spec(
    tag_spec_set: List[Dict[str, Any]], tag_key: str = "Tag"
) -> Dict[str, Dict[str, str]]:
    # IN:   [{"ResourceType": _type, "Tag": [{"Key": k, "Value": v}, ..]}]
    #  (or) [{"ResourceType": _type, "Tags": [{"Key": k, "Value": v}, ..]}] <-- special cfn case
    # OUT:  {_type: {k: v, ..}}
    tags: Dict[str, Dict[str, str]] = {}
    for tag_spec in tag_spec_set:
        if tag_spec["ResourceType"] not in tags:
            tags[tag_spec["ResourceType"]] = {}
        tags[tag_spec["ResourceType"]].update(
            {tag["Key"]: tag["Value"] for tag in tag_spec[tag_key]}
        )
    return tags
