import itertools
from typing import Any, Dict, List, Optional

from ..exceptions import (
    DependencyViolationError,
    InvalidDHCPOptionsIdError,
    InvalidParameterValueError,
    MalformedDHCPOptionsIdError,
)
from ..utils import generic_filter, random_dhcp_option_id
from .core import TaggedEC2Resource


class DHCPOptionsSet(TaggedEC2Resource):
    def __init__(
        self,
        ec2_backend: Any,
        domain_name_servers: Optional[List[str]] = None,
        domain_name: Optional[str] = None,
        ntp_servers: Optional[List[str]] = None,
        netbios_name_servers: Optional[List[str]] = None,
        netbios_node_type: Optional[str] = None,
    ):
        self.ec2_backend = ec2_backend
        self._options = {
            "domain-name-servers": domain_name_servers,
            "domain-name": domain_name,
            "ntp-servers": ntp_servers,
            "netbios-name-servers": netbios_name_servers,
            "netbios-node-type": netbios_node_type,
        }
        self.id = random_dhcp_option_id()
        self.vpc = None

    def get_filter_value(
        self, filter_name: str, method_name: Optional[str] = None
    ) -> Any:
        """
        API Version 2015-10-01 defines the following filters for DescribeDhcpOptions:

        * dhcp-options-id
        * key
        * value
        * tag:key=value
        * tag-key
        * tag-value

        Taken from: http://docs.aws.amazon.com/AWSEC2/latest/APIReference/API_DescribeDhcpOptions.html
        """
        if filter_name == "dhcp-options-id":
            return self.id
        elif filter_name == "key":
            return list(self._options.keys())
        elif filter_name == "value":
            values = [item for item in list(self._options.values()) if item]
            return itertools.chain(*values)
        else:
            return super().get_filter_value(filter_name, "DescribeDhcpOptions")

    @property
    def options(self) -> Dict[str, Any]:  # type: ignore[misc]
        return self._options


class DHCPOptionsSetBackend:
    def __init__(self) -> None:
        self.dhcp_options_sets: Dict[str, DHCPOptionsSet] = {}

    def associate_dhcp_options(self, dhcp_options: DHCPOptionsSet, vpc: Any) -> None:
        dhcp_options.vpc = vpc
        vpc.dhcp_options = dhcp_options

    def disassociate_dhcp_options(self, vpc: Any) -> None:
        if vpc.dhcp_options:
            vpc.dhcp_options.vpc = None
        vpc.dhcp_options = None

    def create_dhcp_options(
        self,
        domain_name_servers: Optional[List[str]] = None,
        domain_name: Optional[str] = None,
        ntp_servers: Optional[List[str]] = None,
        netbios_name_servers: Optional[List[str]] = None,
        netbios_node_type: Optional[str] = None,
    ) -> DHCPOptionsSet:
        NETBIOS_NODE_TYPES = [1, 2, 4, 8]

        for field_value in domain_name_servers, ntp_servers, netbios_name_servers:
            if field_value and len(field_value) > 4:
                raise InvalidParameterValueError(",".join(field_value))

        if netbios_node_type and int(netbios_node_type[0]) not in NETBIOS_NODE_TYPES:
            raise InvalidParameterValueError(netbios_node_type)

        options = DHCPOptionsSet(
            self,
            domain_name_servers,
            domain_name,
            ntp_servers,
            netbios_name_servers,
            netbios_node_type,
        )
        self.dhcp_options_sets[options.id] = options
        return options

    def delete_dhcp_options_set(self, options_id: Optional[str]) -> None:
        if not (options_id and options_id.startswith("dopt-")):
            raise MalformedDHCPOptionsIdError(options_id)

        if options_id in self.dhcp_options_sets:
            if self.dhcp_options_sets[options_id].vpc:
                raise DependencyViolationError("Cannot delete assigned DHCP options.")
            self.dhcp_options_sets.pop(options_id)
        else:
            raise InvalidDHCPOptionsIdError(options_id)

    def describe_dhcp_options(
        self, dhcp_options_ids: Optional[List[str]] = None, filters: Any = None
    ) -> List[DHCPOptionsSet]:
        dhcp_options_sets = list(self.dhcp_options_sets.copy().values())

        if dhcp_options_ids:
            dhcp_options_sets = [
                dhcp_options_set
                for dhcp_options_set in dhcp_options_sets
                if dhcp_options_set.id in dhcp_options_ids
            ]
            if len(dhcp_options_sets) != len(dhcp_options_ids):
                invalid_id = list(
                    set(dhcp_options_ids).difference(
                        set(
                            [
                                dhcp_options_set.id
                                for dhcp_options_set in dhcp_options_sets
                            ]
                        )
                    )
                )[0]
                raise InvalidDHCPOptionsIdError(invalid_id)

        return generic_filter(filters, dhcp_options_sets)
