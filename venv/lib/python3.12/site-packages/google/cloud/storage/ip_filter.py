# Copyright 2014 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""IP Filter configuration for Google Cloud Storage Buckets."""

from typing import Dict, Any, Optional, List

_MODE = "mode"
_PUBLIC_NETWORK_SOURCE = "publicNetworkSource"
_VPC_NETWORK_SOURCES = "vpcNetworkSources"
_ALLOWED_IP_CIDR_RANGES = "allowedIpCidrRanges"
_NETWORK = "network"
_ALLOW_ALL_SERVICE_AGENT_ACCESS = "allowAllServiceAgentAccess"
_ALLOW_CROSS_ORG_VPCS = "allowCrossOrgVpcs"


class PublicNetworkSource:
    """Represents a public network source for a GCS Bucket IP Filter.

    :type allowed_ip_cidr_ranges: list(str) or None
    :param allowed_ip_cidr_ranges: A list of public IPv4 or IPv6 ranges in
                                   CIDR notation that are allowed to access
                                   the bucket.
    """

    def __init__(self, allowed_ip_cidr_ranges: Optional[List[str]] = None):
        self.allowed_ip_cidr_ranges = allowed_ip_cidr_ranges or []

    def _to_api_resource(self) -> Dict[str, Any]:
        """Serializes this object to a dictionary for API requests."""
        return {_ALLOWED_IP_CIDR_RANGES: self.allowed_ip_cidr_ranges}


class VpcNetworkSource:
    """Represents a VPC network source for a GCS Bucket IP Filter.

    :type network: str
    :param network: The resource name of the VPC network.

    :type allowed_ip_cidr_ranges: list(str) or None
    :param allowed_ip_cidr_ranges: A list of IPv4 or IPv6 ranges in CIDR
                                   notation allowed to access the bucket
                                   from this VPC.
    """

    def __init__(
        self, network: str, allowed_ip_cidr_ranges: Optional[List[str]] = None
    ):
        self.network = network
        self.allowed_ip_cidr_ranges = allowed_ip_cidr_ranges or []

    def _to_api_resource(self) -> Dict[str, Any]:
        """Serializes this object to a dictionary for API requests."""
        return {
            _NETWORK: self.network,
            _ALLOWED_IP_CIDR_RANGES: self.allowed_ip_cidr_ranges,
        }


class IPFilter:
    """Represents a GCS Bucket IP Filter configuration.

    This class is a helper for constructing the IP Filter dictionary to be
    assigned to a bucket's ``ip_filter`` property.
    """

    """
    Attributes:
        mode (str): Required. The mode of the IP filter. Can be "Enabled" or "Disabled".
        allow_all_service_agent_access (bool): Required. If True, allows Google
            Cloud service agents to bypass the IP filter.
        public_network_source (PublicNetworkSource): (Optional) The configuration
            for requests from the public internet.
        vpc_network_sources (list(VpcNetworkSource)): (Optional) A list of
            configurations for requests from VPC networks.
        allow_cross_org_vpcs (bool): (Optional) If True, allows VPCs from
            other organizations to be used in the configuration.
    """

    def __init__(self):
        self.mode: Optional[str] = None
        self.public_network_source: Optional[PublicNetworkSource] = None
        self.vpc_network_sources: List[VpcNetworkSource] = []
        self.allow_all_service_agent_access: Optional[bool] = None
        self.allow_cross_org_vpcs: Optional[bool] = None

    @classmethod
    def _from_api_resource(cls, resource: Dict[str, Any]) -> "IPFilter":
        """Factory: creates an IPFilter instance from a server response."""
        ip_filter = cls()
        ip_filter.mode = resource.get(_MODE)
        ip_filter.allow_all_service_agent_access = resource.get(
            _ALLOW_ALL_SERVICE_AGENT_ACCESS, None
        )

        public_network_source_data = resource.get(_PUBLIC_NETWORK_SOURCE, None)
        if public_network_source_data:
            ip_filter.public_network_source = PublicNetworkSource(
                allowed_ip_cidr_ranges=public_network_source_data.get(
                    _ALLOWED_IP_CIDR_RANGES, []
                )
            )

        vns_res_list = resource.get(_VPC_NETWORK_SOURCES, [])
        ip_filter.vpc_network_sources = [
            VpcNetworkSource(
                network=vns.get(_NETWORK),
                allowed_ip_cidr_ranges=vns.get(_ALLOWED_IP_CIDR_RANGES, []),
            )
            for vns in vns_res_list
        ]
        ip_filter.allow_cross_org_vpcs = resource.get(_ALLOW_CROSS_ORG_VPCS, None)
        return ip_filter

    def _to_api_resource(self) -> Dict[str, Any]:
        """Serializes this object to a dictionary for API requests."""
        resource = {
            _MODE: self.mode,
            _ALLOW_ALL_SERVICE_AGENT_ACCESS: self.allow_all_service_agent_access,
        }

        if self.public_network_source:
            resource[
                _PUBLIC_NETWORK_SOURCE
            ] = self.public_network_source._to_api_resource()
        if self.vpc_network_sources is not None:
            resource[_VPC_NETWORK_SOURCES] = [
                vns._to_api_resource() for vns in self.vpc_network_sources
            ]
        if self.allow_cross_org_vpcs is not None:
            resource[_ALLOW_CROSS_ORG_VPCS] = self.allow_cross_org_vpcs
        return resource
