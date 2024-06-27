from typing import Any, Dict, List, Optional

from boto3 import Session

from moto.utilities.utils import filter_resources


class Region:
    def __init__(self, name: str, endpoint: str, opt_in_status: str):
        self.name = name
        self.endpoint = endpoint
        self.opt_in_status = opt_in_status


class Zone:
    def __init__(
        self,
        name: str,
        region_name: str,
        zone_id: str,
        zone_type: str = "availability-zone",
    ):
        self.name = name
        self.region_name = region_name
        self.zone_id = zone_id
        self.zone_type = zone_type
        self.state = "available"


class RegionsAndZonesBackend:
    regions_opt_in_not_required = [
        "af-south-1",
        "ap-northeast-1",
        "ap-northeast-2",
        "ap-northeast-3",
        "ap-south-1",
        "ap-southeast-1",
        "ap-southeast-2",
        "ap-southeast-3",
        "ca-central-1",
        "eu-central-1",
        "eu-north-1",
        "eu-south-1",
        "eu-west-1",
        "eu-west-2",
        "eu-west-3",
        "sa-east-1",
        "us-east-1",
        "us-east-2",
        "us-west-1",
        "us-west-2",
    ]

    regions = []
    for region in Session().get_available_regions("ec2"):
        if region in regions_opt_in_not_required:
            regions.append(
                Region(region, f"ec2.{region}.amazonaws.com", "opt-in-not-required")
            )
        else:
            regions.append(
                Region(region, f"ec2.{region}.amazonaws.com", "not-opted-in")
            )
    for region in Session().get_available_regions("ec2", partition_name="aws-us-gov"):
        regions.append(
            Region(region, f"ec2.{region}.amazonaws.com", "opt-in-not-required")
        )
    for region in Session().get_available_regions("ec2", partition_name="aws-cn"):
        regions.append(
            Region(region, f"ec2.{region}.amazonaws.com.cn", "opt-in-not-required")
        )

    zones = {
        "af-south-1": [
            Zone(region_name="af-south-1", name="af-south-1a", zone_id="afs1-az1"),
            Zone(region_name="af-south-1", name="af-south-1b", zone_id="afs1-az2"),
            Zone(region_name="af-south-1", name="af-south-1c", zone_id="afs1-az3"),
        ],
        "ap-south-1": [
            Zone(region_name="ap-south-1", name="ap-south-1a", zone_id="aps1-az1"),
            Zone(region_name="ap-south-1", name="ap-south-1b", zone_id="aps1-az3"),
            Zone(region_name="ap-south-1", name="ap-south-1c", zone_id="aps1-az2"),
        ],
        "eu-west-3": [
            Zone(region_name="eu-west-3", name="eu-west-3a", zone_id="euw3-az1"),
            Zone(region_name="eu-west-3", name="eu-west-3b", zone_id="euw3-az2"),
            Zone(region_name="eu-west-3", name="eu-west-3c", zone_id="euw3-az3"),
        ],
        "eu-north-1": [
            Zone(region_name="eu-north-1", name="eu-north-1a", zone_id="eun1-az1"),
            Zone(region_name="eu-north-1", name="eu-north-1b", zone_id="eun1-az2"),
            Zone(region_name="eu-north-1", name="eu-north-1c", zone_id="eun1-az3"),
        ],
        "eu-west-2": [
            Zone(region_name="eu-west-2", name="eu-west-2a", zone_id="euw2-az2"),
            Zone(region_name="eu-west-2", name="eu-west-2b", zone_id="euw2-az3"),
            Zone(region_name="eu-west-2", name="eu-west-2c", zone_id="euw2-az1"),
        ],
        "eu-west-1": [
            Zone(region_name="eu-west-1", name="eu-west-1a", zone_id="euw1-az3"),
            Zone(region_name="eu-west-1", name="eu-west-1b", zone_id="euw1-az1"),
            Zone(region_name="eu-west-1", name="eu-west-1c", zone_id="euw1-az2"),
        ],
        "ap-northeast-3": [
            Zone(
                region_name="ap-northeast-3",
                name="ap-northeast-3a",
                zone_id="apne3-az1",
            ),
            Zone(
                region_name="ap-northeast-3",
                name="ap-northeast-3b",
                zone_id="apne3-az2",
            ),
            Zone(
                region_name="ap-northeast-3",
                name="ap-northeast-3c",
                zone_id="apne3-az3",
            ),
        ],
        "ap-northeast-2": [
            Zone(
                region_name="ap-northeast-2",
                name="ap-northeast-2a",
                zone_id="apne2-az1",
            ),
            Zone(
                region_name="ap-northeast-2",
                name="ap-northeast-2b",
                zone_id="apne2-az2",
            ),
            Zone(
                region_name="ap-northeast-2",
                name="ap-northeast-2c",
                zone_id="apne2-az3",
            ),
            Zone(
                region_name="ap-northeast-2",
                name="ap-northeast-2d",
                zone_id="apne2-az4",
            ),
        ],
        "ap-northeast-1": [
            Zone(
                region_name="ap-northeast-1",
                name="ap-northeast-1a",
                zone_id="apne1-az4",
            ),
            Zone(
                region_name="ap-northeast-1",
                name="ap-northeast-1c",
                zone_id="apne1-az1",
            ),
            Zone(
                region_name="ap-northeast-1",
                name="ap-northeast-1d",
                zone_id="apne1-az2",
            ),
        ],
        "ap-east-1": [
            Zone(region_name="ap-east-1", name="ap-east-1a", zone_id="ape1-az1"),
            Zone(region_name="ap-east-1", name="ap-east-1b", zone_id="ape1-az2"),
            Zone(region_name="ap-east-1", name="ap-east-1c", zone_id="ape1-az3"),
        ],
        "sa-east-1": [
            Zone(region_name="sa-east-1", name="sa-east-1a", zone_id="sae1-az1"),
            Zone(region_name="sa-east-1", name="sa-east-1b", zone_id="sae1-az2"),
            Zone(region_name="sa-east-1", name="sa-east-1c", zone_id="sae1-az3"),
        ],
        "ca-central-1": [
            Zone(region_name="ca-central-1", name="ca-central-1a", zone_id="cac1-az1"),
            Zone(region_name="ca-central-1", name="ca-central-1b", zone_id="cac1-az2"),
            Zone(region_name="ca-central-1", name="ca-central-1d", zone_id="cac1-az4"),
        ],
        "ap-southeast-1": [
            Zone(
                region_name="ap-southeast-1",
                name="ap-southeast-1a",
                zone_id="apse1-az1",
            ),
            Zone(
                region_name="ap-southeast-1",
                name="ap-southeast-1b",
                zone_id="apse1-az2",
            ),
            Zone(
                region_name="ap-southeast-1",
                name="ap-southeast-1c",
                zone_id="apse1-az3",
            ),
        ],
        "ap-southeast-2": [
            Zone(
                region_name="ap-southeast-2",
                name="ap-southeast-2a",
                zone_id="apse2-az1",
            ),
            Zone(
                region_name="ap-southeast-2",
                name="ap-southeast-2b",
                zone_id="apse2-az3",
            ),
            Zone(
                region_name="ap-southeast-2",
                name="ap-southeast-2c",
                zone_id="apse2-az2",
            ),
        ],
        "ap-southeast-3": [
            Zone(
                region_name="ap-southeast-3",
                name="ap-southeast-3a",
                zone_id="apse3-az1",
            ),
            Zone(
                region_name="ap-southeast-3",
                name="ap-southeast-3b",
                zone_id="apse3-az2",
            ),
            Zone(
                region_name="ap-southeast-3",
                name="ap-southeast-3c",
                zone_id="apse3-az3",
            ),
        ],
        "eu-central-1": [
            Zone(region_name="eu-central-1", name="eu-central-1a", zone_id="euc1-az2"),
            Zone(region_name="eu-central-1", name="eu-central-1b", zone_id="euc1-az3"),
            Zone(region_name="eu-central-1", name="eu-central-1c", zone_id="euc1-az1"),
        ],
        "eu-central-2": [
            Zone(region_name="eu-central-2", name="eu-central-2a", zone_id="euc2-az1"),
            Zone(region_name="eu-central-2", name="eu-central-2b", zone_id="euc2-az2"),
            Zone(region_name="eu-central-2", name="eu-central-2c", zone_id="euc2-az3"),
        ],
        "eu-south-1": [
            Zone(region_name="eu-south-1", name="eu-south-1a", zone_id="eus1-az1"),
            Zone(region_name="eu-south-1", name="eu-south-1b", zone_id="eus1-az2"),
            Zone(region_name="eu-south-1", name="eu-south-1c", zone_id="eus1-az3"),
        ],
        "us-east-1": [
            Zone(region_name="us-east-1", name="us-east-1a", zone_id="use1-az6"),
            Zone(region_name="us-east-1", name="us-east-1b", zone_id="use1-az1"),
            Zone(region_name="us-east-1", name="us-east-1c", zone_id="use1-az2"),
            Zone(region_name="us-east-1", name="us-east-1d", zone_id="use1-az4"),
            Zone(region_name="us-east-1", name="us-east-1e", zone_id="use1-az3"),
            Zone(region_name="us-east-1", name="us-east-1f", zone_id="use1-az5"),
        ],
        "us-east-2": [
            Zone(region_name="us-east-2", name="us-east-2a", zone_id="use2-az1"),
            Zone(region_name="us-east-2", name="us-east-2b", zone_id="use2-az2"),
            Zone(region_name="us-east-2", name="us-east-2c", zone_id="use2-az3"),
        ],
        "us-west-1": [
            Zone(region_name="us-west-1", name="us-west-1a", zone_id="usw1-az3"),
            Zone(region_name="us-west-1", name="us-west-1b", zone_id="usw1-az1"),
        ],
        "us-west-2": [
            Zone(region_name="us-west-2", name="us-west-2a", zone_id="usw2-az2"),
            Zone(region_name="us-west-2", name="us-west-2b", zone_id="usw2-az1"),
            Zone(region_name="us-west-2", name="us-west-2c", zone_id="usw2-az3"),
            Zone(region_name="us-west-2", name="us-west-2d", zone_id="usw2-az4"),
        ],
        "me-south-1": [
            Zone(region_name="me-south-1", name="me-south-1a", zone_id="mes1-az1"),
            Zone(region_name="me-south-1", name="me-south-1b", zone_id="mes1-az2"),
            Zone(region_name="me-south-1", name="me-south-1c", zone_id="mes1-az3"),
        ],
        "cn-north-1": [
            Zone(region_name="cn-north-1", name="cn-north-1a", zone_id="cnn1-az1"),
            Zone(region_name="cn-north-1", name="cn-north-1b", zone_id="cnn1-az2"),
        ],
        "cn-northwest-1": [
            Zone(
                region_name="cn-northwest-1",
                name="cn-northwest-1a",
                zone_id="cnnw1-az1",
            ),
            Zone(
                region_name="cn-northwest-1",
                name="cn-northwest-1b",
                zone_id="cnnw1-az2",
            ),
            Zone(
                region_name="cn-northwest-1",
                name="cn-northwest-1c",
                zone_id="cnnw1-az3",
            ),
        ],
        "us-gov-west-1": [
            Zone(
                region_name="us-gov-west-1", name="us-gov-west-1a", zone_id="usgw1-az1"
            ),
            Zone(
                region_name="us-gov-west-1", name="us-gov-west-1b", zone_id="usgw1-az2"
            ),
            Zone(
                region_name="us-gov-west-1", name="us-gov-west-1c", zone_id="usgw1-az3"
            ),
        ],
        "us-gov-east-1": [
            Zone(
                region_name="us-gov-east-1", name="us-gov-east-1a", zone_id="usge1-az1"
            ),
            Zone(
                region_name="us-gov-east-1", name="us-gov-east-1b", zone_id="usge1-az2"
            ),
            Zone(
                region_name="us-gov-east-1", name="us-gov-east-1c", zone_id="usge1-az3"
            ),
        ],
    }

    def describe_regions(
        self, region_names: Optional[List[str]] = None
    ) -> List[Region]:
        if not region_names:
            return self.regions
        ret = []
        for name in region_names:
            for region in self.regions:
                if region.name == name:
                    ret.append(region)
        return ret

    def describe_availability_zones(
        self,
        filters: Optional[List[Dict[str, Any]]] = None,
        zone_names: Optional[List[str]] = None,
        zone_ids: Optional[List[str]] = None,
    ) -> List[Zone]:
        """
        The following parameters are supported: ZoneIds, ZoneNames, Filters
        The following filters are supported: zone-id, zone-type, zone-name, region-name, state
        """
        # We might not have any zones for the current region, if it was introduced recently
        zones = self.zones.get(self.region_name, [])  # type: ignore[attr-defined]
        attr_pairs = (
            ("zone-id", "zone_id"),
            ("zone-type", "zone_type"),
            ("zone-name", "name"),
            ("region-name", "region_name"),
            ("state", "state"),
        )
        result = zones
        if filters:
            result = filter_resources(zones, filters, attr_pairs)
        result = [r for r in result if not zone_ids or r.zone_id in zone_ids]
        result = [r for r in result if not zone_names or r.name in zone_names]
        return result

    def get_zone_by_name(self, name: str) -> Optional[Zone]:
        for zone in self.describe_availability_zones():
            if zone.name == name:
                return zone
        return None
