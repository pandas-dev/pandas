from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class AvailabilityZonesAndRegions(EC2BaseResponse):
    def describe_availability_zones(self) -> ActionResult:
        self.error_on_dryrun()
        filters = self._filters_from_querystring()
        zone_names = self._get_multi_param("ZoneName")
        zone_ids = self._get_multi_param("ZoneId")
        zones = self.ec2_backend.describe_availability_zones(
            filters, zone_names=zone_names, zone_ids=zone_ids
        )
        result = {"AvailabilityZones": zones}
        return ActionResult(result)

    def describe_regions(self) -> ActionResult:
        self.error_on_dryrun()
        region_names = self._get_multi_param("RegionName")
        regions = self.ec2_backend.describe_regions(region_names)
        result = {"Regions": regions}
        return ActionResult(result)
