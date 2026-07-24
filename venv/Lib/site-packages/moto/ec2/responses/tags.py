from moto.core.responses import ActionResult, EmptyResult
from moto.ec2.models import validate_resource_ids

from ._base_response import EC2BaseResponse


class TagResponse(EC2BaseResponse):
    def create_tags(self) -> ActionResult:
        resource_ids = self._get_param("Resources", [])
        validate_resource_ids(resource_ids)
        self.ec2_backend.do_resources_exist(resource_ids)
        tags = {tag["Key"]: tag["Value"] for tag in self._get_param("Tags", [])}

        self.error_on_dryrun()

        self.ec2_backend.create_tags(resource_ids, tags)
        return EmptyResult()

    def delete_tags(self) -> ActionResult:
        resource_ids = self._get_param("Resources", [])
        validate_resource_ids(resource_ids)
        tags = {tag["Key"]: tag.get("Value") for tag in self._get_param("Tags", [])}

        self.error_on_dryrun()

        self.ec2_backend.delete_tags(resource_ids, tags)
        return EmptyResult()

    def describe_tags(self) -> ActionResult:
        filters = self._filters_from_querystring()

        self.error_on_dryrun()

        tags = self.ec2_backend.describe_tags(filters=filters)
        result = {"Tags": tags}
        return ActionResult(result)
