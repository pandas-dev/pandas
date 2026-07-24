from moto.core.responses import ActionResult, EmptyResult

from ._base_response import EC2BaseResponse


class InternetGateways(EC2BaseResponse):
    def attach_internet_gateway(self) -> ActionResult:
        igw_id = self._get_param("InternetGatewayId")
        vpc_id = self._get_param("VpcId")

        self.error_on_dryrun()

        self.ec2_backend.attach_internet_gateway(igw_id, vpc_id)
        return EmptyResult()

    def create_internet_gateway(self) -> ActionResult:
        self.error_on_dryrun()

        tags = self._get_param("TagSpecifications", [])
        if tags:
            tags = tags[0].get("Tags") or []
        igw = self.ec2_backend.create_internet_gateway(tags=tags)
        return ActionResult({"InternetGateway": igw})

    def delete_internet_gateway(self) -> ActionResult:
        igw_id = self._get_param("InternetGatewayId")
        self.error_on_dryrun()

        self.ec2_backend.delete_internet_gateway(igw_id)
        return EmptyResult()

    def describe_internet_gateways(self) -> ActionResult:
        filter_dict = self._filters_from_querystring()
        igw_ids = self._get_param("InternetGatewayIds", None)
        igws = self.ec2_backend.describe_internet_gateways(igw_ids, filters=filter_dict)
        return ActionResult({"InternetGateways": igws})

    def detach_internet_gateway(self) -> ActionResult:
        # TODO validate no instances with EIPs in VPC before detaching
        # raise else DependencyViolationError()
        igw_id = self._get_param("InternetGatewayId")
        vpc_id = self._get_param("VpcId")
        self.error_on_dryrun()

        self.ec2_backend.detach_internet_gateway(igw_id, vpc_id)
        return EmptyResult()
