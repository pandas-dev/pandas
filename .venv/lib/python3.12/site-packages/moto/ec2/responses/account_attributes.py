from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class AccountAttributes(EC2BaseResponse):
    def describe_account_attributes(self) -> ActionResult:
        attribute_defaults: dict[str, list[str]] = {
            "vpc-max-security-groups-per-interface": ["5"],
            "max-instances": ["20"],
            "supported-platforms": ["EC2", "VPC"],
            "default-vpc": ["none"],
            "max-elastic-ips": ["5"],
            "vpc-max-elastic-ips": ["5"],
        }
        result = {
            "AccountAttributes": [
                {
                    "AttributeName": k,
                    "AttributeValues": [{"AttributeValue": v} for v in values],
                }
                for k, values in attribute_defaults.items()
            ]
        }
        return ActionResult(result)
