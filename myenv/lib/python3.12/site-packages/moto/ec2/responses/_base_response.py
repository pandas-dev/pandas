from typing import Any, Dict, Optional

from moto.core.responses import BaseResponse

from ..exceptions import EmptyTagSpecError, InvalidParameter
from ..utils import convert_tag_spec


class EC2BaseResponse(BaseResponse):
    @property
    def ec2_backend(self) -> Any:  # type: ignore[misc]
        from moto.ec2.models import ec2_backends

        return ec2_backends[self.current_account][self.region]

    def _filters_from_querystring(self) -> Dict[str, str]:
        # [{"Name": x1, "Value": y1}, ..]
        _filters = self._get_multi_param("Filter.")
        # return {x1: y1, ...}
        return {f["Name"]: f["Value"] for f in _filters}

    def _parse_tag_specification(
        self, expected_type: Optional[str] = None
    ) -> Dict[str, Dict[str, str]]:
        # [{"ResourceType": _type, "Tag": [{"Key": k, "Value": v}, ..]}]
        tag_spec_set = self._get_multi_param(
            "TagSpecification", skip_result_conversion=True
        )
        if not tag_spec_set:
            tag_spec_set = self._get_multi_param(
                "TagSpecifications", skip_result_conversion=True
            )
        if not tag_spec_set:
            return {}

        tags_dict = (
            tag_spec_set[0] if isinstance(tag_spec_set, list) else tag_spec_set
        )  # awscli allows for a json string to be passed and it should be allowed
        if "ResourceType" not in tags_dict:
            raise InvalidParameter("Tag specification resource type must have a value")
        if expected_type and tags_dict["ResourceType"] != expected_type:
            raise InvalidParameter(
                f"'{tags_dict['ResourceType']}' is not a valid taggable resource type for this operation."
            )
        if "Tag" not in tags_dict:
            if tags_dict.get("ResourceType") == "subnet":
                raise InvalidParameter("Tag specification must have at least one tag")
            raise EmptyTagSpecError

        return convert_tag_spec(tag_spec_set)
