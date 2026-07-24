from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class KeyPairs(EC2BaseResponse):
    def create_key_pair(self) -> ActionResult:
        key_name = self._get_param("KeyName")
        key_type = self._get_param("KeyType")
        tags = self._parse_tag_specification("key-pair").get("key-pair", {})
        self.error_on_dryrun()
        key_pair = self.ec2_backend.create_key_pair(key_name, key_type, tags)
        return ActionResult(key_pair)

    def delete_key_pair(self) -> ActionResult:
        key_name = self._get_param("KeyName")
        self.error_on_dryrun()
        key_pair = self.ec2_backend.delete_key_pair(key_name)
        result = {"Return": True, "KeyPairId": key_pair.id if key_pair else None}
        return ActionResult(result)

    def describe_key_pairs(self) -> ActionResult:
        names = self._get_param("KeyNames", [])
        filters = self._filters_from_querystring()
        key_pairs = self.ec2_backend.describe_key_pairs(names, filters)
        include_public_key = self._get_bool_param("IncludePublicKey", False)
        public_key_attribute_path = (
            "DescribeKeyPairsResult.KeyPairs.KeyPairInfo.PublicKey"
        )
        if include_public_key:
            self._include_in_response(public_key_attribute_path)
        else:
            self._exclude_from_response(public_key_attribute_path)
        result = {"KeyPairs": key_pairs}
        return ActionResult(result)

    def import_key_pair(self) -> ActionResult:
        key_name = self._get_param("KeyName")
        material = self._get_param("PublicKeyMaterial")
        tags = self._parse_tag_specification("key-pair").get("key-pair", {})
        self.error_on_dryrun()
        key_pair = self.ec2_backend.import_key_pair(key_name, material, tags)
        return ActionResult(key_pair)
