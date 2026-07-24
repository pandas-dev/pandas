from moto.core.responses import ActionResult, BaseResponse, EmptyResult

from .models import CloudFrontBackend, cloudfront_backends


class CloudFrontResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="cloudfront")
        self.automated_parameter_parsing = True

    @property
    def backend(self) -> CloudFrontBackend:
        return cloudfront_backends[self.current_account][self.partition]

    def _get_action(self) -> str:
        # This is needed because the uri matcher doesn't take queryargs into account
        action = super()._get_action()
        if action == "CreateDistribution" and "WithTags" in self.querystring:
            action = "CreateDistributionWithTags"
        elif action is None and "Operation" in self.querystring:
            op_to_action = {"Tag": "TagResource", "Untag": "UntagResource"}
            operation = self.querystring.get("Operation")[0]
            action = op_to_action.get(operation, action)
        return action

    def create_distribution(self) -> ActionResult:
        distribution_config = self._get_param("DistributionConfig", {})
        distribution, location, e_tag = self.backend.create_distribution(
            distribution_config=distribution_config,
            tags=[],
        )
        result = {"Distribution": distribution, "ETag": e_tag, "Location": location}
        return ActionResult(result)

    def create_distribution_with_tags(self) -> ActionResult:
        distribution_config = self._get_param(
            "DistributionConfigWithTags.DistributionConfig", {}
        )
        tags = self._get_param("DistributionConfigWithTags.Tags.Items", [])
        distribution, location, e_tag = self.backend.create_distribution(
            distribution_config=distribution_config,
            tags=tags,
        )
        result = {"Distribution": distribution, "ETag": e_tag, "Location": location}
        return ActionResult(result)

    def list_distributions(self) -> ActionResult:
        distributions = self.backend.list_distributions()
        result = {
            "DistributionList": {
                "Marker": "",
                "MaxItems": 100,
                "IsTruncated": False,
                "Quantity": len(distributions),
                "Items": distributions if distributions else None,
            }
        }
        return ActionResult(result)

    def delete_distribution(self) -> ActionResult:
        distribution_id = self._get_param("Id")
        if_match = self._get_param("IfMatch")
        self.backend.delete_distribution(distribution_id, if_match)
        return EmptyResult()

    def get_distribution(self) -> ActionResult:
        distribution_id = self._get_param("Id")
        dist, etag = self.backend.get_distribution(distribution_id)
        result = {"Distribution": dist, "ETag": etag}
        return ActionResult(result)

    def get_distribution_config(self) -> ActionResult:
        dist_id = self._get_param("Id")
        distribution_config, etag = self.backend.get_distribution_config(dist_id)
        result = {"DistributionConfig": distribution_config, "ETag": etag}
        return ActionResult(result)

    def update_distribution(self) -> ActionResult:
        dist_id = self._get_param("Id")
        dist_config = self._get_param("DistributionConfig", {})
        if_match = self._get_param("IfMatch")
        dist, location, e_tag = self.backend.update_distribution(
            dist_config=dist_config,
            _id=dist_id,
            if_match=if_match,
        )
        result = {"Distribution": dist, "ETag": e_tag, "Location": location}
        return ActionResult(result)

    def create_invalidation(self) -> ActionResult:
        dist_id = self._get_param("DistributionId")
        paths = self._get_param("InvalidationBatch.Paths.Items", [])
        caller_ref = self._get_param("InvalidationBatch.CallerReference")
        invalidation = self.backend.create_invalidation(dist_id, paths, caller_ref)
        result = {"Invalidation": invalidation, "Location": invalidation.location}
        return ActionResult(result)

    def list_invalidations(self) -> ActionResult:
        dist_id = self._get_param("DistributionId")
        invalidations = self.backend.list_invalidations(dist_id)
        result = {
            "InvalidationList": {
                "MaxItems": 100,
                "IsTruncated": False,
                "Quantity": len(invalidations),
                "Items": invalidations if invalidations else None,
            }
        }
        return ActionResult(result)

    def get_invalidation(self) -> ActionResult:
        invalidation_id = self._get_param("Id")
        dist_id = self._get_param("DistributionId")
        invalidation = self.backend.get_invalidation(dist_id, invalidation_id)
        result = {"Invalidation": invalidation}
        return ActionResult(result)

    def list_tags_for_resource(self) -> ActionResult:
        resource = self._get_param("Resource")
        tags = self.backend.list_tags_for_resource(resource=resource)["Tags"]
        result = {"Tags": {"Items": tags}}
        return ActionResult(result)

    def tag_resource(self) -> ActionResult:
        resource = self._get_param("Resource")
        tags = self._get_param("Tags.Items", []) or []
        tags = {tag["Key"]: tag.get("Value") for tag in tags}
        self.backend.tag_resource(resource, tags)
        return EmptyResult()

    def untag_resource(self) -> ActionResult:
        resource = self._get_param("Resource")
        tag_keys_data = self._get_param("TagKeys.Items", []) or []
        self.backend.untag_resource(resource, tag_keys_data)
        return EmptyResult()

    def create_origin_access_control(self) -> ActionResult:
        config = self._get_param("OriginAccessControlConfig", {})
        control = self.backend.create_origin_access_control(config)
        result = {
            "OriginAccessControl": {
                "Id": control.id,
                "OriginAccessControlConfig": control,
            },
            "ETag": control.etag,
        }
        return ActionResult(result)

    def get_origin_access_control(self) -> ActionResult:
        control_id = self._get_param("Id")
        control = self.backend.get_origin_access_control(control_id)
        result = {
            "OriginAccessControl": {
                "Id": control.id,
                "OriginAccessControlConfig": control,
            },
            "ETag": control.etag,
        }
        return ActionResult(result)

    def list_origin_access_controls(self) -> ActionResult:
        controls = self.backend.list_origin_access_controls()
        result = {
            "OriginAccessControlList": {
                "MaxItems": 100,
                "IsTruncated": False,
                "Quantity": len(controls),
                "Items": controls,
            }
        }
        return ActionResult(result)

    def update_origin_access_control(self) -> ActionResult:
        control_id = self._get_param("Id")
        config = self._get_param("OriginAccessControlConfig", {})
        control = self.backend.update_origin_access_control(control_id, config)
        result = {
            "OriginAccessControl": {
                "Id": control.id,
                "OriginAccessControlConfig": control,
            },
            "ETag": control.etag,
        }
        return ActionResult(result)

    def delete_origin_access_control(self) -> ActionResult:
        control_id = self._get_param("Id")
        self.backend.delete_origin_access_control(control_id)
        return EmptyResult()

    def create_public_key(self) -> ActionResult:
        config = self._get_param("PublicKeyConfig")
        caller_ref = config["CallerReference"]
        name = config["Name"]
        encoded_key = config["EncodedKey"]
        public_key = self.backend.create_public_key(
            caller_ref=caller_ref, name=name, encoded_key=encoded_key
        )
        result = {
            "PublicKey": public_key,
            "Location": public_key.location,
            "ETag": public_key.etag,
        }
        return ActionResult(result)

    def get_public_key(self) -> ActionResult:
        key_id = self._get_param("Id")
        public_key = self.backend.get_public_key(key_id=key_id)
        result = {"PublicKey": public_key, "ETag": public_key.etag}
        return ActionResult(result)

    def delete_public_key(self) -> ActionResult:
        key_id = self._get_param("Id")
        self.backend.delete_public_key(key_id=key_id)
        return EmptyResult()

    def list_public_keys(self) -> ActionResult:
        keys = self.backend.list_public_keys()
        result = {
            "PublicKeyList": {
                "MaxItems": 100,
                "Quantity": len(keys),
                "Items": keys if keys else None,
            }
        }
        return ActionResult(result)

    def create_key_group(self) -> ActionResult:
        name = self._get_param("KeyGroupConfig.Name")
        items = self._get_param("KeyGroupConfig.Items", [])
        key_group = self.backend.create_key_group(name=name, items=items)
        result = {
            "KeyGroup": key_group,
            "Location": key_group.location,
            "ETag": key_group.etag,
        }
        return ActionResult(result)

    def get_key_group(self) -> ActionResult:
        group_id = self._get_param("Id")
        key_group = self.backend.get_key_group(group_id=group_id)
        result = {"KeyGroup": key_group, "ETag": key_group.etag}
        return ActionResult(result)

    def list_key_groups(self) -> ActionResult:
        groups = self.backend.list_key_groups()
        result = {
            "KeyGroupList": {
                "Quantity": len(groups),
                "Items": [{"KeyGroup": key_group} for key_group in groups],
            }
        }
        return ActionResult(result)
