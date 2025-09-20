"""Handles incoming shield requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import ShieldBackend, shield_backends


class ShieldResponse(BaseResponse):
    """Handler for Shield requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="shield")

    @property
    def shield_backend(self) -> ShieldBackend:
        """Return backend instance specific for this region."""
        return shield_backends[self.current_account][self.region]

    def create_protection(self) -> str:
        params = json.loads(self.body)
        name = params.get("Name")
        resource_arn = params.get("ResourceArn")
        tags = params.get("Tags")
        protection_id = self.shield_backend.create_protection(
            name=name,
            resource_arn=resource_arn,
            tags=tags,
        )
        return json.dumps(dict(ProtectionId=protection_id))

    def describe_protection(self) -> str:
        params = json.loads(self.body)
        protection_id = params.get("ProtectionId")
        resource_arn = params.get("ResourceArn")
        protection = self.shield_backend.describe_protection(
            protection_id=protection_id,
            resource_arn=resource_arn,
        )
        return json.dumps(dict(Protection=protection.to_dict()))

    def list_protections(self) -> str:
        params = json.loads(self.body)
        inclusion_filters = params.get("InclusionFilters")
        protections = self.shield_backend.list_protections(
            inclusion_filters=inclusion_filters,
        )
        return json.dumps(
            dict(Protections=list(protection.to_dict() for protection in protections))
        )

    def delete_protection(self) -> str:
        params = json.loads(self.body)
        protection_id = params.get("ProtectionId")
        self.shield_backend.delete_protection(
            protection_id=protection_id,
        )
        return "{}"

    def list_tags_for_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tags = self.shield_backend.list_tags_for_resource(
            resource_arn=resource_arn,
        )
        return json.dumps(dict(Tags=tags))

    def tag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tags = params.get("Tags")
        self.shield_backend.tag_resource(
            resource_arn=resource_arn,
            tags=tags,
        )
        return "{}"

    def untag_resource(self) -> str:
        params = json.loads(self.body)
        resource_arn = params.get("ResourceARN")
        tag_keys = params.get("TagKeys")
        self.shield_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,
        )
        return "{}"

    def create_subscription(self) -> str:
        self.shield_backend.create_subscription()
        return "{}"

    def describe_subscription(self) -> str:
        subscription = self.shield_backend.describe_subscription()
        return json.dumps(dict(Subscription=subscription.to_dict()))
