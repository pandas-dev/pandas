"""Handles incoming signer requests, invokes methods, returns responses."""

import json
from typing import Any
from urllib.parse import unquote

from moto.core.responses import BaseResponse

from .models import SignerBackend, signer_backends


class SignerResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="signer")

    @property
    def signer_backend(self) -> SignerBackend:
        """Return backend instance specific for this region."""
        return signer_backends[self.current_account][self.region]

    def cancel_signing_profile(self) -> str:
        profile_name = self.path.split("/")[-1]
        self.signer_backend.cancel_signing_profile(profile_name=profile_name)
        return "{}"

    def get_signing_profile(self) -> str:
        profile_name = self.path.split("/")[-1]
        profile = self.signer_backend.get_signing_profile(profile_name=profile_name)
        return json.dumps(profile.to_dict())

    def put_signing_profile(self) -> str:
        params = json.loads(self.body)
        profile_name = self.path.split("/")[-1]
        signature_validity_period = params.get("signatureValidityPeriod")
        platform_id = params.get("platformId")
        tags = params.get("tags")
        signing_material = params.get("signingMaterial")
        profile = self.signer_backend.put_signing_profile(
            profile_name=profile_name,
            signature_validity_period=signature_validity_period,
            platform_id=platform_id,
            signing_material=signing_material,
            tags=tags,
        )
        return json.dumps(profile.to_dict(full=False))

    def list_signing_platforms(self) -> str:
        platforms = self.signer_backend.list_signing_platforms()
        return json.dumps(dict(platforms=platforms))

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        return json.dumps(
            {"tags": self.signer_backend.list_tags_for_resource(resource_arn)}
        )

    def tag_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tags = self._get_param("tags")
        self.signer_backend.tag_resource(resource_arn, tags)
        return "{}"

    def untag_resource(self) -> str:
        resource_arn = unquote(self.path.split("/tags/")[-1])
        tag_keys = self.querystring.get("tagKeys")
        self.signer_backend.untag_resource(resource_arn, tag_keys)  # type: ignore
        return "{}"

    def tags(self, request: Any, full_url: str, headers: Any) -> str:  # type: ignore[return]
        self.setup_class(request, full_url, headers)
        if request.method == "GET":
            return self.list_tags_for_resource()
        if request.method == "POST":
            return self.tag_resource()
        if request.method == "DELETE":
            return self.untag_resource()
