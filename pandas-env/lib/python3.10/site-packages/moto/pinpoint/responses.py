"""Handles incoming pinpoint requests, invokes methods, returns responses."""

import json
from urllib.parse import unquote

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse

from .models import PinpointBackend, pinpoint_backends


class PinpointResponse(BaseResponse):
    """Handler for Pinpoint requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="pinpoint")

    @property
    def pinpoint_backend(self) -> PinpointBackend:
        """Return backend instance specific for this region."""
        return pinpoint_backends[self.current_account][self.region]

    def create_app(self) -> TYPE_RESPONSE:
        params = json.loads(self.body)
        name = params.get("Name")
        tags = params.get("tags", {})
        app = self.pinpoint_backend.create_app(name=name, tags=tags)
        return 201, {"status": 201}, json.dumps(app.to_json())

    def delete_app(self) -> str:
        application_id = self.path.split("/")[-1]
        app = self.pinpoint_backend.delete_app(application_id=application_id)
        return json.dumps(app.to_json())

    def get_app(self) -> str:
        application_id = self.path.split("/")[-1]
        app = self.pinpoint_backend.get_app(application_id=application_id)
        return json.dumps(app.to_json())

    def get_apps(self) -> str:
        apps = self.pinpoint_backend.get_apps()
        resp = {"Item": [a.to_json() for a in apps]}
        return json.dumps(resp)

    def update_application_settings(self) -> str:
        application_id = self.path.split("/")[-2]
        settings = json.loads(self.body)
        app_settings = self.pinpoint_backend.update_application_settings(
            application_id=application_id, settings=settings
        )
        response = app_settings.to_json()
        response["ApplicationId"] = application_id
        return json.dumps(response)

    def get_application_settings(self) -> str:
        application_id = self.path.split("/")[-2]
        app_settings = self.pinpoint_backend.get_application_settings(
            application_id=application_id
        )
        response = app_settings.to_json()
        response["ApplicationId"] = application_id
        return json.dumps(response)

    def list_tags_for_resource(self) -> str:
        resource_arn = unquote(self.path).split("/tags/")[-1]
        tags = self.pinpoint_backend.list_tags_for_resource(resource_arn=resource_arn)
        return json.dumps(tags)

    def tag_resource(self) -> str:
        resource_arn = unquote(self.path).split("/tags/")[-1]
        tags = json.loads(self.body).get("tags", {})
        self.pinpoint_backend.tag_resource(resource_arn=resource_arn, tags=tags)
        return "{}"

    def untag_resource(self) -> str:
        resource_arn = unquote(self.path).split("/tags/")[-1]
        tag_keys = self.querystring.get("tagKeys")
        self.pinpoint_backend.untag_resource(
            resource_arn=resource_arn,
            tag_keys=tag_keys,  # type: ignore[arg-type]
        )
        return "{}"

    def put_event_stream(self) -> str:
        application_id = self.path.split("/")[-2]
        params = json.loads(self.body)
        stream_arn = params.get("DestinationStreamArn")
        role_arn = params.get("RoleArn")
        event_stream = self.pinpoint_backend.put_event_stream(
            application_id=application_id, stream_arn=stream_arn, role_arn=role_arn
        )
        resp = event_stream.to_json()
        resp["ApplicationId"] = application_id
        return json.dumps(resp)

    def get_event_stream(self) -> str:
        application_id = self.path.split("/")[-2]
        event_stream = self.pinpoint_backend.get_event_stream(
            application_id=application_id
        )
        resp = event_stream.to_json()
        resp["ApplicationId"] = application_id
        return json.dumps(resp)

    def delete_event_stream(self) -> str:
        application_id = self.path.split("/")[-2]
        event_stream = self.pinpoint_backend.delete_event_stream(
            application_id=application_id
        )
        resp = event_stream.to_json()
        resp["ApplicationId"] = application_id
        return json.dumps(resp)
