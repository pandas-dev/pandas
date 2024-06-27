from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import unix_time
from moto.moto_api._internal import mock_random
from moto.utilities.tagging_service import TaggingService
from moto.utilities.utils import get_partition

from .exceptions import ApplicationNotFound, EventStreamNotFound


class App(BaseModel):
    def __init__(self, account_id: str, region_name: str, name: str):
        self.application_id = str(mock_random.uuid4()).replace("-", "")
        self.arn = f"arn:{get_partition(region_name)}:mobiletargeting:us-east-1:{account_id}:apps/{self.application_id}"
        self.name = name
        self.created = unix_time()
        self.settings = AppSettings()
        self.event_stream: Optional[EventStream] = None

    def get_settings(self) -> "AppSettings":
        return self.settings

    def update_settings(self, settings: Dict[str, Any]) -> "AppSettings":
        self.settings.update(settings)
        return self.settings

    def delete_event_stream(self) -> "EventStream":
        stream = self.event_stream
        self.event_stream = None
        return stream  # type: ignore

    def get_event_stream(self) -> "EventStream":
        if self.event_stream is None:
            raise EventStreamNotFound()
        return self.event_stream

    def put_event_stream(self, stream_arn: str, role_arn: str) -> "EventStream":
        self.event_stream = EventStream(stream_arn, role_arn)
        return self.event_stream

    def to_json(self) -> Dict[str, Any]:
        return {
            "Arn": self.arn,
            "Id": self.application_id,
            "Name": self.name,
            "CreationDate": self.created,
        }


class AppSettings(BaseModel):
    def __init__(self) -> None:
        self.settings: Dict[str, Any] = dict()
        self.last_modified = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%fZ")

    def update(self, settings: Dict[str, Any]) -> None:
        self.settings = settings
        self.last_modified = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%fZ")

    def to_json(self) -> Dict[str, Any]:
        return {
            "CampaignHook": self.settings.get("CampaignHook", {}),
            "CloudWatchMetricsEnabled": self.settings.get(
                "CloudWatchMetricsEnabled", False
            ),
            "LastModifiedDate": self.last_modified,
            "Limits": self.settings.get("Limits", {}),
            "QuietTime": self.settings.get("QuietTime", {}),
        }


class EventStream(BaseModel):
    def __init__(self, stream_arn: str, role_arn: str):
        self.stream_arn = stream_arn
        self.role_arn = role_arn
        self.last_modified = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%fZ")

    def to_json(self) -> Dict[str, Any]:
        return {
            "DestinationStreamArn": self.stream_arn,
            "RoleArn": self.role_arn,
            "LastModifiedDate": self.last_modified,
        }


class PinpointBackend(BaseBackend):
    """Implementation of Pinpoint APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.apps: Dict[str, App] = {}
        self.tagger = TaggingService()

    def create_app(self, name: str, tags: Dict[str, str]) -> App:
        app = App(self.account_id, self.region_name, name)
        self.apps[app.application_id] = app
        tag_list = self.tagger.convert_dict_to_tags_input(tags)
        self.tagger.tag_resource(app.arn, tag_list)
        return app

    def delete_app(self, application_id: str) -> App:
        self.get_app(application_id)
        return self.apps.pop(application_id)

    def get_app(self, application_id: str) -> App:
        if application_id not in self.apps:
            raise ApplicationNotFound()
        return self.apps[application_id]

    def get_apps(self) -> Iterable[App]:
        """
        Pagination is not yet implemented
        """
        return self.apps.values()

    def update_application_settings(
        self, application_id: str, settings: Dict[str, Any]
    ) -> AppSettings:
        app = self.get_app(application_id)
        return app.update_settings(settings)

    def get_application_settings(self, application_id: str) -> AppSettings:
        app = self.get_app(application_id)
        return app.get_settings()

    def list_tags_for_resource(self, resource_arn: str) -> Dict[str, Dict[str, str]]:
        tags = self.tagger.get_tag_dict_for_resource(resource_arn)
        return {"tags": tags}

    def tag_resource(self, resource_arn: str, tags: Dict[str, str]) -> None:
        tag_list = TaggingService.convert_dict_to_tags_input(tags)
        self.tagger.tag_resource(resource_arn, tag_list)

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def put_event_stream(
        self, application_id: str, stream_arn: str, role_arn: str
    ) -> EventStream:
        app = self.get_app(application_id)
        return app.put_event_stream(stream_arn, role_arn)

    def get_event_stream(self, application_id: str) -> EventStream:
        app = self.get_app(application_id)
        return app.get_event_stream()

    def delete_event_stream(self, application_id: str) -> EventStream:
        app = self.get_app(application_id)
        return app.delete_event_stream()


pinpoint_backends = BackendDict(PinpointBackend, "pinpoint")
