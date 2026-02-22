"""mediapackagev2Backend class with methods for supported APIs."""

from datetime import datetime

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.utilities.paginator import paginate

from .exceptions import ChannelGroupNotEmpty, ChannelGroupNotFound, ChannelNotFound
from .utils import PAGINATION_MODEL


class Channel(BaseModel):
    def __init__(self, group: "ChannelGroup", channel_name: str):
        self.group = group
        self.channel_group_name = group.channel_group_name
        self.channel_name = channel_name
        self.arn = f"{group.arn}/channel/{channel_name}"

        ingest_domain1 = group.egress_domain.replace(".egress", "-1.ingest")
        ingest_domain2 = group.egress_domain.replace(".egress", "-2.ingest")
        self.ingest_endpoints = [
            {
                "Id": "1",
                "Url": f"https://{ingest_domain1}/in/v1/{group.channel_group_name}/1/{channel_name}/index",
            },
            {
                "Id": "2",
                "Url": f"https://{ingest_domain2}/in/v1/{group.channel_group_name}/2/{channel_name}/index",
            },
        ]

        self.e_tag = mock_random.get_random_hex(10)
        self.tags: dict[str, str] = {}

        self.input_switch_configuration = {"MQCSInputSwitching": False}
        self.output_header_configuration = {"PublishMQCS": False}

        self.input_type = "HLS"


class ChannelGroup(BaseModel):
    def __init__(self, account_id: str, region_name: str, channel_group_name: str):
        self.channel_group_name = channel_group_name
        self.arn = f"arn:aws:mediapackagev2:{region_name}:{account_id}:channelGroup/{channel_group_name}"
        self.egress_domain = f"{mock_random.get_random_hex(6)}.egress.{mock_random.get_random_hex(6)}.mediapackagev2.{region_name}.amazonaws.com"
        self.created_at = datetime.now()
        self.modified_at = datetime.now()
        self.e_tag = mock_random.get_random_hex(10)
        self.tags: dict[str, str] = {}

        self.channels: dict[str, Channel] = {}

    def create_channel(self, channel_name: str) -> Channel:
        channel = Channel(group=self, channel_name=channel_name)
        self.channels[channel_name] = channel
        return channel


class MediaPackagev2Backend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.channel_groups: dict[str, ChannelGroup] = {}

    def create_channel_group(self, channel_group_name: str) -> ChannelGroup:
        group = ChannelGroup(
            account_id=self.account_id,
            region_name=self.region_name,
            channel_group_name=channel_group_name,
        )
        self.channel_groups[channel_group_name] = group
        return group

    def create_channel(self, channel_group_name: str, channel_name: str) -> Channel:
        group = self.get_channel_group(channel_group_name)
        return group.create_channel(channel_name)

    def get_channel_group(self, channel_group_name: str) -> ChannelGroup:
        if channel_group_name not in self.channel_groups:
            raise ChannelGroupNotFound
        return self.channel_groups[channel_group_name]

    @paginate(PAGINATION_MODEL)
    def list_channel_groups(self) -> list[ChannelGroup]:
        return list(self.channel_groups.values())

    def get_channel(self, channel_group_name: str, channel_name: str) -> Channel:
        group = self.get_channel_group(channel_group_name)
        if channel_name not in group.channels:
            raise ChannelNotFound
        return group.channels[channel_name]

    def delete_channel_group(self, channel_group_name: str) -> None:
        if channel := self.channel_groups.get(channel_group_name):
            if channel.channels:
                raise ChannelGroupNotEmpty
        self.channel_groups.pop(channel_group_name, None)

    def delete_channel(self, channel_group_name: str, channel_name: str) -> None:
        group = self.get_channel_group(channel_group_name)
        group.channels.pop(channel_name, None)


mediapackagev2_backends = BackendDict(MediaPackagev2Backend, "mediapackagev2")
