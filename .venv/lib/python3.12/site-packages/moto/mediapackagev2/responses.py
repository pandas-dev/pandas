from moto.core.responses import ActionResult, BaseResponse, EmptyResult
from moto.core.serialize import never_return

from .models import MediaPackagev2Backend, mediapackagev2_backends


class mediapackagev2Response(BaseResponse):
    """Handler for mediapackagev2 requests and responses."""

    RESPONSE_KEY_PATH_TO_TRANSFORMER = {
        # InputType is only returned on GET
        "CreateChannel.CreateChannelResponse.InputType": never_return
    }

    def __init__(self) -> None:
        super().__init__(service_name="mediapackagev2")

    @property
    def mediapackagev2_backend(self) -> MediaPackagev2Backend:
        return mediapackagev2_backends[self.current_account][self.region]

    def create_channel_group(self) -> ActionResult:
        channel_group_name = self._get_param("ChannelGroupName")
        group = self.mediapackagev2_backend.create_channel_group(
            channel_group_name=channel_group_name,
        )
        return ActionResult(result=group)

    def create_channel(self) -> ActionResult:
        channel_group_name = self._get_param("ChannelGroupName")
        channel_name = self._get_param("ChannelName")
        channel = self.mediapackagev2_backend.create_channel(
            channel_group_name=channel_group_name,
            channel_name=channel_name,
        )
        return ActionResult(result=channel)

    def get_channel_group(self) -> ActionResult:
        channel_group_name = self._get_param("ChannelGroupName")
        group = self.mediapackagev2_backend.get_channel_group(
            channel_group_name=channel_group_name,
        )
        return ActionResult(result=group)

    def list_channel_groups(self) -> ActionResult:
        max_results = self._get_int_param("maxResults")
        next_token = self._get_param("nextToken")

        groups, next_token = self.mediapackagev2_backend.list_channel_groups(
            max_results=max_results,
            next_token=next_token,
        )
        return ActionResult(result={"Items": groups, "NextToken": next_token})

    def get_channel(self) -> ActionResult:
        channel_group_name = self._get_param("ChannelGroupName")
        channel_name = self._get_param("ChannelName")
        channel = self.mediapackagev2_backend.get_channel(
            channel_group_name=channel_group_name,
            channel_name=channel_name,
        )
        return ActionResult(result=channel)

    def delete_channel_group(self) -> ActionResult:
        channel_group_name = self._get_param("ChannelGroupName")
        self.mediapackagev2_backend.delete_channel_group(
            channel_group_name=channel_group_name
        )
        return EmptyResult()

    def delete_channel(self) -> ActionResult:
        channel_group_name = self._get_param("ChannelGroupName")
        channel_name = self._get_param("ChannelName")
        self.mediapackagev2_backend.delete_channel(
            channel_group_name=channel_group_name,
            channel_name=channel_name,
        )
        return EmptyResult()
