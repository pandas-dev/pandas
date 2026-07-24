"""Exceptions raised by the mediapackagev2 service."""

from moto.core.exceptions import JsonRESTError


class ChannelGroupNotEmpty(JsonRESTError):
    def __init__(self) -> None:
        msg = "The channel group you tried to delete has channels attached to it. If you want to delete this channel group, you must first delete the channels attached to it"
        super().__init__(error_type="ConflictException", message=msg)


class ChannelNotFound(JsonRESTError):
    def __init__(self) -> None:
        msg = "MediaPackage can't process your request because we can't find your channel. Verify your channel name or add a channel and then try again"
        super().__init__(error_type="ResourceNotFoundException", message=msg)


class ChannelGroupNotFound(JsonRESTError):
    def __init__(self) -> None:
        msg = "MediaPackage can't process your request because we can't find your channel group. Verify your channel group name or add a channel group and then try again"
        super().__init__(error_type="ResourceNotFoundException", message=msg)
