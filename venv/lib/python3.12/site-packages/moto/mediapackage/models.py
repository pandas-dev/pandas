from collections import OrderedDict
from typing import Any, Dict, List

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.utils import get_partition

from .exceptions import ClientError


class Channel(BaseModel):
    def __init__(self, **kwargs: Any):
        self.arn = kwargs.get("arn")
        self.channel_id = kwargs.get("channel_id")
        self.description = kwargs.get("description")
        self.tags = kwargs.get("tags")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "arn": self.arn,
            "id": self.channel_id,
            "description": self.description,
            "tags": self.tags,
        }


class OriginEndpoint(BaseModel):
    def __init__(self, **kwargs: Any):
        self.arn = kwargs.get("arn")
        self.authorization = kwargs.get("authorization")
        self.channel_id = kwargs.get("channel_id")
        self.cmaf_package = kwargs.get("cmaf_package")
        self.dash_package = kwargs.get("dash_package")
        self.description = kwargs.get("description")
        self.hls_package = kwargs.get("hls_package")
        self.id = kwargs.get("endpoint_id")
        self.manifest_name = kwargs.get("manifest_name")
        self.mss_package = kwargs.get("mss_package")
        self.origination = kwargs.get("origination")
        self.startover_window_seconds = kwargs.get("startover_window_seconds")
        self.tags = kwargs.get("tags")
        self.time_delay_seconds = kwargs.get("time_delay_seconds")
        self.url = kwargs.get("url")
        self.whitelist = kwargs.get("whitelist")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "arn": self.arn,
            "authorization": self.authorization,
            "channelId": self.channel_id,
            "cmafPackage": self.cmaf_package,
            "dashPackage": self.dash_package,
            "description": self.description,
            "hlsPackage": self.hls_package,
            "id": self.id,
            "manifestName": self.manifest_name,
            "mssPackage": self.mss_package,
            "origination": self.origination,
            "startoverWindowSeconds": self.startover_window_seconds,
            "tags": self.tags,
            "timeDelaySeconds": self.time_delay_seconds,
            "url": self.url,
            "whitelist": self.whitelist,
        }


class MediaPackageBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._channels: Dict[str, Channel] = OrderedDict()
        self._origin_endpoints: Dict[str, OriginEndpoint] = OrderedDict()

    def create_channel(
        self, description: str, channel_id: str, tags: Dict[str, str]
    ) -> Channel:
        arn = f"arn:{get_partition(self.region_name)}:mediapackage:channel:{channel_id}"
        channel = Channel(
            arn=arn,
            description=description,
            channel_id=channel_id,
            tags=tags,
        )
        self._channels[channel_id] = channel
        return channel

    def list_channels(self) -> List[Dict[str, Any]]:
        return [c.to_dict() for c in self._channels.values()]

    def describe_channel(self, channel_id: str) -> Channel:
        try:
            return self._channels[channel_id]
        except KeyError:
            raise ClientError(
                "NotFoundException", f"channel with id={channel_id} not found"
            )

    def delete_channel(self, channel_id: str) -> Channel:
        if channel_id in self._channels:
            return self._channels.pop(channel_id)
        raise ClientError(
            "NotFoundException", f"channel with id={channel_id} not found"
        )

    def create_origin_endpoint(
        self,
        authorization: Dict[str, str],
        channel_id: str,
        cmaf_package: Dict[str, Any],
        dash_package: Dict[str, Any],
        description: str,
        hls_package: Dict[str, Any],
        endpoint_id: str,
        manifest_name: str,
        mss_package: Dict[str, Any],
        origination: str,
        startover_window_seconds: int,
        tags: Dict[str, str],
        time_delay_seconds: int,
        whitelist: List[str],
    ) -> OriginEndpoint:
        arn = f"arn:{get_partition(self.region_name)}:mediapackage:origin_endpoint:{endpoint_id}"
        url = f"https://origin-endpoint.mediapackage.{self.region_name}.amazonaws.com/{endpoint_id}"
        origin_endpoint = OriginEndpoint(
            arn=arn,
            authorization=authorization,
            channel_id=channel_id,
            cmaf_package=cmaf_package,
            dash_package=dash_package,
            description=description,
            hls_package=hls_package,
            endpoint_id=endpoint_id,
            manifest_name=manifest_name,
            mss_package=mss_package,
            origination=origination,
            startover_window_seconds=startover_window_seconds,
            tags=tags,
            time_delay_seconds=time_delay_seconds,
            url=url,
            whitelist=whitelist,
        )
        self._origin_endpoints[endpoint_id] = origin_endpoint
        return origin_endpoint

    def describe_origin_endpoint(self, endpoint_id: str) -> OriginEndpoint:
        try:
            return self._origin_endpoints[endpoint_id]
        except KeyError:
            raise ClientError(
                "NotFoundException", f"origin endpoint with id={endpoint_id} not found"
            )

    def list_origin_endpoints(self) -> List[Dict[str, Any]]:
        return [o.to_dict() for o in self._origin_endpoints.values()]

    def delete_origin_endpoint(self, endpoint_id: str) -> OriginEndpoint:
        if endpoint_id in self._origin_endpoints:
            return self._origin_endpoints.pop(endpoint_id)
        raise ClientError(
            "NotFoundException", f"origin endpoint with id={endpoint_id} not found"
        )

    def update_origin_endpoint(
        self,
        authorization: Dict[str, str],
        cmaf_package: Dict[str, Any],
        dash_package: Dict[str, Any],
        description: str,
        hls_package: Dict[str, Any],
        endpoint_id: str,
        manifest_name: str,
        mss_package: Dict[str, Any],
        origination: str,
        startover_window_seconds: int,
        time_delay_seconds: int,
        whitelist: List[str],
    ) -> OriginEndpoint:
        try:
            origin_endpoint = self._origin_endpoints[endpoint_id]
            origin_endpoint.authorization = authorization
            origin_endpoint.cmaf_package = cmaf_package
            origin_endpoint.dash_package = dash_package
            origin_endpoint.description = description
            origin_endpoint.hls_package = hls_package
            origin_endpoint.manifest_name = manifest_name
            origin_endpoint.mss_package = mss_package
            origin_endpoint.origination = origination
            origin_endpoint.startover_window_seconds = startover_window_seconds
            origin_endpoint.time_delay_seconds = time_delay_seconds
            origin_endpoint.whitelist = whitelist
            return origin_endpoint
        except KeyError:
            raise ClientError(
                "NotFoundException", f"origin endpoint with id={endpoint_id} not found"
            )


mediapackage_backends = BackendDict(MediaPackageBackend, "mediapackage")
