"""DirectConnectBackend class with methods for supported APIs."""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.utilities.tagging_service import TaggingService

from .enums import (
    ConnectionStateType,
    EncryptionModeType,
    LagStateType,
    MacSecKeyStateType,
    PortEncryptionStatusType,
)
from .exceptions import (
    ConnectionIdMissing,
    ConnectionNotFound,
    LAGNotFound,
    MacSecKeyNotFound,
)


@dataclass
class MacSecKey(BaseModel):
    secret_arn: Optional[str]
    ckn: Optional[str]
    state: MacSecKeyStateType
    start_on: str
    cak: Optional[str] = None

    def to_dict(self) -> Dict[str, str]:
        return {
            "secretARN": self.secret_arn or "",
            "ckn": self.ckn or "",
            "state": self.state,
            "startOn": self.start_on,
        }


@dataclass
class Connection(BaseModel):
    aws_device_v2: str
    aws_device: str
    aws_logical_device_id: str
    bandwidth: str
    connection_name: str
    connection_state: ConnectionStateType
    encryption_mode: EncryptionModeType
    has_logical_redundancy: bool
    jumbo_frame_capable: bool
    lag_id: Optional[str]
    loa_issue_time: str
    location: str
    mac_sec_capable: Optional[bool]
    mac_sec_keys: List[MacSecKey]
    owner_account: str
    partner_name: str
    port_encryption_status: PortEncryptionStatusType
    provider_name: Optional[str]
    region: str
    tags: Optional[List[Dict[str, str]]]
    vlan: int
    connection_id: str = field(default="", init=False)
    backend: "DirectConnectBackend"

    def __post_init__(self) -> None:
        if self.connection_id == "":
            self.connection_id = f"arn:aws:directconnect:{self.region}:{self.owner_account}:dx-con/dx-moto-{self.connection_name}-{datetime.now().strftime('%Y%m%d%H%M%S')}"

    def to_dict(
        self,
    ) -> Dict[str, Any]:
        return {
            "awsDevice": self.aws_device,
            "awsDeviceV2": self.aws_device_v2,
            "awsLogicalDeviceId": self.aws_logical_device_id,
            "bandwidth": self.bandwidth,
            "connectionId": self.connection_id,
            "connectionName": self.connection_name,
            "connectionState": self.connection_state,
            "encryptionMode": self.encryption_mode,
            "hasLogicalRedundancy": self.has_logical_redundancy,
            "jumboFrameCapable": self.jumbo_frame_capable,
            "lagId": self.lag_id,
            "loaIssueTime": self.loa_issue_time,
            "location": self.location,
            "macSecCapable": self.mac_sec_capable,
            "macSecKeys": [key.to_dict() for key in self.mac_sec_keys],
            "partnerName": self.partner_name,
            "portEncryptionStatus": self.port_encryption_status,
            "providerName": self.provider_name,
            "region": self.region,
            "tags": self.backend.list_tags_for_resource(self.connection_id),
            "vlan": self.vlan,
        }


@dataclass
class LAG(BaseModel):
    aws_device_v2: str
    aws_device: str
    aws_logical_device_id: str
    connections_bandwidth: str
    number_of_connections: int
    minimum_links: int
    connections: List[Connection]
    lag_name: str
    lag_state: LagStateType
    encryption_mode: EncryptionModeType
    has_logical_redundancy: bool
    jumbo_frame_capable: bool
    location: str
    mac_sec_capable: Optional[bool]
    mac_sec_keys: List[MacSecKey]
    owner_account: str
    provider_name: Optional[str]
    region: str
    tags: Optional[List[Dict[str, str]]]
    lag_id: str = field(default="", init=False)
    backend: "DirectConnectBackend"

    def __post_init__(self) -> None:
        if self.lag_id == "":
            self.lag_id = f"arn:aws:directconnect:{self.region}:{self.owner_account}:dxlag/dxlag-moto-{self.lag_name}-{datetime.now().strftime('%Y%m%d%H%M%S')}"

    def to_dict(
        self,
    ) -> Dict[str, Any]:
        return {
            "awsDevice": self.aws_device,
            "awsDeviceV2": self.aws_device_v2,
            "awsLogicalDeviceId": self.aws_logical_device_id,
            "connectionsBandwidth": self.connections_bandwidth,
            "numberOfConnections": self.number_of_connections,
            "minimumLinks": self.minimum_links,
            "connections": [conn.to_dict() for conn in self.connections],
            "lagId": self.lag_id,
            "lagName": self.lag_name,
            "lagState": self.lag_state,
            "encryptionMode": self.encryption_mode,
            "hasLogicalRedundancy": self.has_logical_redundancy,
            "jumboFrameCapable": self.jumbo_frame_capable,
            "location": self.location,
            "macSecCapable": self.mac_sec_capable,
            "macSecKeys": [key.to_dict() for key in self.mac_sec_keys],
            "providerName": self.provider_name,
            "region": self.region,
            "tags": self.backend.list_tags_for_resource(self.lag_id),
        }


class DirectConnectBackend(BaseBackend):
    """Implementation of DirectConnect APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.connections: Dict[str, Connection] = {}
        self.lags: Dict[str, LAG] = {}
        self.tagger = TaggingService(key_name="key", value_name="value")

    def describe_connections(self, connection_id: Optional[str]) -> List[Connection]:
        if connection_id and connection_id not in self.connections:
            raise ConnectionNotFound(connection_id, self.region_name)
        if connection_id:
            connection = self.connections.get(connection_id)
            return [] if not connection else [connection]
        return list(self.connections.values())

    def create_connection(
        self,
        location: str,
        bandwidth: str,
        connection_name: str,
        lag_id: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        provider_name: Optional[str],
        request_mac_sec: Optional[bool],
    ) -> Connection:
        encryption_mode = EncryptionModeType.NO
        mac_sec_keys = []
        if request_mac_sec:
            encryption_mode = EncryptionModeType.MUST
            mac_sec_keys = [
                MacSecKey(
                    secret_arn="mock_secret_arn",
                    ckn="mock_ckn",
                    state=MacSecKeyStateType.ASSOCIATED,
                    start_on=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                )
            ]
        connection = Connection(
            aws_device_v2="mock_device_v2",
            aws_device="mock_device",
            aws_logical_device_id="mock_logical_device_id",
            bandwidth=bandwidth,
            connection_name=connection_name,
            connection_state=ConnectionStateType.AVAILABLE,
            encryption_mode=encryption_mode,
            has_logical_redundancy=False,
            jumbo_frame_capable=False,
            lag_id=lag_id,
            loa_issue_time=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            location=location,
            mac_sec_capable=request_mac_sec,
            mac_sec_keys=mac_sec_keys,
            owner_account=self.account_id,
            partner_name="mock_partner",
            port_encryption_status=PortEncryptionStatusType.DOWN,
            provider_name=provider_name,
            region=self.region_name,
            tags=tags,
            vlan=0,
            backend=self,
        )
        if tags:
            self.tag_resource(connection.connection_id, tags)
        self.connections[connection.connection_id] = connection
        return connection

    def tag_resource(self, resource_arn: str, tags: List[Dict[str, str]]) -> None:
        self.tagger.tag_resource(
            resource_arn,
            tags=tags if tags else [],
        )

    def list_tags_for_resource(self, resource_arn: str) -> List[Dict[str, str]]:
        tags = self.tagger.get_tag_dict_for_resource(resource_arn)
        if not tags:
            return []
        return [{"key": k, "value": v} for (k, v) in tags.items()]

    def list_tags_for_resources(self, resource_arns: List[str]) -> Dict[str, List[Any]]:
        response: Dict[str, List[Any]] = {"resourceTags": []}
        for resource_arn in resource_arns:
            response["resourceTags"].append(
                {
                    "resourceArn": resource_arn,
                    "tags": self.list_tags_for_resource(resource_arn),
                }
            )

        return response

    def untag_resource(self, resource_arn: str, tag_keys: List[str]) -> None:
        self.tagger.untag_resource_using_names(resource_arn, tag_keys)

    def delete_connection(self, connection_id: str) -> Connection:
        if not connection_id:
            raise ConnectionIdMissing()
        connection = self.connections.get(connection_id)
        if connection:
            self.connections[
                connection_id
            ].connection_state = ConnectionStateType.DELETED
            return connection
        raise ConnectionNotFound(connection_id, self.region_name)

    def update_connection(
        self,
        connection_id: str,
        new_connection_name: Optional[str],
        new_encryption_mode: Optional[EncryptionModeType],
    ) -> Connection:
        if not connection_id:
            raise ConnectionIdMissing()
        connection = self.connections.get(connection_id)
        if connection:
            if new_connection_name:
                self.connections[connection_id].connection_name = new_connection_name
            if new_encryption_mode:
                self.connections[connection_id].encryption_mode = new_encryption_mode
            return connection
        raise ConnectionNotFound(connection_id, self.region_name)

    def associate_mac_sec_key(
        self,
        connection_id: str,
        secret_arn: Optional[str],
        ckn: Optional[str],
        cak: Optional[str],
    ) -> Tuple[str, List[MacSecKey]]:
        mac_sec_key = MacSecKey(
            secret_arn=secret_arn or "mock_secret_arn",
            ckn=ckn,
            cak=cak,
            state=MacSecKeyStateType.ASSOCIATED,
            start_on=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

        if "dxlag-" in connection_id:
            return self._associate_mac_sec_key_with_lag(
                lag_id=connection_id, mac_sec_key=mac_sec_key
            )

        return self._associate_mac_sec_key_with_connection(
            connection_id=connection_id, mac_sec_key=mac_sec_key
        )

    def _associate_mac_sec_key_with_lag(
        self, lag_id: str, mac_sec_key: MacSecKey
    ) -> Tuple[str, List[MacSecKey]]:
        lag = self.lags.get(lag_id) or None
        if not lag:
            raise LAGNotFound(lag_id, self.region_name)

        lag.mac_sec_keys.append(mac_sec_key)
        for connection in lag.connections:
            connection.mac_sec_keys = lag.mac_sec_keys

        return lag_id, lag.mac_sec_keys

    def _associate_mac_sec_key_with_connection(
        self, connection_id: str, mac_sec_key: MacSecKey
    ) -> Tuple[str, List[MacSecKey]]:
        connection = self.connections.get(connection_id) or None
        if not connection:
            raise ConnectionNotFound(connection_id, self.region_name)

        self.connections[connection_id].mac_sec_keys.append(mac_sec_key)
        return connection_id, self.connections[connection_id].mac_sec_keys

    def create_lag(
        self,
        number_of_connections: int,
        location: str,
        connections_bandwidth: str,
        lag_name: str,
        connection_id: Optional[str],
        tags: Optional[List[Dict[str, str]]],
        child_connection_tags: Optional[List[Dict[str, str]]],
        provider_name: Optional[str],
        request_mac_sec: Optional[bool],
    ) -> LAG:
        if connection_id:
            raise NotImplementedError(
                "creating a lag with a connection_id is not currently supported by moto"
            )
        encryption_mode = EncryptionModeType.NO
        mac_sec_keys = []
        if request_mac_sec:
            encryption_mode = EncryptionModeType.MUST
            mac_sec_keys = [
                MacSecKey(
                    secret_arn="mock_secret_arn",
                    ckn="mock_ckn",
                    state=MacSecKeyStateType.ASSOCIATED,
                    start_on=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                )
            ]
        lag = LAG(
            aws_device_v2="mock_device_v2",
            aws_device="mock_device",
            aws_logical_device_id="mock_logical_device_id",
            connections_bandwidth=connections_bandwidth,
            lag_name=lag_name,
            lag_state=LagStateType.AVAILABLE,
            minimum_links=0,
            encryption_mode=encryption_mode,
            has_logical_redundancy=False,
            jumbo_frame_capable=False,
            number_of_connections=number_of_connections,
            connections=[],
            location=location,
            mac_sec_capable=request_mac_sec,
            mac_sec_keys=mac_sec_keys,
            owner_account=self.account_id,
            provider_name=provider_name,
            region=self.region_name,
            tags=tags,
            backend=self,
        )
        for i in range(number_of_connections):
            connection = self.create_connection(
                location=location,
                bandwidth=connections_bandwidth,
                connection_name=f"Requested Connection {i + 1} for Lag {lag.lag_id}",
                lag_id=lag.lag_id,
                tags=child_connection_tags,
                request_mac_sec=False,
                provider_name=provider_name,
            )
            if request_mac_sec:
                connection.mac_sec_capable = True
                connection.mac_sec_keys = mac_sec_keys
                connection.encryption_mode = encryption_mode
            lag.connections.append(connection)

        if tags:
            self.tag_resource(lag.lag_id, tags)
        self.lags[lag.lag_id] = lag
        return lag

    def describe_lags(self, lag_id: Optional[str]) -> List[LAG]:
        if lag_id and lag_id not in self.lags:
            raise LAGNotFound(lag_id, self.region_name)
        if lag_id:
            lag = self.lags.get(lag_id)
            return [] if not lag else [lag]
        return list(self.lags.values())

    def disassociate_mac_sec_key(
        self, connection_id: str, secret_arn: str
    ) -> Tuple[str, MacSecKey]:
        mac_sec_keys: List[MacSecKey] = []
        if "dxlag-" in connection_id:
            if connection_id in self.lags:
                mac_sec_keys = self.lags[connection_id].mac_sec_keys
        elif connection_id in self.connections:
            mac_sec_keys = self.connections[connection_id].mac_sec_keys
        if not mac_sec_keys:
            raise ConnectionNotFound(connection_id, self.region_name)

        arn_casefold = secret_arn.casefold()
        for i, mac_sec_key in enumerate(mac_sec_keys):
            if str(mac_sec_key.secret_arn).casefold() == arn_casefold:
                mac_sec_key.state = MacSecKeyStateType.DISASSOCIATED
                return connection_id, mac_sec_keys.pop(i)

        raise MacSecKeyNotFound(secret_arn=secret_arn, connection_id=connection_id)


directconnect_backends = BackendDict(DirectConnectBackend, "directconnect")
