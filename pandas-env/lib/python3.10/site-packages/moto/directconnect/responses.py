"""Handles incoming directconnect requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import LAG, Connection, DirectConnectBackend, directconnect_backends


class DirectConnectResponse(BaseResponse):
    """Handler for DirectConnect requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="directconnect")

    @property
    def directconnect_backend(self) -> DirectConnectBackend:
        return directconnect_backends[self.current_account][self.region]

    def describe_connections(self) -> str:
        params = json.loads(self.body)
        connections = self.directconnect_backend.describe_connections(
            connection_id=params.get("connectionId"),
        )
        return json.dumps(
            dict(connections=[connection.to_dict() for connection in connections])
        )

    def create_connection(self) -> str:
        params = json.loads(self.body)
        connection: Connection = self.directconnect_backend.create_connection(
            location=params.get("location"),
            bandwidth=params.get("bandwidth"),
            connection_name=params.get("connectionName"),
            lag_id=params.get("lagId"),
            tags=params.get("tags"),
            provider_name=params.get("providerName"),
            request_mac_sec=params.get("requestMACSec"),
        )
        return json.dumps(connection.to_dict())

    def delete_connection(self) -> str:
        params = json.loads(self.body)
        connection: Connection = self.directconnect_backend.delete_connection(
            connection_id=params.get("connectionId"),
        )
        return json.dumps(connection.to_dict())

    def update_connection(self) -> str:
        params = json.loads(self.body)
        connection: Connection = self.directconnect_backend.update_connection(
            connection_id=params.get("connectionId"),
            new_connection_name=params.get("connectionName"),
            new_encryption_mode=params.get("encryptionMode"),
        )
        return json.dumps(connection.to_dict())

    def associate_mac_sec_key(self) -> str:
        params = json.loads(self.body)
        connection_id = params.get("connectionId")
        secret_arn = params.get("secretARN")
        ckn = params.get("ckn")
        cak = params.get("cak")
        connection_id, mac_sec_keys = self.directconnect_backend.associate_mac_sec_key(
            connection_id=connection_id,
            secret_arn=secret_arn,
            ckn=ckn,
            cak=cak,
        )
        return json.dumps(
            dict(
                connectionId=connection_id,
                macSecKeys=[mac_sec_key.to_dict() for mac_sec_key in mac_sec_keys],
            )
        )

    def create_lag(self) -> str:
        params = json.loads(self.body)
        number_of_connections = params.get("numberOfConnections")
        location = params.get("location")
        connections_bandwidth = params.get("connectionsBandwidth")
        lag_name = params.get("lagName")
        connection_id = params.get("connectionId")
        tags = params.get("tags")
        child_connection_tags = params.get("childConnectionTags")
        provider_name = params.get("providerName")
        request_mac_sec = params.get("requestMACSec")
        lag: LAG = self.directconnect_backend.create_lag(
            number_of_connections=number_of_connections,
            location=location,
            connections_bandwidth=connections_bandwidth,
            lag_name=lag_name,
            connection_id=connection_id,
            tags=tags,
            child_connection_tags=child_connection_tags,
            provider_name=provider_name,
            request_mac_sec=request_mac_sec,
        )
        return json.dumps(lag.to_dict())

    def describe_lags(self) -> str:
        params = json.loads(self.body)
        lags = self.directconnect_backend.describe_lags(
            lag_id=params.get("lagId"),
        )
        return json.dumps(dict(lags=[lag.to_dict() for lag in lags]))

    def disassociate_mac_sec_key(self) -> str:
        params = json.loads(self.body)
        connection_id = params.get("connectionId")
        secret_arn = params.get("secretARN")
        connection_id, mac_sec_key = (
            self.directconnect_backend.disassociate_mac_sec_key(
                connection_id=connection_id,
                secret_arn=secret_arn,
            )
        )
        return json.dumps(
            dict(
                connectionId=connection_id,
                macSecKeys=[mac_sec_key.to_dict()],
            )
        )
