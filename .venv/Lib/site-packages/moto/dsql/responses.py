"""Handles incoming dsql requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import AuroraDSQLBackend, dsql_backends


class AuroraDSQLResponse(BaseResponse):
    """Handler for AuroraDSQL requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="dsql")

    @property
    def dsql_backend(self) -> AuroraDSQLBackend:
        """Return backend instance specific for this region."""
        return dsql_backends[self.current_account][self.region]

    def create_cluster(self) -> str:
        params = self._get_params()
        deletion_protection_enabled = params.get("deletionProtectionEnabled", True)
        tags = params.get("tags")
        client_token = params.get("clientToken")
        cluster = self.dsql_backend.create_cluster(
            deletion_protection_enabled=deletion_protection_enabled,
            tags=tags,
            client_token=client_token,
        )

        return json.dumps(dict(cluster.to_dict()))

    def get_cluster(self) -> str:
        identifier = self.path.split("/")[-1]
        cluster = self.dsql_backend.get_cluster(identifier=identifier)

        return json.dumps(dict(cluster.to_dict()))
