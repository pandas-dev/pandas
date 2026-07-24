"""AWS Config integration for SNS resources."""

import json
from typing import Any

from moto.core.common_models import ConfigQueryModel
from moto.sns import sns_backends
from moto.sns.models import SNSBackend


class SNSConfigQuery(ConfigQueryModel[SNSBackend]):
    """Config query model for SNS topics."""

    def list_config_service_resources(
        self,
        account_id: str,
        partition: str,
        resource_ids: list[str] | None,
        resource_name: str | None,
        limit: int,
        next_token: str | None,
        backend_region: str | None = None,
        resource_region: str | None = None,
        aggregator: dict[str, Any] | None = None,
    ) -> tuple[list[dict[str, Any]], str | None]:
        region = resource_region or backend_region or "us-east-1"
        backend = self.backends[account_id][region]

        topics, new_token = backend.list_config_service_resources(
            resource_ids=resource_ids,
            resource_name=resource_name,
            limit=limit,
            next_token=next_token,
        )

        return topics, new_token

    def get_config_resource(
        self,
        account_id: str,
        partition: str,
        resource_id: str,
        resource_name: str | None = None,
        backend_region: str | None = None,
        resource_region: str | None = None,
    ) -> dict[str, Any] | None:
        region = resource_region or backend_region or "us-east-1"
        backend = self.backends[account_id][region]

        config_data = backend.get_config_resource(resource_id)

        if not config_data:
            return None

        if resource_name and config_data.get("resourceName") != resource_name:
            return None

        if "configuration" in config_data and not isinstance(
            config_data["configuration"], str
        ):
            config_data["configuration"] = json.dumps(config_data["configuration"])

        if "supplementaryConfiguration" in config_data:
            for field, value in config_data["supplementaryConfiguration"].items():
                if not isinstance(value, str):
                    config_data["supplementaryConfiguration"][field] = json.dumps(value)

        return config_data


sns_config_query = SNSConfigQuery(sns_backends)
