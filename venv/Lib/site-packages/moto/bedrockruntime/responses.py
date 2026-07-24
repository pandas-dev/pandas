"""Handles incoming bedrockruntime requests, invokes methods, returns responses."""

import json
from typing import Any

from moto.core.responses import ActionResult, BaseResponse
from moto.utilities.constants import APPLICATION_JSON

from .models import BedrockRuntimeBackend, bedrockruntime_backends


class BedrockRuntimeResponse(BaseResponse):
    """Handler for BedrockRuntime requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="bedrock-runtime")
        self.automated_parameter_parsing = True

    @property
    def bedrockruntime_backend(self) -> BedrockRuntimeBackend:
        """Return backend instance specific for this region."""
        return bedrockruntime_backends[self.current_account][self.region]

    def invoke_model(self) -> ActionResult:
        payload = self._get_param("body")
        content_type = self._get_param("contentType", APPLICATION_JSON)
        if content_type == APPLICATION_JSON:
            payload = json.loads(payload)
        accept = self._get_param("accept", APPLICATION_JSON)
        model_id = self._get_param("modelId")
        performance_config_latency = self._get_param(
            "performanceConfigLatency", "standard"
        )
        service_tier = self._get_param("serviceTier", "default")
        inference_result = self.bedrockruntime_backend.invoke_model(
            payload=payload,
            model_id=model_id,
        )
        body: Any = inference_result
        if accept == APPLICATION_JSON:
            body = json.dumps(inference_result)
        result = {
            "body": body,
            "contentType": "application/json",
            "performanceConfigLatency": performance_config_latency,
            "serviceTier": service_tier,
        }
        return ActionResult(result)
