"""BedrockRuntimeBackend class with methods for supported APIs."""

from typing import Any

from moto.core.base_backend import BackendDict, BaseBackend


class BedrockRuntimeBackend(BaseBackend):
    """Implementation of BedrockRuntime APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)

    def invoke_model(
        self,
        payload: dict[str, Any],
        model_id: str,
    ) -> dict[str, Any]:
        assert payload is not None
        assert model_id is not None
        inference_result: dict[str, Any] = {}
        return inference_result


# Using `ec2` for the service name to work around lack of regions for bedrock-runtime in Botocore.
# https://github.com/getmoto/moto/issues/7745
bedrockruntime_backends = BackendDict(BedrockRuntimeBackend, "ec2")
