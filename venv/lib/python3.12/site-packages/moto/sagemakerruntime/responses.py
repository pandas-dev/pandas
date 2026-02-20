import base64
import json

from moto.core.common_types import TYPE_RESPONSE
from moto.core.responses import BaseResponse
from moto.moto_api._internal import mock_random as random

from .models import SageMakerRuntimeBackend, sagemakerruntime_backends


class SageMakerRuntimeResponse(BaseResponse):
    """Handler for SageMakerRuntime requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="sagemaker-runtime")

    @property
    def sagemakerruntime_backend(self) -> SageMakerRuntimeBackend:
        """Return backend instance specific for this region."""
        return sagemakerruntime_backends[self.current_account][self.region]

    def invoke_endpoint(self) -> TYPE_RESPONSE:
        params = self._get_params()
        unique_repr = {
            key: value
            for key, value in self.headers.items()
            if key.lower().startswith("x-amzn-sagemaker")
        }
        unique_repr["Accept"] = self.headers.get("Accept")
        unique_repr["Body"] = self.body
        endpoint_name = params.get("EndpointName")
        (
            body,
            content_type,
            invoked_production_variant,
            custom_attributes,
        ) = self.sagemakerruntime_backend.invoke_endpoint(
            endpoint_name=endpoint_name,  # type: ignore[arg-type]
            unique_repr=base64.b64encode(json.dumps(unique_repr).encode("utf-8")),
        )
        headers = {"Content-Type": content_type}
        if invoked_production_variant:
            headers["x-Amzn-Invoked-Production-Variant"] = invoked_production_variant
        if custom_attributes:
            headers["X-Amzn-SageMaker-Custom-Attributes"] = custom_attributes
        return 200, headers, body

    def invoke_endpoint_async(self) -> TYPE_RESPONSE:
        endpoint_name = self.path.split("/")[2]
        input_location = self.headers.get("X-Amzn-SageMaker-InputLocation")
        inference_id = self.headers.get("X-Amzn-SageMaker-Inference-Id")
        output_location, failure_location = (
            self.sagemakerruntime_backend.invoke_endpoint_async(
                endpoint_name, input_location
            )
        )
        resp = {"InferenceId": inference_id or str(random.uuid4())}
        headers = {
            "X-Amzn-SageMaker-OutputLocation": output_location,
            "X-Amzn-SageMaker-FailureLocation": failure_location,
        }
        return 200, headers, json.dumps(resp)
