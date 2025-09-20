"""Handles incoming sagemakermetrics requests, invokes methods, returns responses."""

import json

from moto.core.responses import BaseResponse

from .models import SageMakerMetricsBackend, sagemakermetrics_backends


class SageMakerMetricsResponse(BaseResponse):
    """Handler for SageMakerMetrics requests and responses."""

    def __init__(self) -> None:
        super().__init__(service_name="sagemaker-metrics")

    @property
    def sagemakermetrics_backend(self) -> SageMakerMetricsBackend:
        """Return backend instance specific for this region."""
        return sagemakermetrics_backends[self.current_account][self.region]

    def batch_put_metrics(self) -> str:
        trial_component_name = self._get_param("TrialComponentName")
        metric_data = self._get_param("MetricData")
        errors = self.sagemakermetrics_backend.batch_put_metrics(
            trial_component_name=trial_component_name,
            metric_data=metric_data,
        )
        return json.dumps(errors)
