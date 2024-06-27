"""SageMakerMetricsBackend class with methods for supported APIs."""

from datetime import datetime
from typing import Dict, List, Union, cast

from moto.core.base_backend import BackendDict, BaseBackend
from moto.sagemaker import sagemaker_backends
from moto.sagemaker.models import METRIC_STEP_TYPE

RESPONSE_TYPE = Dict[str, List[Dict[str, Union[str, int]]]]


class SageMakerMetricsBackend(BaseBackend):
    """Implementation of SageMakerMetrics APIs."""

    def __init__(self, region_name: str, account_id: str) -> None:
        super().__init__(region_name, account_id)
        self.sagemaker_backend = sagemaker_backends[account_id][region_name]

    def batch_put_metrics(
        self,
        trial_component_name: str,
        metric_data: List[Dict[str, Union[str, int, float, datetime]]],
    ) -> RESPONSE_TYPE:
        return_response: RESPONSE_TYPE = {"Errors": []}

        if trial_component_name not in self.sagemaker_backend.trial_components:
            return_response["Errors"].append(
                {"Code": "VALIDATION_ERROR", "MetricIndex": 0}
            )
            return return_response

        trial_component = self.sagemaker_backend.trial_components[trial_component_name]
        for metric in metric_data:
            metric_step: int = cast(int, metric["Step"])
            metric_name: str = cast(str, metric["MetricName"])
            if metric_name not in trial_component.metrics:
                metric_timestamp: int = cast(int, metric["Timestamp"])
                values_dict: Dict[int, Dict[str, Union[str, int, float, datetime]]] = {}
                new_metric: Dict[str, Union[str, int, METRIC_STEP_TYPE]] = {
                    "MetricName": metric_name,
                    "Timestamp": metric_timestamp,
                    "Values": values_dict,
                }
                trial_component.metrics[metric_name] = new_metric
            new_step: METRIC_STEP_TYPE = {metric_step: metric}
            trial_component_metric_values: METRIC_STEP_TYPE = cast(
                METRIC_STEP_TYPE, trial_component.metrics[metric_name]["Values"]
            )
            trial_component_metric_values.update(new_step)  # type ignore
        return return_response


sagemakermetrics_backends = BackendDict(SageMakerMetricsBackend, "sagemaker-metrics")
