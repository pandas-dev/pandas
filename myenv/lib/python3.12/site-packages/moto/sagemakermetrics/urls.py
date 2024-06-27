"""sagemakermetrics base URL and path."""

from .responses import SageMakerMetricsResponse

url_bases = [
    r"https?://metrics.sagemaker\.(.+)\.amazonaws\.com",
]

url_paths = {
    "{0}/$": SageMakerMetricsResponse.dispatch,
    "{0}/BatchPutMetrics$": SageMakerMetricsResponse.dispatch,
}
