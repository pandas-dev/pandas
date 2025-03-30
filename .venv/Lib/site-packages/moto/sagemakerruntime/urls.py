"""sagemakerruntime base URL and path."""

from .responses import SageMakerRuntimeResponse

url_bases = [
    r"https?://runtime\.sagemaker\.(.+)\.amazonaws\.com",
]


response = SageMakerRuntimeResponse()


url_paths = {
    "{0}/endpoints/(?P<name>[^/]+)/async-invocations$": response.dispatch,
    "{0}/endpoints/(?P<name>[^/]+)/invocations$": response.dispatch,
}
