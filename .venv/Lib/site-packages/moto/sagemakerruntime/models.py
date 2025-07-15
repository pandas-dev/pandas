import json
from typing import Dict, List, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.moto_api._internal import mock_random as random


class SageMakerRuntimeBackend(BaseBackend):
    """Implementation of SageMakerRuntime APIs."""

    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.async_results: Dict[str, Dict[str, Tuple[str, str]]] = {}
        self.results: Dict[str, Dict[bytes, Tuple[str, str, str, str]]] = {}
        self.results_queue: List[Tuple[str, str, str, str]] = []
        self.async_results_queue: List[Tuple[bool, str]] = []

    def invoke_endpoint(
        self, endpoint_name: str, unique_repr: bytes
    ) -> Tuple[str, str, str, str]:
        """
        This call will return static data by default.

        You can use a dedicated API to override this, by configuring a queue of expected results.

        A request to `get_query_results` will take the first result from that queue. Subsequent requests using the same details will return the same result. Other requests using a different QueryExecutionId will take the next result from the queue, or return static data if the queue is empty.

        Configuring this queue by making an HTTP request to `/moto-api/static/sagemaker/endpoint-results`. An example invocation looks like this:

        .. sourcecode:: python

            expected_results = {
                "account_id": "123456789012",  # This is the default - can be omitted
                "region": "us-east-1",  # This is the default - can be omitted
                "results": [
                    {
                         "Body": "first body",
                         "ContentType": "text/xml",
                         "InvokedProductionVariant": "prod",
                         "CustomAttributes": "my_attr",
                     },
                    # other results as required
                ],
            }
            requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/sagemaker/endpoint-results",
                json=expected_results,
            )

            client = boto3.client("sagemaker-runtime", region_name="us-east-1")
            details = client.invoke_endpoint(EndpointName="asdf", Body="qwer")

        """
        if endpoint_name not in self.results:
            self.results[endpoint_name] = {}
        if unique_repr in self.results[endpoint_name]:
            return self.results[endpoint_name][unique_repr]
        if self.results_queue:
            self.results[endpoint_name][unique_repr] = self.results_queue.pop(0)
        else:
            self.results[endpoint_name][unique_repr] = (
                "body",
                "content_type",
                "invoked_production_variant",
                "custom_attributes",
            )
        return self.results[endpoint_name][unique_repr]

    def invoke_endpoint_async(
        self, endpoint_name: str, input_location: str
    ) -> Tuple[str, str]:
        """
        This call will return static data by default.

        You can use a dedicated API to override this, by configuring a queue of expected results.
        Subsequent requests using the same `InputLocation` will return the same result.
        Other requests using a different `InputLocation` will take the next result from the queue.
        Please note that for backward compatibility, if the async queue is empty, the sync queue will be used, if configured.

        Configuring this queue by making an HTTP request to `/moto-api/static/sagemaker/async-endpoint-results`.
        An example invocation looks like this:

        .. sourcecode:: python

            expected_results = {
                "account_id": "123456789012",  # This is the default - can be omitted
                "region": "us-east-1",  # This is the default - can be omitted
                "results": [
                    {
                        "data": json.dumps({"first": "output"}),
                    },
                    {
                        "is_failure": True,
                        "data": "second inference failed",
                    },
                    # other results as required
                ],
            }
            requests.post(
                "http://motoapi.amazonaws.com/moto-api/static/sagemaker/async-endpoint-results",
                json=expected_results,
            )

            client = boto3.client("sagemaker-runtime", region_name="us-east-1")
            details = client.invoke_endpoint_async(EndpointName="asdf", InputLocation="qwer")

        """
        if endpoint_name not in self.async_results:
            self.async_results[endpoint_name] = {}
        if input_location in self.async_results[endpoint_name]:
            return self.async_results[endpoint_name][input_location]
        if self.async_results_queue:
            is_failure, data = self.async_results_queue.pop(0)
        elif self.results_queue:  # Backward compatibility
            is_failure = False
            body, _type, variant, attrs = self.results_queue.pop(0)
            data = json.dumps(
                {
                    "Body": body,
                    "ContentType": _type,
                    "InvokedProductionVariant": variant,
                    "CustomAttributes": attrs,
                }
            )
        else:
            is_failure = False
            data = json.dumps({"default": "response"})

        output_bucket = f"sagemaker-output-{random.uuid4()}"
        output_location = "response.json"
        failure_location = "failure.json"
        result_location = failure_location if is_failure else output_location
        from moto.s3.models import s3_backends

        s3_backend = s3_backends[self.account_id][self.partition]
        s3_backend.create_bucket(output_bucket, region_name=self.region_name)
        s3_backend.put_object(
            output_bucket, result_location, value=data.encode("utf-8")
        )

        self.async_results[endpoint_name][input_location] = (
            f"s3://{output_bucket}/{output_location}",
            f"s3://{output_bucket}/{failure_location}",
        )
        return self.async_results[endpoint_name][input_location]


sagemakerruntime_backends = BackendDict(SageMakerRuntimeBackend, "sagemaker-runtime")
