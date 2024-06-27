import json
import threading
from typing import Any, Dict

from moto import settings
from moto.awslambda import lambda_backends
from moto.core.common_models import CloudFormationModel
from moto.moto_api._internal import mock_random


class CustomModel(CloudFormationModel):
    def __init__(
        self, region_name: str, request_id: str, logical_id: str, resource_name: str
    ):
        self.region_name = region_name
        self.request_id = request_id
        self.logical_id = logical_id
        self.resource_name = resource_name
        self.data: Dict[str, Any] = dict()
        self._finished = False

    def set_data(self, data: Dict[str, Any]) -> None:
        self.data = data
        self._finished = True

    def is_created(self) -> bool:
        return self._finished

    @property
    def physical_resource_id(self) -> str:
        return self.resource_name

    @staticmethod
    def cloudformation_type() -> str:
        return "?"

    @classmethod
    def create_from_cloudformation_json(  # type: ignore[misc]
        cls,
        resource_name: str,
        cloudformation_json: Dict[str, Any],
        account_id: str,
        region_name: str,
        **kwargs: Any,
    ) -> "CustomModel":
        logical_id = kwargs["LogicalId"]
        stack_id = kwargs["StackId"]
        resource_type = kwargs["ResourceType"]
        properties = cloudformation_json["Properties"]
        service_token = properties["ServiceToken"]

        backend = lambda_backends[account_id][region_name]
        fn = backend.get_function(service_token)

        request_id = str(mock_random.uuid4())

        custom_resource = CustomModel(
            region_name, request_id, logical_id, resource_name
        )

        from moto.cloudformation import cloudformation_backends

        stack = cloudformation_backends[account_id][region_name].get_stack(stack_id)
        stack.add_custom_resource(custom_resource)

        # A request will be send to this URL to indicate success/failure
        # This request will be coming from inside a Docker container
        # Note that, in order to reach the Moto host, the Moto-server should be listening on 0.0.0.0
        #
        # Alternative: Maybe we should let the user pass in a container-name where Moto is running?
        # Similar to how we know for sure that the container in our CI is called 'motoserver'
        host = f"{settings.moto_server_host()}:{settings.moto_server_port()}"
        response_url = (
            f"{host}/cloudformation_{region_name}/cfnresponse?stack={stack_id}"
        )

        event = {
            "RequestType": "Create",
            "ServiceToken": service_token,
            "ResponseURL": response_url,
            "StackId": stack_id,
            "RequestId": request_id,
            "LogicalResourceId": logical_id,
            "ResourceType": resource_type,
            "ResourceProperties": properties,
        }

        invoke_thread = threading.Thread(
            target=fn.invoke, args=(json.dumps(event), {}, {})
        )
        invoke_thread.start()

        return custom_resource

    @classmethod
    def has_cfn_attr(cls, attr: str) -> bool:  # pylint: disable=unused-argument
        # We don't know which attributes are supported for third-party resources
        return True

    def get_cfn_attribute(self, attribute_name: str) -> Any:
        if attribute_name in self.data:
            return self.data[attribute_name]
        return None
