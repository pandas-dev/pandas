from typing import Tuple, Union

import requests

from ..models import Integration
from . import IntegrationParser


class TypeAwsParser(IntegrationParser):
    def invoke(
        self, request: requests.PreparedRequest, integration: Integration
    ) -> Tuple[int, Union[str, bytes]]:
        # integration.uri = arn:aws:apigateway:{region}:{subdomain.service|service}:path|action/{service_api}
        # example value = 'arn:aws:apigateway:us-west-2:dynamodb:action/PutItem'
        try:
            # We need a better way to support services automatically
            # This is how AWS does it though - sending a new HTTP request to the target service
            arn, action = integration.uri.split("/")
            _, _, _, region, service, path_or_action = arn.split(":")
            if service == "dynamodb" and path_or_action == "action":
                target_url = f"https://dynamodb.{region}.amazonaws.com/"
                headers = {"X-Amz-Target": f"DynamoDB_20120810.{action}"}
                res = requests.post(target_url, request.body, headers=headers)
                return res.status_code, res.content
            else:
                return (
                    400,
                    f"Integration for service {service} / {path_or_action} is not yet supported",
                )
        except Exception as e:
            return 400, str(e)
