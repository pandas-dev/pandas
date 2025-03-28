import json
import string
from typing import Any, Dict, List, Union

import yaml

from moto.moto_api._internal import mock_random as random
from moto.utilities.id_generator import ResourceIdentifier, Tags, generate_str_id


class ApigwIdentifier(ResourceIdentifier):
    service = "apigateway"

    def __init__(self, account_id: str, region: str, name: str):
        super().__init__(account_id, region, name)

    def generate(
        self, existing_ids: Union[List[str], None] = None, tags: Tags = None
    ) -> str:
        return generate_str_id(
            resource_identifier=self,
            existing_ids=existing_ids,
            tags=tags,
            length=10,
            include_digits=True,
            lower_case=True,
        )


class ApigwApiKeyIdentifier(ApigwIdentifier):
    resource = "api_key"

    def __init__(self, account_id: str, region: str, value: str):
        super().__init__(account_id, region, value)


class ApigwAuthorizerIdentifier(ApigwIdentifier):
    resource = "authorizer"


class ApigwDeploymentIdentifier(ApigwIdentifier):
    resource = "deployment"

    def __init__(self, account_id: str, region: str, stage_name: str):
        super().__init__(account_id, region, stage_name)


class ApigwModelIdentifier(ApigwIdentifier):
    resource = "model"


class ApigwRequestValidatorIdentifier(ApigwIdentifier):
    resource = "request_validator"


class ApigwResourceIdentifier(ApigwIdentifier):
    resource = "resource"

    def __init__(
        self, account_id: str, region: str, parent_id: str = "", path_name: str = "/"
    ):
        super().__init__(
            account_id,
            region,
            ".".join((parent_id, path_name)),
        )


class ApigwRestApiIdentifier(ApigwIdentifier):
    resource = "rest_api"


class ApigwUsagePlanIdentifier(ApigwIdentifier):
    resource = "usage_plan"


class ApigwVpcLinkIdentifier(ApigwIdentifier):
    resource = "vpc_link"


def create_id() -> str:
    size = 10
    chars = list(range(10)) + list(string.ascii_lowercase)
    return "".join(str(random.choice(chars)) for x in range(size))


def deserialize_body(body: str) -> Dict[str, Any]:
    try:
        api_doc = json.loads(body)
    except json.JSONDecodeError:
        api_doc = yaml.safe_load(body)

    if "openapi" in api_doc or "swagger" in api_doc:
        return api_doc

    return {}


def to_path(prop: str) -> str:
    return "/" + prop
