from typing import Any

from moto.core.exceptions import RESTError, ServiceException

ERROR_WITH_ACCESS_POINT_NAME = """{% extends 'wrapped_single_error' %}
{% block extra %}<AccessPointName>{{ name }}</AccessPointName>{% endblock %}
"""


ERROR_WITH_ACCESS_POINT_POLICY = """{% extends 'wrapped_single_error' %}
{% block extra %}<AccessPointName>{{ name }}</AccessPointName>{% endblock %}
"""

ERROR_WITH_MRAP_NAME = """{% extends 'wrapped_single_error' %}
{% block extra %}<Name>{{ name }}</Name>{% endblock %}
"""


ERROR_WITH_REQUEST_TOKEN = """{% extends 'wrapped_single_error' %}
{% block extra %}<RequestTokenARN>{{ request_token }}</RequestTokenARN>{% endblock %}
"""


class S3ControlError(RESTError):
    extended_templates = {
        "ap_not_found": ERROR_WITH_ACCESS_POINT_NAME,
        "apf_not_found": ERROR_WITH_ACCESS_POINT_POLICY,
        "mrap_not_found": ERROR_WITH_MRAP_NAME,
        "mrap_policy_not_found": ERROR_WITH_MRAP_NAME,
        "mrap_operation_not_found": ERROR_WITH_REQUEST_TOKEN,
    }
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "single_error")
        super().__init__(*args, **kwargs)


class AccessPointNotFound(S3ControlError):
    code = 404

    def __init__(self, name: str, **kwargs: Any):
        kwargs.setdefault("template", "ap_not_found")
        kwargs["name"] = name
        super().__init__(
            "NoSuchAccessPoint", "The specified accesspoint does not exist", **kwargs
        )


class AccessPointPolicyNotFound(S3ControlError):
    code = 404

    def __init__(self, name: str, **kwargs: Any):
        kwargs.setdefault("template", "apf_not_found")
        kwargs["name"] = name
        super().__init__(
            "NoSuchAccessPointPolicy",
            "The specified accesspoint policy does not exist",
            **kwargs,
        )


class MultiRegionAccessPointNotFound(S3ControlError):
    code = 404

    def __init__(self, name: str, **kwargs: Any):
        kwargs.setdefault("template", "mrap_not_found")
        kwargs["name"] = name
        super().__init__(
            "NoSuchMultiRegionAccessPoint",
            "The specified multi-region access point does not exist",
            **kwargs,
        )


class MultiRegionAccessPointPolicyNotFound(S3ControlError):
    code = 404

    def __init__(self, name: str, **kwargs: Any):
        kwargs.setdefault("template", "mrap_policy_not_found")
        kwargs["name"] = name
        super().__init__(
            "NoSuchMultiRegionAccessPointPolicy",
            "The specified multi-region access point policy does not exist",
            **kwargs,
        )


class MultiRegionAccessPointOperationNotFound(S3ControlError):
    code = 404

    def __init__(self, request_token: str, **kwargs: Any):
        kwargs.setdefault("template", "mrap_operation_not_found")
        kwargs["request_token"] = request_token
        super().__init__(
            "NoSuchAsyncRequest",
            "The specified async request does not exist",
            **kwargs,
        )


class NoSuchPublicAccessBlockConfiguration(ServiceException):
    # Note that this exception is in the different format then the S3 exception with the same name
    # This exception should return a nested response `<ErrorResponse><Error>..`
    # The S3 variant uses a flat `<Error>`-response
    code = "NoSuchPublicAccessBlockConfiguration"

    def __init__(self) -> None:
        super().__init__("The public access block configuration was not found")


class InvalidRequestException(S3ControlError):
    code = 400

    def __init__(self, message: str, **kwargs: Any):
        super().__init__(
            "InvalidRequest",
            message,
            **kwargs,
        )


class StorageLensConfigurationNotFound(S3ControlError):
    code = 404

    def __init__(self, config_id: str, **kwargs: Any):
        super().__init__(
            "NoSuchConfiguration",
            f"The specified configuration does not exist: {config_id}",
            **kwargs,
        )
