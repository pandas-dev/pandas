import json
from typing import Any, Optional

from jinja2 import DictLoader, Environment
from werkzeug.exceptions import HTTPException

from moto.core.mime_types import APP_XML


class ServiceException(Exception):
    """
    The base class for all (serializable) Moto service exceptions.

    Attributes:
        code (str): AWS service error code, e.g. ``InvalidParameterCombination``.
        message (str): A descriptive error message.

    The ``code`` and ``message`` attributes can be set as class attributes or
    provided at initialization, or a combination thereof.

    A single argument will explicitly set the message attribute:
    >>> raise ServiceException("A specific error has occurred.")

    Both class attributes overridden at initialization:
    >>> raise ServiceException("ErrorCode", "Error message")

    Notes:
       * The ``code`` value should match an exception ShapeID in the AWS model
         specification for a given service.  When the exception is serialized as
         part of a Moto server response, additional metadata from the model will
         be included (e.g. an HTTP status code).
       * If the AWS error model expects specific attributes in addition to ``message``,
         they can be set directly on the ``ServiceException`` (or subclass) object as
         class or instance attributes.
    """

    code = "UnspecifiedErrorCode"
    message = "An unspecified service error occurred"

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        if len(args) == 1:
            msg = args[0]
        elif len(args) == 2:
            self.code = args[0]
            msg = args[1]
        else:
            msg = self.message.format(**kwargs)
        Exception.__init__(self, msg)
        self.message = msg

    def __str__(self) -> str:
        return f"{self.code}: {self.message}"


SINGLE_ERROR_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<Error>
    <Code>{{error_type}}</Code>
    <Message><![CDATA[{{message}}]]></Message>
    {% block extra %}{% endblock %}
    <{{request_id_tag}}>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</{{request_id_tag}}>
</Error>
"""

WRAPPED_SINGLE_ERROR_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<ErrorResponse{% if xmlns is defined %} xmlns="{{xmlns}}"{% endif %}>
    <Error>
        <Code>{{error_type}}</Code>
        <Message><![CDATA[{{message}}]]></Message>
        {% if include_type_sender %}
        <Type>Sender</Type>
        {% endif %}
        {% block extra %}{% endblock %}
        <{{request_id_tag}}>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</{{request_id_tag}}>
    </Error>
</ErrorResponse>"""

ERROR_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
  <ErrorResponse>
    <Errors>
      <Error>
        <Code>{{error_type}}</Code>
        <Message><![CDATA[{{message}}]]></Message>
        {% block extra %}{% endblock %}
      </Error>
    </Errors>
  <{{request_id_tag}}>7a62c49f-347e-4fc4-9331-6e8eEXAMPLE</{{request_id_tag}}>
</ErrorResponse>
"""


class RESTError(HTTPException):
    code = 400
    # most APIs use <RequestId>, but some APIs (including EC2, S3) use <RequestID>
    request_id_tag_name = "RequestId"

    # When this field is set, the `Type` field will be included in the response
    # This indicates that the fault lies with the client
    include_type_sender = True

    templates = {
        "single_error": SINGLE_ERROR_RESPONSE,
        "wrapped_single_error": WRAPPED_SINGLE_ERROR_RESPONSE,
        "error": ERROR_RESPONSE,
    }
    env = Environment(loader=DictLoader(templates))

    def __init__(
        self, error_type: str, message: str, template: str = "error", **kwargs: Any
    ):
        super().__init__()
        self.error_type = error_type
        self.message = message

        if template in self.env.list_templates():
            self.description: str = self.__class__.env.get_template(template).render(
                error_type=error_type,
                message=message,
                request_id_tag=self.request_id_tag_name,
                include_type_sender=self.include_type_sender,
                **kwargs,
            )
            self.content_type = APP_XML

    def get_headers(self, *args: Any, **kwargs: Any) -> list[tuple[str, str]]:
        return [
            ("X-Amzn-ErrorType", self.relative_error_type or "UnknownError"),
            ("Content-Type", self.content_type),
        ]

    @property
    def relative_error_type(self) -> str:
        return self.error_type

    def get_body(self, *args: Any, **kwargs: Any) -> str:
        return self.description

    def to_json(self) -> "JsonRESTError":
        err = JsonRESTError(error_type=self.error_type, message=self.message)
        err.code = self.code
        return err

    @classmethod
    def extended_environment(cls, extended_templates: dict[str, str]) -> Environment:
        templates = cls.templates | extended_templates
        return Environment(loader=DictLoader(templates))


class JsonRESTError(RESTError):
    def __init__(
        self, error_type: str, message: str, template: str = "error_json", **kwargs: Any
    ):
        super().__init__(error_type, message, template, **kwargs)
        self.description: str = json.dumps(
            {"__type": self.error_type, "message": self.message}
        )
        self.content_type = "application/json"

    @property
    def relative_error_type(self) -> str:
        # https://smithy.io/2.0/aws/protocols/aws-json-1_1-protocol.html
        # If a # character is present, then take only the contents after the first # character in the value
        return (self.error_type.split("#")[-1]) if self.error_type else ""

    def get_body(self, *args: Any, **kwargs: Any) -> str:
        return self.description


class SignatureDoesNotMatchError(RESTError):
    code = 403

    def __init__(self) -> None:
        super().__init__(
            "SignatureDoesNotMatch",
            "The request signature we calculated does not match the signature you provided. Check your AWS Secret Access Key and signing method. Consult the service documentation for details.",
        )


class InvalidClientTokenIdError(RESTError):
    code = 403

    def __init__(self) -> None:
        super().__init__(
            "InvalidClientTokenId",
            "The security token included in the request is invalid.",
        )


class AccessDeniedError(RESTError):
    code = 403

    def __init__(self, user_arn: str, action: str):
        super().__init__(
            "AccessDenied", f"User: {user_arn} is not authorized to perform: {action}"
        )


class AuthFailureError(RESTError):
    code = 401

    def __init__(self) -> None:
        super().__init__(
            "AuthFailure",
            "AWS was not able to validate the provided access credentials",
        )


class AWSError(JsonRESTError):
    TYPE: Optional[str] = None
    STATUS = 400

    def __init__(
        self,
        message: str,
        exception_type: Optional[str] = None,
        status: Optional[int] = None,
    ):
        super().__init__(exception_type or self.TYPE, message)  # type: ignore[arg-type]
        self.code = status or self.STATUS


class InvalidNextTokenException(JsonRESTError):
    """For AWS Config resource listing. This will be used by many different resource types, and so it is in moto.core."""

    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidNextTokenException", "The nextToken provided is invalid"
        )


class InvalidToken(AWSError):
    code = 400

    def __init__(self, message: str = "Invalid token"):
        super().__init__(f"Invalid Token: {message}", "InvalidToken")


class ServiceNotWhitelisted(Exception):
    def __init__(self, service_name: str):
        from moto.settings import default_user_config

        services_whitelisted = default_user_config.get("core", {}).get(
            "service_whitelist"
        )
        super().__init__(
            f"Service {service_name} not whitelisted. Only services {services_whitelisted} are allowed."
        )
