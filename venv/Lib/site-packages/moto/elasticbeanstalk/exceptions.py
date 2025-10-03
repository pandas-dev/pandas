from typing import Any

from moto.core.exceptions import RESTError

EXCEPTION_RESPONSE = """<?xml version="1.0"?>
<ErrorResponse xmlns="http://elasticache.amazonaws.com/doc/2015-02-02/">
  <Error>
    <Type>Sender</Type>
    <Code>{{ error_type }}</Code>
    <Message>{{ message }}</Message>
  </Error>
  <{{ request_id_tag }}>30c0dedb-92b1-4e2b-9be4-1188e3ed86ab</{{ request_id_tag }}>
</ErrorResponse>"""


class ElasticBeanstalkException(RESTError):
    code = 400
    extended_templates = {"ecerror": EXCEPTION_RESPONSE}
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, code: str, message: str, **kwargs: Any):
        kwargs.setdefault("template", "ecerror")
        super().__init__(code, message)


class InvalidParameterValueError(RESTError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterValue", message)


class ResourceNotFoundException(RESTError):
    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class ApplicationNotFound(ElasticBeanstalkException):
    code = 404

    def __init__(self, application_name: str):
        super().__init__(
            "ApplicationNotFound",
            message=f"Elastic Beanstalk application {application_name} not found.",
        )
