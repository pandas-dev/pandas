from typing import Optional

from jinja2 import Template

from moto.core.exceptions import RESTError


class UnformattedGetAttTemplateException(Exception):
    description = (
        "Template error: resource {0} does not support attribute type {1} in Fn::GetAtt"
    )
    status_code = 400


class ValidationError(RESTError):
    def __init__(self, name_or_id: Optional[str] = None, message: Optional[str] = None):
        if message is None:
            message = f"Stack with id {name_or_id} does not exist"

        template = Template(ERROR_RESPONSE)
        super().__init__(error_type="ValidationError", message=message)
        self.description = template.render(code="ValidationError", message=message)


class MissingParameterError(RESTError):
    def __init__(self, parameter_name: str):
        template = Template(ERROR_RESPONSE)
        message = f"Missing parameter {parameter_name}"
        super().__init__(error_type="ValidationError", message=message)
        self.description = template.render(code="Missing Parameter", message=message)


class ExportNotFound(RESTError):
    """Exception to raise if a template tries to import a non-existent export"""

    def __init__(self, export_name: str):
        template = Template(ERROR_RESPONSE)
        message = f"No export named {export_name} found."
        super().__init__(error_type="ExportNotFound", message=message)
        self.description = template.render(code="ExportNotFound", message=message)


class StackSetNotEmpty(RESTError):
    def __init__(self) -> None:
        template = Template(ERROR_RESPONSE)
        message = "StackSet is not empty"
        super().__init__(error_type="StackSetNotEmptyException", message=message)
        self.description = template.render(
            code="StackSetNotEmptyException", message=message
        )


class StackSetNotFoundException(RESTError):
    def __init__(self, name: str):
        template = Template(ERROR_RESPONSE)
        message = f"StackSet {name} not found"
        super().__init__(error_type="StackSetNotFoundException", message=message)
        self.description = template.render(
            code="StackSetNotFoundException", message=message
        )


class UnsupportedAttribute(ValidationError):
    def __init__(self, resource: str, attr: str):
        template = Template(ERROR_RESPONSE)
        super().__init__()
        self.description = template.render(
            code="ValidationError",
            message=f"Template error: resource {resource} does not support attribute type {attr} in Fn::GetAtt",
        )


ERROR_RESPONSE = """<ErrorResponse xmlns="http://cloudformation.amazonaws.com/doc/2010-05-15/">
  <Error>
    <Type>Sender</Type>
    <Code>{{ code }}</Code>
    <Message>{{ message }}</Message>
  </Error>
  <RequestId>cf4c737e-5ae2-11e4-a7c9-ad44eEXAMPLE</RequestId>
</ErrorResponse>
"""
