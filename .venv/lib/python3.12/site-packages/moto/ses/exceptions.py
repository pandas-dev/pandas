from moto.core.exceptions import JsonRESTError, ServiceException


class SesError(ServiceException):
    pass


class MessageRejectedError(SesError):
    code = "MessageRejected"


class ConfigurationSetDoesNotExist(SesError):
    code = "ConfigurationSetDoesNotExist"


class ConfigurationSetAlreadyExists(SesError):
    code = "ConfigurationSetAlreadyExists"


class EventDestinationAlreadyExists(SesError):
    code = "EventDestinationAlreadyExists"


class TemplateNameAlreadyExists(SesError):
    code = "AlreadyExists"


class ValidationError(SesError):
    code = "ValidationError"


class InvalidParameterValue(SesError):
    code = "InvalidParameterValue"


class InvalidRenderingParameterException(SesError):
    code = "InvalidRenderingParameterException"


class TemplateDoesNotExist(SesError):
    code = "TemplateDoesNotExist"


class AlreadyExists(SesError):
    code = "AlreadyExists"


class RuleSetDoesNotExist(SesError):
    code = "RuleSetDoesNotExist"


class RuleDoesNotExist(SesError):
    code = "RuleDoesNotExist"


class CannotDelete(SesError):
    code = "CannotDelete"


class InvalidS3ConfigurationException(SesError):
    code = "InvalidS3Configuration"


class InvalidSnsTopicException(SesError):
    code = "InvalidSnsTopic"


class InvalidLambdaFunctionException(SesError):
    code = "InvalidLambdaFunction"


class MissingRenderingAttributeException(SesError):
    code = "MissingRenderingAttributeException"

    def __init__(self, var: str):
        super().__init__(f"Attribute '{var}' is not present in the rendering data.")


class NotFoundException(JsonRESTError):
    code = 404

    def __init__(self, message: str):
        super().__init__("NotFoundException", message)
