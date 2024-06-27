from typing import Any

from moto.core.exceptions import AWSError, JsonRESTError, RESTError

ERROR_WITH_MODEL_NAME = """{% extends 'single_error' %}
{% block extra %}<ModelName>{{ model }}</ModelName>{% endblock %}
"""


class SagemakerClientError(RESTError):
    extended_templates = {"model_error": ERROR_WITH_MODEL_NAME}
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "single_error")
        super().__init__(*args, **kwargs)


class ModelError(RESTError):
    extended_templates = {"model_error": ERROR_WITH_MODEL_NAME}
    env = RESTError.extended_environment(extended_templates)

    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("template", "model_error")
        super().__init__(*args, **kwargs)


class MissingModel(ModelError):
    code = 404

    def __init__(self, model: str):
        super().__init__("NoSuchModel", "Could not find model", model=model)


class ValidationError(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class AWSValidationException(AWSError):
    TYPE = "ValidationException"


class ResourceInUseException(AWSError):
    TYPE = "ResourceInUse"


class ResourceNotFound(JsonRESTError):
    def __init__(self, message: str):
        super().__init__(__class__.__name__, message)  # type: ignore
