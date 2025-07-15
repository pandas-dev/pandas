from typing import TYPE_CHECKING, Any, Dict, List, Optional

from moto.core.exceptions import JsonRESTError

if TYPE_CHECKING:
    from .models.generic_type import GenericType


class SWFClientError(JsonRESTError):
    code = 400


class SWFUnknownResourceFault(SWFClientError):
    def __init__(self, resource_type: str, resource_name: Optional[str] = None):
        if resource_name:
            message = f"Unknown {resource_type}: {resource_name}"
        else:
            message = f"Unknown {resource_type}"
        super().__init__("com.amazonaws.swf.base.model#UnknownResourceFault", message)


class SWFDomainAlreadyExistsFault(SWFClientError):
    def __init__(self, domain_name: str):
        super().__init__(
            "com.amazonaws.swf.base.model#DomainAlreadyExistsFault", domain_name
        )


class SWFDomainDeprecatedFault(SWFClientError):
    def __init__(self, domain_name: str):
        super().__init__(
            "com.amazonaws.swf.base.model#DomainDeprecatedFault", domain_name
        )


class SWFSerializationException(SWFClientError):
    def __init__(self, value: Any):
        message = "class java.lang.Foo can not be converted to an String "
        message += f" (not a real SWF exception ; happened on: {value})"
        __type = "com.amazonaws.swf.base.model#SerializationException"
        super().__init__(__type, message)


class SWFTypeAlreadyExistsFault(SWFClientError):
    def __init__(self, _type: "GenericType"):
        super().__init__(
            "com.amazonaws.swf.base.model#TypeAlreadyExistsFault",
            f"{_type.__class__.__name__}=[name={_type.name}, version={_type.version}]",
        )


class SWFTypeDeprecatedFault(SWFClientError):
    def __init__(self, _type: "GenericType"):
        super().__init__(
            "com.amazonaws.swf.base.model#TypeDeprecatedFault",
            f"{_type.__class__.__name__}=[name={_type.name}, version={_type.version}]",
        )


class SWFWorkflowExecutionAlreadyStartedFault(SWFClientError):
    def __init__(self) -> None:
        super().__init__(
            "com.amazonaws.swf.base.model#WorkflowExecutionAlreadyStartedFault",
            "Already Started",
        )


class SWFDefaultUndefinedFault(SWFClientError):
    def __init__(self, key: str):
        # TODO: move that into moto.core.utils maybe?
        words = key.split("_")
        key_camel_case = words.pop(0)
        for word in words:
            key_camel_case += word.capitalize()
        super().__init__(
            "com.amazonaws.swf.base.model#DefaultUndefinedFault", key_camel_case
        )


class SWFValidationException(SWFClientError):
    def __init__(self, message: str):
        super().__init__("com.amazon.coral.validate#ValidationException", message)


class SWFDecisionValidationException(SWFClientError):
    def __init__(self, problems: List[Dict[str, Any]]):
        # messages
        messages = []
        for pb in problems:
            if pb["type"] == "null_value":
                messages.append(
                    f"Value null at '{pb['where']}' failed to satisfy constraint: Member must not be null"
                )
            elif pb["type"] == "bad_decision_type":
                messages.append(
                    f"Value '{pb['value']}' at '{pb['where']}' failed to satisfy constraint: "
                    f"Member must satisfy enum value set: [{pb['possible_values']}]"
                )
            else:
                raise ValueError(f"Unhandled decision constraint type: {pb['type']}")
        # prefix
        count = len(problems)
        if count < 2:
            prefix = "{0} validation error detected: "
        else:
            prefix = "{0} validation errors detected: "
        super().__init__(
            "com.amazon.coral.validate#ValidationException",
            prefix.format(count) + "; ".join(messages),
        )


class SWFWorkflowExecutionClosedError(Exception):
    def __str__(self) -> str:
        return repr("Cannot change this object because the WorkflowExecution is closed")
