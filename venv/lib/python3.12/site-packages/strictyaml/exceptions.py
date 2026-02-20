from strictyaml.ruamel.error import MarkedYAMLError
from strictyaml.ruamel.dumper import RoundTripDumper
from strictyaml.ruamel import dump

try:
    from strictyaml.ruamel.error import Mark as StringMark
except ImportError:
    from strictyaml.ruamel.error import StringMark


class StrictYAMLError(MarkedYAMLError):
    pass


class InvalidValidatorError(StrictYAMLError):
    pass


class CannotBuildDocumentFromInvalidData(StrictYAMLError):
    pass


class CannotBuildDocumentsFromEmptyDictOrList(StrictYAMLError):
    pass


class YAMLSerializationError(StrictYAMLError):
    pass


class InvalidOptionalDefault(YAMLSerializationError):
    pass


class YAMLValidationError(StrictYAMLError):
    def __init__(self, context, problem, chunk):
        self.context = context
        self.problem = problem
        self._chunk = chunk
        self.note = None

    @property
    def context_mark(self):
        context_line = self._chunk.start_line() - 1
        str_document = dump(self._chunk.whole_document, Dumper=RoundTripDumper)
        context_index = len("\n".join(str_document.split("\n")[:context_line]))
        return StringMark(
            self._chunk.label,
            context_index,
            context_line,
            0,
            str_document,
            context_index + 1,
        )

    @property
    def problem_mark(self):
        problem_line = self._chunk.end_line() - 1
        str_document = dump(self._chunk.whole_document, Dumper=RoundTripDumper)
        problem_index = len("\n".join(str_document.split("\n")[:problem_line]))
        return StringMark(
            self._chunk.label,
            problem_index,
            problem_line,
            0,
            str_document,
            problem_index + 1,
        )


class DisallowedToken(StrictYAMLError):
    MESSAGE = "Disallowed token"


class TagTokenDisallowed(DisallowedToken):
    MESSAGE = "Tag tokens not allowed"


class FlowMappingDisallowed(DisallowedToken):
    MESSAGE = "Flow mapping tokens not allowed"


class AnchorTokenDisallowed(DisallowedToken):
    MESSAGE = "Anchor tokens not allowed"


class DuplicateKeysDisallowed(DisallowedToken):
    MESSAGE = "Duplicate keys not allowed"


class InconsistentIndentationDisallowed(DisallowedToken):
    MESSAGE = "Inconsistent indentation not allowed"


def raise_type_error(yaml_object, to_type, alternatives):
    raise TypeError(
        ("Cannot cast {0} to {1}.\n" "Use {2} instead.").format(
            repr(yaml_object), to_type, alternatives
        )
    )
