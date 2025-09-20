from strictyaml.ruamel.comments import CommentedSeq, CommentedMap
from strictyaml.compound import FixedSeq, Map
from strictyaml.validators import Validator
from strictyaml.exceptions import YAMLSerializationError
from strictyaml.scalar import Bool, EmptyDict, EmptyList, Float, Int, Str


def schema_from_document(document):
    if isinstance(document, CommentedMap):
        return Map(
            {key: schema_from_document(value) for key, value in document.items()}
        )
    elif isinstance(document, CommentedSeq):
        return FixedSeq([schema_from_document(item) for item in document])
    else:
        return Str()


def schema_from_data(data, allow_empty):
    if isinstance(data, dict):
        if len(data) == 0:
            if allow_empty:
                return EmptyDict()
            raise YAMLSerializationError(
                "Empty dicts are not serializable to StrictYAML unless schema is used."
            )
        return Map(
            {key: schema_from_data(value, allow_empty) for key, value in data.items()}
        )
    elif isinstance(data, list):
        if len(data) == 0:
            if allow_empty:
                return EmptyList()
            raise YAMLSerializationError(
                "Empty lists are not serializable to StrictYAML unless schema is used."
            )
        return FixedSeq([schema_from_data(item, allow_empty) for item in data])
    elif isinstance(data, bool):
        return Bool()
    elif isinstance(data, int):
        return Int()
    elif isinstance(data, float):
        return Float()
    else:
        return Str()


class Any(Validator):
    """
    Validates any YAML and returns simple dicts/lists of strings.
    """

    def validate(self, chunk):
        return schema_from_document(chunk.contents)(chunk)

    def to_yaml(self, data, allow_empty=False):
        """
        Args:
            allow_empty (bool): True to allow EmptyDict and EmptyList in the
                    schema generated from the data.
        """
        return schema_from_data(data, allow_empty=allow_empty).to_yaml(data)

    @property
    def key_validator(self):
        return Str()
