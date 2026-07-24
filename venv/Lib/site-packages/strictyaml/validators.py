from strictyaml.exceptions import YAMLValidationError, YAMLSerializationError
from strictyaml.exceptions import InvalidValidatorError
from strictyaml.representation import YAML
from strictyaml import utils
import sys


if sys.version_info[0] == 3:
    unicode = str


class Validator(object):
    def __or__(self, other):
        return OrValidator(self, other)

    def __call__(self, chunk):
        self.validate(chunk)
        return YAML(chunk, validator=self)

    def __repr__(self):
        return "{0}()".format(self.__class__.__name__)


class MapValidator(Validator):
    def _should_be_mapping(self, data):
        if not isinstance(data, dict):
            raise YAMLSerializationError("Expected a dict, found '{}'".format(data))
        if len(data) == 0:
            raise YAMLSerializationError(
                (
                    "Expected a non-empty dict, found an empty dict.\n"
                    "Use EmptyDict validator to serialize empty dicts."
                )
            )


class SeqValidator(Validator):
    def _should_be_list(self, data):
        if not isinstance(data, list):
            raise YAMLSerializationError("Expected a list, found '{}'".format(data))
        if len(data) == 0:
            raise YAMLSerializationError(
                (
                    "Expected a non-empty list, found an empty list.\n"
                    "Use EmptyList validator to serialize empty lists."
                )
            )


class OrValidator(Validator):
    def __init__(self, validator_a, validator_b):
        assert isinstance(validator_a, Validator), "validator_a must be a Validator"
        assert isinstance(validator_b, Validator), "validator_b must be a Validator"

        self._validator_a = validator_a
        self._validator_b = validator_b

        def unpacked(validator):
            if isinstance(validator, OrValidator):
                return [
                    unpacked(validator._validator_a),
                    unpacked(validator._validator_b),
                ]
            else:
                return [validator]

        map_validator_count = len(
            [
                validator
                for validator in list(utils.flatten(unpacked(self)))
                if isinstance(validator, MapValidator)
            ]
        )

        if map_validator_count > 1:
            raise InvalidValidatorError(
                (
                    "You tried to Or ('|') together {} Map validators. "
                    "Try using revalidation instead."
                ).format(map_validator_count)
            )

        seq_validator_count = len(
            [
                validator
                for validator in list(utils.flatten(unpacked(self)))
                if isinstance(validator, SeqValidator)
            ]
        )

        if seq_validator_count > 1:
            raise InvalidValidatorError(
                (
                    "You tried to Or ('|') together {} Seq validators. "
                    "Try using revalidation instead."
                ).format(seq_validator_count)
            )

    def to_yaml(self, value):
        try:
            return self._validator_a.to_yaml(value)
        except YAMLSerializationError:
            return self._validator_b.to_yaml(value)

    def __call__(self, chunk):
        try:
            result = self._validator_a(chunk)
            result._selected_validator = result._validator
            result._validator = self
            return result
        except YAMLValidationError:
            result = self._validator_b(chunk)
            result._selected_validator = result._validator
            result._validator = self
            return result

    def __repr__(self):
        return "{0} | {1}".format(repr(self._validator_a), repr(self._validator_b))
