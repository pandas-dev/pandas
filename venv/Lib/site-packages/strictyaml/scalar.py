import math

from strictyaml.exceptions import YAMLSerializationError
from strictyaml.validators import Validator
from strictyaml.representation import YAML
from strictyaml import constants
from strictyaml import utils
from datetime import datetime
import dateutil.parser
import decimal
import sys
import re
import urllib.parse
from strictyaml.ruamel.scalarstring import PreservedScalarString


if sys.version_info[0] == 3:
    unicode = str


class ScalarValidator(Validator):
    @property
    def rule_description(self):
        return "a {0}".format(self.__class__.__name__.lower())

    def __call__(self, chunk):
        chunk.expect_scalar(self.rule_description)
        return YAML(chunk, validator=self)

    def validate(self, chunk):
        return self.validate_scalar(chunk)

    def should_be_string(self, data, message):
        if not utils.is_string(data):
            raise YAMLSerializationError(
                "{0} got '{1}' of type {2}.".format(message, data, type(data).__name__)
            )

    def validate_scalar(self, chunk):
        raise NotImplementedError("validate_scalar(self, chunk) must be implemented")


class Enum(ScalarValidator):
    def __init__(self, restricted_to, item_validator=None):
        self._item_validator = Str() if item_validator is None else item_validator
        assert isinstance(
            self._item_validator, ScalarValidator
        ), "item validator must be scalar too"
        self._restricted_to = restricted_to

    def validate_scalar(self, chunk):
        val = self._item_validator(chunk)
        val._validator = self
        if val.scalar not in self._restricted_to:
            chunk.expecting_but_found(
                "when expecting one of: {0}".format(
                    ", ".join(map(str, self._restricted_to))
                )
            )
        else:
            return val

    def to_yaml(self, data):
        if data not in self._restricted_to:
            raise YAMLSerializationError(
                "Got '{0}' when  expecting one of: {1}".format(
                    data, ", ".join(map(str, self._restricted_to))
                )
            )
        return self._item_validator.to_yaml(data)

    def __repr__(self):
        # TODO : item_validator
        return "Enum({0})".format(repr(self._restricted_to))


class CommaSeparated(ScalarValidator):
    def __init__(self, item_validator):
        self._item_validator = item_validator
        assert isinstance(
            self._item_validator, ScalarValidator
        ), "item validator must be scalar too"

    def validate_scalar(self, chunk):
        if chunk.contents == "":
            return []
        return [
            self._item_validator.validate_scalar(
                chunk.textslice(positions[0], positions[1])
            )
            for positions in utils.comma_separated_positions(chunk.contents)
        ]

    def to_yaml(self, data):
        if isinstance(data, list):
            return ", ".join([self._item_validator.to_yaml(item) for item in data])
        elif utils.is_string(data):
            for item in data.split(","):
                self._item_validator.to_yaml(item)
            return data
        else:
            raise YAMLSerializationError(
                "expected string or list, got '{}' of type '{}'".format(
                    data, type(data).__name__
                )
            )

    def __repr__(self):
        return "CommaSeparated({0})".format(self._item_validator)


class Regex(ScalarValidator):
    def __init__(self, regular_expression):
        """
        Give regular expression, e.g. u'[0-9]'
        """
        self._regex = regular_expression
        # re.fullmatch is only available in Python 3.4+ so append "$" if needed
        if not regular_expression.endswith(r"$"):
            regular_expression += r"$"
        self._fullmatch = re.compile(regular_expression).match
        self._matching_message = "when expecting string matching {0}".format(
            self._regex
        )

    def validate_scalar(self, chunk):
        if self._fullmatch(chunk.contents) is None:
            chunk.expecting_but_found(
                self._matching_message, "found non-matching string"
            )
        return chunk.contents

    def to_yaml(self, data):
        self.should_be_string(data, self._matching_message)
        if self._fullmatch(data) is None:
            raise YAMLSerializationError(
                "{} found '{}'".format(self._matching_message, data)
            )
        return data


class Email(Regex):
    def __init__(self):
        super(Email, self).__init__(constants.REGEXES["email"])
        self._matching_message = "when expecting an email address"


class Url(ScalarValidator):
    def __is_absolute_url(self, raw):
        try:
            ret = urllib.parse.urlparse(raw)
            return ret.scheme != "" and ret.netloc != ""
        except ValueError:
            return False

    def validate_scalar(self, chunk):
        if not self.__is_absolute_url(chunk.contents):
            chunk.expecting_but_found("when expecting a URL")
        return chunk.contents

    def to_yaml(self, data):
        self.should_be_string(data, "expected a URL,")
        if not self.__is_absolute_url(data):
            raise YAMLSerializationError("'{}' is not a URL".format(data))
        return data


class Str(ScalarValidator):
    def validate_scalar(self, chunk):
        return chunk.contents

    def to_yaml(self, data):
        if not utils.is_string(data):
            raise YAMLSerializationError("'{}' is not a string".format(data))
        if "\n" in data:
            return PreservedScalarString(data)
        return data


class Int(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if not utils.is_integer(val):
            chunk.expecting_but_found("when expecting an integer")
        else:
            # Only Python 3.6+ supports underscores in numeric literals
            return int(val.replace("_", ""))

    def to_yaml(self, data):
        if utils.is_string(data) or isinstance(data, int):
            if utils.is_integer(str(data)):
                return str(data)
        raise YAMLSerializationError("'{}' not an integer.".format(data))


class HexInt(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if not utils.is_hexadecimal(val):
            chunk.expecting_but_found("when expecting a hexadecimal integer")
        return int(val, 16)

    def to_yaml(self, data):
        if utils.is_hexadecimal(data):
            if isinstance(data, int):
                return hex(data)
            else:
                return data
        raise YAMLSerializationError("'{}' not a hexademial integer.".format(data))


class Bool(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if unicode(val).lower() not in constants.BOOL_VALUES:
            chunk.expecting_but_found(
                """when expecting a boolean value (one of "{0}")""".format(
                    '", "'.join(constants.BOOL_VALUES)
                )
            )
        else:
            if val.lower() in constants.TRUE_VALUES:
                return True
            else:
                return False

    def to_yaml(self, data):
        if not isinstance(data, bool):
            if str(data).lower() in constants.BOOL_VALUES:
                return data
            else:
                raise YAMLSerializationError("Not a boolean")
        else:
            return "yes" if data else "no"


class Float(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if utils.is_infinity(val) or utils.is_not_a_number(val):
            val = val.replace(".", "")
        elif not utils.is_decimal(val):
            chunk.expecting_but_found("when expecting a float")
        # Only Python 3.6+ supports underscores in numeric literals
        return float(val.replace("_", ""))

    def to_yaml(self, data):
        if utils.has_number_type(data):
            if math.isnan(data):
                return "nan"
            if data == float("inf"):
                return "inf"
            if data == float("-inf"):
                return "-inf"
            return str(data)
        if utils.is_string(data) and utils.is_decimal(data):
            return data
        raise YAMLSerializationError("when expecting a float, got '{}'".format(data))


class Decimal(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if not utils.is_decimal(val):
            chunk.expecting_but_found("when expecting a decimal")
        else:
            return decimal.Decimal(val)


class Datetime(ScalarValidator):
    def validate_scalar(self, chunk):
        try:
            return dateutil.parser.parse(chunk.contents)
        except ValueError:
            chunk.expecting_but_found("when expecting a datetime")

    def to_yaml(self, data):
        if isinstance(data, datetime):
            return data.isoformat()
        if utils.is_string(data):
            try:
                dateutil.parser.parse(data)
                return data
            except ValueError:
                raise YAMLSerializationError(
                    "expected a datetime, got '{}'".format(data)
                )
        raise YAMLSerializationError(
            "expected a datetime, got '{}' of type '{}'".format(
                data, type(data).__name__
            )
        )


class NullNone(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if val.lower() != "null":
            chunk.expecting_but_found(
                "when expecting a 'null', got '{}' instead.".format(val)
            )
        else:
            return self.empty(chunk)

    def empty(self, chunk):
        return None

    def to_yaml(self, data):
        if data is None:
            return "null"
        raise YAMLSerializationError("expected None, got '{}'")


class EmptyNone(ScalarValidator):
    def validate_scalar(self, chunk):
        val = chunk.contents
        if val != "":
            chunk.expecting_but_found("when expecting an empty value")
        else:
            return self.empty(chunk)

    def empty(self, chunk):
        return None

    def to_yaml(self, data):
        if data is None:
            return ""
        raise YAMLSerializationError("expected None, got '{}'")


class EmptyDict(EmptyNone):
    def empty(self, chunk):
        return {}

    def to_yaml(self, data):
        if data == {}:
            return ""
        raise YAMLSerializationError("Not an empty dict")


class EmptyList(EmptyNone):
    def empty(self, chunk):
        return []

    def to_yaml(self, data):
        if data == []:
            return ""
        raise YAMLSerializationError("expected empty list, got '{}'")
