# mypy: ignore-errors
from collections import OrderedDict
from collections.abc import Mapping, MutableMapping

from botocore import xform_name
from botocore.utils import parse_timestamp

UNDEFINED = object()  # Sentinel to signal the absence of a field in the input


class QueryParser:
    TIMESTAMP_FORMAT = "iso8601"

    MAP_TYPE = dict

    def __init__(self, timestamp_parser=None, map_type=None):
        if timestamp_parser is None:
            timestamp_parser = parse_timestamp
        self._timestamp_parser = timestamp_parser
        if map_type is not None:
            self.MAP_TYPE = map_type

    def parse(self, request_dict, operation_model):
        shape = operation_model.input_shape
        parsed = self._do_parse(request_dict, shape)
        return parsed if parsed is not UNDEFINED else {}

    def _do_parse(self, request_dict, shape):
        parsed = self._parse_shape(shape, request_dict["query_params"])
        return parsed if parsed is not UNDEFINED else {}

    def _parse_shape(self, shape, node, prefix=""):
        handler = getattr(self, "_handle_%s" % shape.type_name, self._default_handle)
        return handler(shape, node, prefix)

    def _gonna_recurse(self, query_params, prefix):
        if prefix == "":
            return False
        return not any([param_key.startswith(prefix) for param_key in query_params])

    def _handle_structure(self, shape, query_params, prefix=""):
        if self._gonna_recurse(query_params, prefix):
            return UNDEFINED
        parsed = self.MAP_TYPE()
        members = shape.members
        for member_name in members:
            member_shape = members[member_name]
            member_prefix = self._get_serialized_name(member_shape, member_name)
            if prefix:
                member_prefix = "%s.%s" % (prefix, member_prefix)
            value = self._parse_shape(member_shape, query_params, member_prefix)
            parsed_key = self._parsed_key_name(member_name)
            if value is not UNDEFINED:
                parsed[parsed_key] = value
        return parsed if parsed != {} else UNDEFINED

    def _handle_list(self, shape, node, prefix=""):
        # The query protocol serializes empty lists as an empty string.
        value = self._parse_shape(shape.member, node, prefix)
        if value == "":
            return []

        list_name = shape.member.serialization.get("name", "member")
        list_prefix = f"{prefix}.{list_name}"
        parsed_list = []
        i = 1
        while True:
            element_name = f"{list_prefix}.{i}"
            element_shape = shape.member
            value = self._parse_shape(element_shape, node, element_name)
            if value is UNDEFINED:
                break
            parsed_list.append(value)
            i += 1
        return parsed_list if parsed_list != [] else UNDEFINED

    def _handle_timestamp(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        return value if value is UNDEFINED else self._timestamp_parser(value)

    def _handle_boolean(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        try:
            return value.lower() == "true"
        except AttributeError:
            pass
        return UNDEFINED

    def _handle_integer(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        return value if value is UNDEFINED else int(value)

    def _handle_float(self, shape, query_params, prefix=""):
        value = self._default_handle(shape, query_params, prefix)
        return value if value is UNDEFINED else float(value)

    _handle_double = _handle_float
    _handle_long = _handle_integer

    def _default_handle(self, shape, value, prefix=""):
        default_value = shape.metadata.get("default", UNDEFINED)
        return value.get(prefix, default_value)

    def _get_serialized_name(self, shape, default_name):
        return shape.serialization.get("name", default_name)

    def _parsed_key_name(self, member_name):
        key_name = member_name
        return key_name


class XFormedDict(MutableMapping):
    """
    A Pascal/Snake case-insensitive  ``dict``-like object.

        xfd = XFormedDict()
        xfd['DBInstanceIdentifier'] = 'identifier'
        xfd['DBInstanceIdentifier'] == 'identifier'  # True
        xfd['db_instance_identifier'] == 'identifier'  # True
        list(xfd) == ['db_instance_identifier']  # True

    """

    def __init__(self, data=None, **kwargs):
        self._xform_cache = {}
        self._store = OrderedDict()
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def _xformed(self, key: str):
        return xform_name(key, _xform_cache=self._xform_cache)

    def __setitem__(self, key, value):
        # Use the xformed key for lookups, but store the actual
        # key alongside the value.
        self._store[self._xformed(key)] = (key, value)

    def __getitem__(self, key: str):
        return self._store[self._xformed(key)][1]

    def __delitem__(self, key):
        del self._store[self._xformed(key)]

    def __iter__(self):
        return (key for key in self._store.keys())

    def __len__(self):
        return len(self._store)

    def original_items(self):
        """Like iteritems(), but with all PascalCase keys."""
        return ((keyval[0], keyval[1]) for (_, keyval) in self._store.items())

    def __eq__(self, other):
        if isinstance(other, Mapping):
            other = XFormedDict(other)
        else:
            return NotImplemented
        # Compare xformed
        return dict(self.items()) == dict(other.items())

    def copy(self):
        return XFormedDict(self._store.values())

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, dict(self.items()))
