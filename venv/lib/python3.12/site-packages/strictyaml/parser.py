"""
Parsing code for strictyaml.
"""
import sys

from strictyaml import ruamel as ruamelyaml
from strictyaml import exceptions
from strictyaml.ruamel.comments import CommentedSeq, CommentedMap

from strictyaml.any_validator import Any
from strictyaml.yamllocation import YAMLChunk
from strictyaml import utils

from strictyaml.ruamel.reader import Reader
from strictyaml.ruamel.scanner import RoundTripScanner
from strictyaml.ruamel.parser import RoundTripParser
from strictyaml.ruamel.composer import Composer
from strictyaml.ruamel.constructor import RoundTripConstructor
from strictyaml.ruamel.resolver import VersionedResolver
from strictyaml.ruamel.nodes import MappingNode
from strictyaml.ruamel.compat import PY2
from strictyaml.ruamel.constructor import ConstructorError

if sys.version_info[:2] > (3, 4):
    from collections.abc import Hashable
else:
    from collections import Hashable


# StrictYAMLConstructor is mostly taken from RoundTripConstructor ruamel/yaml/constructor.py
# Differences:
#  * If a duplicate key is found, an exception is raised


class StrictYAMLConstructor(RoundTripConstructor):
    yaml_constructors = {}

    def construct_mapping(self, node, maptyp, deep=False):
        if not isinstance(node, MappingNode):
            raise ConstructorError(
                None,
                None,
                "expected a mapping node, but found %s" % node.id,
                node.start_mark,
            )
        merge_map = self.flatten_mapping(node)

        # mapping = {}
        if node.comment:
            maptyp._yaml_add_comment(node.comment[:2])
            if len(node.comment) > 2:
                maptyp.yaml_end_comment_extend(node.comment[2], clear=True)
        if node.anchor:
            from strictyaml.ruamel.serializer import templated_id

            if not templated_id(node.anchor):
                maptyp.yaml_set_anchor(node.anchor)
        for key_node, value_node in node.value:
            # keys can be list -> deep
            key = self.construct_object(key_node, deep=True)
            # lists are not hashable, but tuples are
            if not isinstance(key, Hashable):
                if isinstance(key, list):
                    key = tuple(key)
            if PY2:
                try:
                    hash(key)
                except TypeError as exc:
                    raise ConstructorError(
                        "while constructing a mapping",
                        node.start_mark,
                        "found unacceptable key (%s)" % exc,
                        key_node.start_mark,
                    )
            else:
                if not isinstance(key, Hashable):
                    raise ConstructorError(
                        "while constructing a mapping",
                        node.start_mark,
                        "found unhashable key",
                        key_node.start_mark,
                    )
            value = self.construct_object(value_node, deep=deep)
            if key_node.comment:
                maptyp._yaml_add_comment(key_node.comment, key=key)
            if value_node.comment:
                maptyp._yaml_add_comment(value_node.comment, value=key)
            maptyp._yaml_set_kv_line_col(
                key,
                [
                    key_node.start_mark.line,
                    key_node.start_mark.column,
                    value_node.start_mark.line,
                    value_node.start_mark.column,
                ],
            )
            if key in maptyp:
                key_node.start_mark.name = self.label
                key_node.end_mark.name = self.label

                raise exceptions.DuplicateKeysDisallowed(
                    "While parsing",
                    key_node.start_mark,
                    "Duplicate key '{0}' found".format(key),
                    key_node.end_mark,
                )
            maptyp[key] = value
        # do this last, or <<: before a key will prevent insertion in instances
        # of collections.OrderedDict (as they have no __contains__
        if merge_map:
            maptyp.add_yaml_merge(merge_map)

        # Don't verify Mapping indentation when allowing flow,
        # as that disallows:
        #   short_key: { x = 1 }
        #   very_long_key: { x = 1 }
        if not self.allow_flow_style:
            previous_indentation = None

            for node in [
                nodegroup[1]
                for nodegroup in node.value
                if isinstance(nodegroup[1], ruamelyaml.nodes.MappingNode)
            ]:
                if previous_indentation is None:
                    previous_indentation = node.start_mark.column
                if node.start_mark.column != previous_indentation:
                    raise exceptions.InconsistentIndentationDisallowed(
                        "While parsing",
                        node.start_mark,
                        "Found mapping with indentation "
                        "inconsistent with previous mapping",
                        node.end_mark,
                    )


StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:null", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:bool", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:int", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:float", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:binary", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:timestamp", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:omap", RoundTripConstructor.construct_yaml_omap
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:pairs", RoundTripConstructor.construct_yaml_pairs
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:set", RoundTripConstructor.construct_yaml_set
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:str", RoundTripConstructor.construct_yaml_str
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:seq", RoundTripConstructor.construct_yaml_seq
)

StrictYAMLConstructor.add_constructor(
    "tag:yaml.org,2002:map", RoundTripConstructor.construct_yaml_map
)

StrictYAMLConstructor.add_constructor(None, RoundTripConstructor.construct_undefined)


# StrictYAMLScanner is mostly taken from RoundTripScanner in ruamel/yaml/scanner.py
# Differences:
#  * Tokens are checked for disallowed tokens.


class StrictYAMLScanner(RoundTripScanner):
    def check_token(self, *choices):
        # Check if the next token is one of the given types.
        while self.need_more_tokens():
            self.fetch_more_tokens()
        self._gather_comments()
        if self.tokens:
            if not choices:
                return True
            for choice in choices:
                if isinstance(self.tokens[0], choice):
                    token = self.tokens[0]
                    token.start_mark.name = self.label
                    token.end_mark.name = self.label

                    if isinstance(token, ruamelyaml.tokens.TagToken):
                        raise exceptions.TagTokenDisallowed(
                            "While scanning",
                            token.end_mark,
                            "Found disallowed tag tokens "
                            "(do not specify types in markup)",
                            token.start_mark,
                        )
                    if not self.allow_flow_style:
                        if isinstance(
                            token, ruamelyaml.tokens.FlowMappingStartToken
                        ) or isinstance(
                            token, ruamelyaml.tokens.FlowSequenceStartToken
                        ):
                            raise exceptions.FlowMappingDisallowed(
                                "While scanning",
                                token.start_mark,
                                "Found ugly disallowed JSONesque flow mapping "
                                "(surround with ' and ' to make text appear literally)",
                                token.end_mark,
                            )
                    if isinstance(token, ruamelyaml.tokens.AnchorToken):
                        raise exceptions.AnchorTokenDisallowed(
                            "While scanning",
                            token.start_mark,
                            "Found confusing disallowed anchor token "
                            "(surround with ' and ' to make text appear literally)",
                            token.end_mark,
                        )
                    return True
        return False


class StrictYAMLLoader(
    Reader,
    StrictYAMLScanner,
    RoundTripParser,
    Composer,
    StrictYAMLConstructor,
    VersionedResolver,
):
    def __init__(self, stream, version=None, preserve_quotes=None):
        Reader.__init__(self, stream, loader=self)
        StrictYAMLScanner.__init__(self, loader=self)
        RoundTripParser.__init__(self, loader=self)
        Composer.__init__(self, loader=self)
        StrictYAMLConstructor.__init__(
            self, preserve_quotes=preserve_quotes, loader=self
        )
        VersionedResolver.__init__(self, version, loader=self)


def as_document(data, schema=None, label="<unicode string>"):
    """
    Translate dicts/lists and scalar (string/bool/float/int/etc.) values into a
    YAML object which can be dumped out.
    """
    if schema is None:
        schema = Any()

    return schema(YAMLChunk(schema.to_yaml(data), label=label))


def generic_load(
    yaml_string, schema=None, label="<unicode string>", allow_flow_style=False
):
    if not utils.is_string(yaml_string):
        raise TypeError("StrictYAML can only read a string of valid YAML.")

    # We manufacture a class that has the label we want
    DynamicStrictYAMLLoader = type(
        "DynamicStrictYAMLLoader",
        (StrictYAMLLoader,),
        {"label": label, "allow_flow_style": allow_flow_style},
    )

    try:
        document = ruamelyaml.load(yaml_string, Loader=DynamicStrictYAMLLoader)
    except ruamelyaml.YAMLError as parse_error:
        if parse_error.context_mark is not None:
            parse_error.context_mark.name = label
        if parse_error.problem_mark is not None:
            parse_error.problem_mark.name = label

        raise parse_error

    # Document is just a (string, int, etc.)
    if type(document) not in (CommentedMap, CommentedSeq):
        document = yaml_string

    if schema is None:
        schema = Any()

    return schema(YAMLChunk(document, label=label))


def dirty_load(
    yaml_string, schema=None, label="<unicode string>", allow_flow_style=False
):
    """
    Parse the first YAML document in a string
    and produce corresponding YAML object.

    If allow_flow_style is set to True, then flow style is allowed.
    """
    return generic_load(
        yaml_string, schema=schema, label=label, allow_flow_style=allow_flow_style
    )


def load(yaml_string, schema=None, label="<unicode string>"):
    """
    Parse the first YAML document in a string
    and produce corresponding YAML object.
    """
    return generic_load(yaml_string, schema=schema, label=label)
