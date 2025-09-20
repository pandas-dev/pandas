# coding: utf-8

from __future__ import absolute_import

from strictyaml.ruamel.representer import RoundTripRepresenter
from strictyaml.ruamel.scalarstring import ScalarString
from strictyaml.ruamel.emitter import Emitter
from strictyaml.ruamel.serializer import Serializer
from strictyaml.ruamel.resolver import BaseResolver
import sys

if sys.version_info[0] == 3:
    RoundTripRepresenter.add_representer(
        ScalarString, RoundTripRepresenter.represent_str
    )
else:
    RoundTripRepresenter.add_representer(
        ScalarString, RoundTripRepresenter.represent_unicode
    )


class StrictYAMLResolver(BaseResolver):
    def __init__(self, version=None, loader=None):
        BaseResolver.__init__(self, loader)


class StrictYAMLDumper(Emitter, Serializer, RoundTripRepresenter, StrictYAMLResolver):
    def __init__(
        self,
        stream,
        default_style=None,
        default_flow_style=None,
        canonical=None,
        indent=None,
        width=None,
        allow_unicode=None,
        line_break=None,
        encoding=None,
        explicit_start=None,
        explicit_end=None,
        version=None,
        tags=None,
        block_seq_indent=None,
        top_level_colon_align=None,
        prefix_colon=None,
    ):
        # type: (Any, StreamType, Any, bool, Union[None, int], Union[None, int], bool, Any, Any, Union[None, bool], Union[None, bool], Any, Any, Any, Any, Any) -> None  # NOQA
        Emitter.__init__(
            self,
            stream,
            canonical=canonical,
            indent=indent,
            width=width,
            allow_unicode=allow_unicode,
            line_break=line_break,
            block_seq_indent=block_seq_indent,
            top_level_colon_align=top_level_colon_align,
            prefix_colon=prefix_colon,
            dumper=self,
        )
        Serializer.__init__(
            self,
            encoding=encoding,
            explicit_start=explicit_start,
            explicit_end=explicit_end,
            version=version,
            tags=tags,
            dumper=self,
        )
        RoundTripRepresenter.__init__(
            self,
            default_style=default_style,
            default_flow_style=default_flow_style,
            dumper=self,
        )
        StrictYAMLResolver.__init__(self, loader=self)
