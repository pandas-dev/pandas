from _typeshed import ReadableBuffer, SupportsRead, SupportsWrite
from collections.abc import Callable, Container, Generator, Mapping
from typing import Any, overload
from typing_extensions import TypeAlias

class ParsingInterrupted(Exception): ...

# dict as attribute value is exclusive to xmlns: https://github.com/bigpick/xmltodict/commit/22541b4874365cb8d2397f23087a866b3081fd9c
_AttrValue: TypeAlias = str | dict[str, str]
_AttrDict: TypeAlias = dict[str, _AttrValue]

class _DictSAXHandler:
    path: list[tuple[str, _AttrDict | None]]
    stack: list[tuple[_AttrDict | None, list[str]]]
    data: list[str]
    item: _AttrDict | None
    item_depth: int
    xml_attribs: bool
    item_callback: Callable[[list[tuple[str, _AttrDict | None]], str | _AttrDict | None], bool]
    attr_prefix: str
    cdata_key: str
    force_cdata: bool | Container[str] | Callable[[tuple[str, _AttrDict | None], str, str], bool]
    cdata_separator: str
    postprocessor: Callable[[list[tuple[str, _AttrDict | None]], str, _AttrValue], tuple[str, _AttrValue]] | None
    dict_constructor: type
    strip_whitespace: bool
    namespace_separator: str
    namespaces: Mapping[str, str | None] | None
    namespace_declarations: dict[str, str]
    force_list: bool | Container[str] | Callable[[tuple[str, _AttrDict | None], str, str], bool] | None
    comment_key: str
    def __init__(
        self,
        item_depth: int = 0,
        item_callback: Callable[[list[tuple[str, _AttrDict | None]], str | _AttrDict | None], bool] = ...,
        xml_attribs: bool = True,
        attr_prefix: str = "@",
        cdata_key: str = "#text",
        force_cdata: bool | Container[str] | Callable[[tuple[str, _AttrDict | None], str, str], bool] = False,
        cdata_separator: str = "",
        postprocessor: Callable[[list[tuple[str, _AttrDict | None]], str, _AttrValue], tuple[str, _AttrValue]] | None = None,
        dict_constructor: type = ...,
        strip_whitespace: bool = True,
        namespace_separator: str = ":",
        namespaces: Mapping[str, str | None] | None = None,
        force_list: bool | Container[str] | Callable[[tuple[str, _AttrDict | None], str, str], bool] | None = None,
        comment_key: str = "#comment",
    ) -> None: ...
    def startNamespaceDecl(self, prefix: str, uri: str) -> None: ...
    def startElement(self, full_name: str, attrs: dict[str, str] | list[str]) -> None: ...
    def endElement(self, full_name: str) -> None: ...
    def characters(self, data: str) -> None: ...
    def comments(self, data: str) -> None: ...
    def push_data(self, item: _AttrDict | None, key: str, data: str) -> _AttrDict: ...

def parse(
    xml_input: str | ReadableBuffer | SupportsRead[bytes] | Generator[ReadableBuffer],
    encoding: str | None = None,
    expat: Any = ...,
    process_namespaces: bool = False,
    namespace_separator: str = ":",
    disable_entities: bool = True,
    process_comments: bool = False,
    *,
    item_depth: int = 0,
    item_callback: Callable[[list[tuple[str, _AttrDict | None]], str | _AttrDict | None], bool] = ...,
    xml_attribs: bool = True,
    attr_prefix: str = "@",
    cdata_key: str = "#text",
    force_cdata: bool | Container[str] | Callable[[tuple[str, _AttrDict | None], str, str], bool] = False,
    cdata_separator: str = "",
    postprocessor: Callable[[list[tuple[str, _AttrDict | None]], str, _AttrValue], tuple[str, _AttrValue]] | None = None,
    dict_constructor: type = ...,
    strip_whitespace: bool = True,
    namespaces: Mapping[str, str | None] | None = None,
    force_list: bool | Container[str] | Callable[[tuple[str, _AttrDict | None], str, str], bool] | None = None,
    comment_key: str = "#comment",
) -> dict[str, Any]: ...
@overload
def unparse(
    input_dict: Mapping[str, Any],
    output: SupportsWrite[bytes] | SupportsWrite[str],
    encoding: str = "utf-8",
    full_document: bool = True,
    short_empty_elements: bool = False,
    comment_key: str = "#comment",
    *,
    attr_prefix: str = "@",
    cdata_key: str = "#text",
    depth: int = 0,
    # preprocessor is called like (preprocessor(key, value) for key, value in input_dict.items()).
    # It is expected to return its input, or a modification thereof
    preprocessor: Callable[[str, Any], tuple[str, Any]] | None = None,
    pretty: bool = False,
    newl: str = "\n",
    indent: str | int = "\t",
    namespace_separator: str = ":",
    namespaces: Mapping[str, str | None] | None = None,
    expand_iter: str | None = None,
) -> None: ...
@overload
def unparse(
    input_dict: Mapping[str, Any],
    output: None = None,
    encoding: str = "utf-8",
    full_document: bool = True,
    short_empty_elements: bool = False,
    comment_key: str = "#comment",
    *,
    attr_prefix: str = "@",
    cdata_key: str = "#text",
    depth: int = 0,
    # preprocessor is called like (preprocessor(key, value) for key, value in input_dict.items()).
    # It is expected to return its input, or a modification thereof
    preprocessor: Callable[[str, Any], tuple[str, Any]] | None = None,
    pretty: bool = False,
    newl: str = "\n",
    indent: str | int = "\t",
    namespace_separator: str = ":",
    namespaces: Mapping[str, str | None] | None = None,
    expand_iter: str | None = None,
) -> str: ...
