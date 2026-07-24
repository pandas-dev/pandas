import xml.etree.ElementTree as ET
from _typeshed import Unused
from collections.abc import Callable
from typing import Any, Literal, overload

def current_global_nsmap() -> dict[str, str]: ...

class IncrementalTree(ET.ElementTree):
    def write(  # type: ignore[override]
        self,
        file_or_filename: ET._FileWrite,
        encoding: str | None = None,
        xml_declaration: bool | None = None,
        default_namespace: str | None = None,
        method: Literal["xml", "html", "text"] | None = None,  # does not accept 'c14n', unlike parent method
        *,
        short_empty_elements: bool = True,
        nsmap: dict[str, str] | None = None,
        root_ns_only: bool = False,
        minimal_ns_only: bool = False,
    ) -> None: ...

def process_attribs(
    elem: ET.Element[Any],
    is_nsmap_scope_changed: bool | None,
    default_ns_attr_prefix: str | None,
    nsmap_scope: dict[str, str],
    global_nsmap: dict[str, str],
    new_namespace_prefixes: set[str],
    uri_to_prefix: dict[str, str],
) -> tuple[list[tuple[str, str]], str | None, dict[str, str]]: ...
def write_elem_start(
    write: Callable[..., None],
    elem: ET.Element[Any],
    nsmap_scope: dict[str, str],
    global_nsmap: dict[str, str],
    short_empty_elements: bool | None,
    is_html: bool | None,
    is_root: bool = False,
    uri_to_prefix: dict[str, str] | None = None,
    default_ns_attr_prefix: str | None = None,
    new_nsmap: dict[str, str] | None = None,
    **kwargs: Unused,
) -> tuple[str | None, dict[str, str], str | None, dict[str, str] | None, bool]: ...
@overload
def tostring(
    element: ET.Element[Any],
    encoding: None = None,
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = False,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> bytes: ...
@overload
def tostring(
    element: ET.Element[Any],
    encoding: Literal["unicode"],
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = False,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> str: ...
@overload
def tostring(
    element: ET.Element[Any],
    encoding: str,
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = False,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> Any: ...
@overload
def tostringlist(
    element: ET.Element[Any],
    encoding: None = None,
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = False,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> list[bytes]: ...
@overload
def tostringlist(
    element: ET.Element[Any],
    encoding: Literal["unicode"],
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = False,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> list[str]: ...
@overload
def tostringlist(
    element: ET.Element[Any],
    encoding: str,
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = False,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> list[Any]: ...
@overload
def compat_tostring(
    element: ET.Element[Any],
    encoding: None = None,
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = True,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> bytes: ...
@overload
def compat_tostring(
    element: ET.Element[Any],
    encoding: Literal["unicode"],
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = True,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> str: ...
@overload
def compat_tostring(
    element: ET.Element[Any],
    encoding: str,
    method: Literal["xml", "html", "text"] | None = None,
    *,
    xml_declaration: bool | None = None,
    default_namespace: str | None = None,
    short_empty_elements: bool = True,
    nsmap: dict[str, str] | None = None,
    root_ns_only: bool = True,
    minimal_ns_only: bool = False,
    tree_cls: type[ET.ElementTree] = ...,
) -> Any: ...
