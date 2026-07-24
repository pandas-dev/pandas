import abc
from decimal import Decimal
from typing import Any, Literal, overload
from xml.etree.ElementTree import Element

from . import base
from .styles.moduledrawers.base import QRModuleDrawer

class SvgFragmentImage(base.BaseImageWithDrawer, metaclass=abc.ABCMeta):
    kind: str
    allowed_kinds: tuple[Literal["SVG"]]
    default_drawer_class: type[QRModuleDrawer]
    unit_size: Decimal | str
    @overload
    def units(self, pixels: int | Decimal, text: Literal[False]) -> Decimal: ...
    @overload
    def units(self, pixels: int | Decimal, text: Literal[True] = True) -> str: ...
    # to_string is delegated to ET.Element.tostring, which dictates the overload
    # options here.
    @overload
    def to_string(
        self,
        *,
        encoding: None = None,
        method: str | None = None,
        xml_declaration: bool | None = None,
        default_namespace: str | None = None,
        short_empty_elements: bool = True,
    ) -> bytes: ...
    @overload
    def to_string(
        self,
        *,
        encoding: Literal["unicode"],
        method: str | None = None,
        xml_declaration: bool | None = None,
        default_namespace: str | None = None,
        short_empty_elements: bool = True,
    ) -> str: ...
    @overload
    def to_string(
        self,
        *,
        encoding: str,
        method: str | None = None,
        xml_declaration: bool | None = None,
        default_namespace: str | None = None,
        short_empty_elements: bool = True,
    ) -> Any: ...
    # the new_image method accepts arbitrary keyword arguments to accommodate
    # subclasses with additional arguments.
    def new_image(self, **kwargs: Any) -> Element: ...

class SvgImage(SvgFragmentImage, metaclass=abc.ABCMeta):
    background: str | None
    drawer_aliases: base.DrawerAliases

class SvgPathImage(SvgImage, metaclass=abc.ABCMeta):
    QR_PATH_STYLE: dict[str, str]
    path: Element | None

class SvgFillImage(SvgImage, metaclass=abc.ABCMeta):
    background: str

class SvgPathFillImage(SvgPathImage, metaclass=abc.ABCMeta):
    background: str
