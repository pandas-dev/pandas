from _typeshed import SupportsWrite
from collections.abc import Iterable, Iterator, Mapping
from decimal import Decimal
from pathlib import Path
from typing import Any, Final, Literal, overload
from typing_extensions import TypeAlias, deprecated

from qrcode.image.svg import SvgPathImage

# NOTE: Since svgwrite doesn't have any stubs we provide some type aliases, that
#       we can choose to refine in the future, e.g. using a Protocol
_SvgDrawing: TypeAlias = Any  # svgwrite.Drawing
_SvgGroup: TypeAlias = Any  # svgwrite.container.Group

# NOTE: Eventually we may want to consider replacing this with typed dicts, even
#       if that means disallowing non-dict arguments. It should allow anything
#       that is valid to pass into `Address.create`.
_AddressDict: TypeAlias = Mapping[str, str | None]

IBAN_ALLOWED_COUNTRIES: list[str]
QR_IID: dict[str, int]
AMOUNT_REGEX: str
MM_TO_UU: float
BILL_HEIGHT: int
RECEIPT_WIDTH: str
PAYMENT_WIDTH: str
MAX_CHARS_PAYMENT_LINE: int
MAX_CHARS_RECEIPT_LINE: int
A4: tuple[str, str]
LABELS: dict[str, dict[str, str]]
SCISSORS_SVG_PATH: str

class Address:
    @overload
    @classmethod
    def create(cls, *, name: str | None = None, line1: str, line2: str | None = None, country: str | None) -> CombinedAddress: ...
    @overload
    @classmethod
    def create(cls, *, name: str | None = None, line1: str | None = None, line2: str, country: str | None) -> CombinedAddress: ...
    @overload
    @classmethod
    def create(
        cls,
        *,
        name: str,
        street: str | None = None,
        house_num: str | None = None,
        pcode: str,
        city: str,
        country: str | None = None,
    ) -> StructuredAddress: ...
    @staticmethod
    def parse_country(country: str | None) -> str: ...

class CombinedAddress(Address):
    combined: Final = True
    name: str
    line1: str
    line2: str
    country: str
    def __init__(
        self, *, name: str | None = None, line1: str | None = None, line2: str | None = None, country: str | None = None
    ) -> None: ...
    def data_list(self) -> list[str]: ...
    def as_paragraph(self, max_chars: int = 72) -> Iterator[str]: ...

class StructuredAddress(Address):
    combined: Final = False
    name: str
    street: str
    house_num: str
    pcode: str
    city: str
    country: str
    def __init__(
        self,
        *,
        name: str | None = None,
        street: str | None = None,
        house_num: str | None = None,
        pcode: str | None = None,
        city: str | None = None,
        country: str | None = None,
    ) -> None: ...
    def data_list(self) -> list[str]: ...
    def as_paragraph(self, max_chars: int = 72) -> Iterator[str]: ...

class QRBill:
    qr_type: str
    version: str
    coding: int
    allowed_currencies: tuple[Literal["CHF"], Literal["EUR"]]
    font_family: str
    creditor: CombinedAddress | StructuredAddress
    final_creditor: CombinedAddress | StructuredAddress | None
    debtor: CombinedAddress | StructuredAddress | None
    ref_type: str
    reference_number: str | None
    account: str
    account_is_qriban: bool
    amount: str | None
    currency: Literal["CHF", "EUR"]
    additional_information: str
    billing_information: str
    @overload
    def __init__(
        self,
        account: str,
        creditor: _AddressDict,
        final_creditor: None = None,
        amount: Decimal | str | None = None,
        currency: Literal["CHF", "EUR"] = "CHF",
        debtor: _AddressDict | None = None,
        ref_number: None = None,
        reference_number: str | None = None,
        extra_infos: Literal[""] = "",
        additional_information: str = "",
        billing_information: str = "",
        alt_procs: list[str] | tuple[()] | tuple[str] | tuple[str, str] = (),
        language: Literal["en", "de", "fr", "it"] = "en",
        top_line: bool = True,
        payment_line: bool = True,
        font_factor: int = 1,
    ) -> None: ...
    @overload
    @deprecated("ref_number is deprecated and replaced by reference_number")
    def __init__(
        self,
        account: str,
        creditor: _AddressDict,
        final_creditor: None = None,
        amount: Decimal | str | None = None,
        currency: Literal["CHF", "EUR"] = "CHF",
        debtor: _AddressDict | None = None,
        *,
        ref_number: str,
        reference_number: None = None,
        extra_infos: str = "",
        additional_information: str = "",
        billing_information: str = "",
        alt_procs: list[str] | tuple[()] | tuple[str] | tuple[str, str] = (),
        language: Literal["en", "de", "fr", "it"] = "en",
        top_line: bool = True,
        payment_line: bool = True,
        font_factor: int = 1,
    ) -> None: ...
    @overload
    @deprecated("extra_infos is deprecated and replaced by additional_information")
    def __init__(
        self,
        account: str,
        creditor: _AddressDict,
        final_creditor: None = None,
        amount: Decimal | str | None = None,
        currency: Literal["CHF", "EUR"] = "CHF",
        debtor: _AddressDict | None = None,
        ref_number: None = None,
        reference_number: str | None = None,
        *,
        extra_infos: str,
        additional_information: str = "",
        billing_information: str = "",
        alt_procs: list[str] | tuple[()] | tuple[str] | tuple[str, str] = (),
        language: Literal["en", "de", "fr", "it"] = "en",
        top_line: bool = True,
        payment_line: bool = True,
        font_factor: int = 1,
    ) -> None: ...
    @property
    def title_font_info(self) -> dict[str, Any]: ...
    @property
    def font_info(self) -> dict[str, Any]: ...
    def head_font_info(self, part: str | None = None) -> dict[str, Any]: ...
    @property
    def proc_font_info(self) -> dict[str, Any]: ...
    def qr_data(self) -> str: ...
    def qr_image(self) -> SvgPathImage: ...
    def draw_swiss_cross(self, dwg: _SvgDrawing, grp: _SvgGroup, origin: tuple[float, float], size: float) -> None: ...
    def draw_blank_rect(self, dwg: _SvgDrawing, grp: _SvgGroup, x: float, y: float, width: float, height: float) -> None: ...
    def label(self, txt: str) -> str: ...
    def as_svg(self, file_out: str | Path | SupportsWrite[str], full_page: bool = False) -> None: ...
    def transform_to_full_page(self, dwg: _SvgDrawing, bill: _SvgGroup) -> None: ...
    def draw_bill(self, dwg: _SvgDrawing, horiz_scissors: bool = True) -> _SvgGroup: ...

def add_mm(*mms: str | float) -> float: ...
def mm(val: str | float) -> float: ...
def format_ref_number(bill: QRBill) -> str: ...
def format_amount(amount_: str | float) -> str: ...
def wrap_infos(infos: Iterable[str]) -> Iterator[str]: ...
def replace_linebreaks(text: str | None) -> str: ...
