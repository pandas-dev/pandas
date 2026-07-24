from lark import Lark, ParseTree, exceptions

from pathlib import Path

from rfc3987_syntax.utils import load_grammar

RFC3987_SYNTAX_PARSER_TYPE: str = "earley"
RFC3987_SYNTAX_GRAMMAR_PATH: Path = Path(__file__).parent / "syntax_rfc3987.lark"
RFC3987_SYNTAX_TERMS: list[str] = [
    "iri",
    "iri_reference",
    "absolute_iri",
    "scheme",
    "irelative_ref",
    "irelative_part"
    "ihier_part",
    "iauthority",
    "iuserinfo",
    "ihost",
    "ireg_name",
    "ipath_abempty",
    "isegment",
    "isegment_nz",
    "isegment_nz_nc",
    "ipchar",
    "iquery",
    "ifragment",
    "iunreserved",
    "ucschar",
    "iprivate",
    "sub_delims",
    "ip_literal",
    "ipvfuture",
    "ipv6address",
    "h16",
    "ls32",
    "ipv4address",
    "dec_octet",
    "digit",
    "non_zero",
    "unreserved",
    "alpha",
    "hexdig",
    "port",
    "pct_encoded",
]

grammar: str = load_grammar(RFC3987_SYNTAX_GRAMMAR_PATH)

syntax_parser = Lark(grammar, start=["iri", "iri_reference", "absolute_iri"], parser=RFC3987_SYNTAX_PARSER_TYPE)


def parse(term: str, value: str) -> ParseTree:
    return syntax_parser.parse(value, start=term)


def is_valid_syntax(term: str, value: str):
    try:
        parse(term=term, value=value)
        return True
    except exceptions.LarkError:
        return False


def make_syntax_validator(rule_name):
    parser = Lark(grammar, start=rule_name, parser=RFC3987_SYNTAX_PARSER_TYPE)

    def syntax_validator(text):
        try:
            parser.parse(text)
            return True
        except exceptions.LarkError:
            return False

    return syntax_validator


is_valid_syntax_iri = make_syntax_validator("iri")

is_valid_syntax_iri_reference = make_syntax_validator("iri_reference")

is_valid_syntax_absolute_iri = make_syntax_validator("absolute_iri")

is_valid_syntax_irelative_ref = make_syntax_validator("irelative_ref")

is_valid_syntax_irelative_part = make_syntax_validator("irelative_part")

is_valid_syntax_ihier_part = make_syntax_validator("ihier_part")

is_valid_syntax_iauthority = make_syntax_validator("iauthority")

is_valid_syntax_iuserinfo = make_syntax_validator("iuserinfo")

is_valid_syntax_ihost = make_syntax_validator("ihost")

is_valid_syntax_ireg_name = make_syntax_validator("ireg_name")

is_valid_syntax_ipath = make_syntax_validator("ipath")

is_valid_syntax_ipath_abempty = make_syntax_validator("ipath_abempty")

is_valid_syntax_ipath_absolute = make_syntax_validator("ipath_absolute")

is_valid_syntax_ipath_noscheme = make_syntax_validator("ipath_noscheme")

is_valid_syntax_ipath_rootless = make_syntax_validator("ipath_rootless")

is_valid_syntax_ipath_empty = make_syntax_validator("ipath_empty")

is_valid_syntax_isegment = make_syntax_validator("isegment")

is_valid_syntax_isegment_nz = make_syntax_validator("isegment_nz")

is_valid_syntax_isegment_nz_nc = make_syntax_validator("isegment_nz_nc")

is_valid_syntax_ipchar = make_syntax_validator("ipchar")

is_valid_syntax_iquery = make_syntax_validator("iquery")

is_valid_syntax_ifragment = make_syntax_validator("ifragment")

is_valid_syntax_iunreserved = make_syntax_validator("iunreserved")

is_valid_syntax_ucschar = make_syntax_validator("ucschar")

is_valid_syntax_iprivate = make_syntax_validator("iprivate")

is_valid_syntax_sub_delims = make_syntax_validator("sub_delims")

is_valid_syntax_ip_literal = make_syntax_validator("ip_literal")

is_valid_syntax_ipvfuture = make_syntax_validator("ipvfuture")

is_valid_syntax_ipv6address = make_syntax_validator("ipv6address")

is_valid_syntax_h16 = make_syntax_validator("h16")

is_valid_syntax_ls32 = make_syntax_validator("ls32")

is_valid_syntax_ipv4address = make_syntax_validator("ipv4address")

is_valid_syntax_dec_octet = make_syntax_validator("dec_octet")

is_valid_syntax_unreserved = make_syntax_validator("unreserved")

is_valid_syntax_alpha = make_syntax_validator("alpha")

is_valid_syntax_digit = make_syntax_validator("digit")

is_valid_syntax_hexdig = make_syntax_validator("hexdig")

is_valid_syntax_port = make_syntax_validator("port")
