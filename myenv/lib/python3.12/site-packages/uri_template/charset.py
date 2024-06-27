"""Character sets."""

from __future__ import annotations


class Charset:
    """Define character sets used in other classes."""

    ALPHA = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    DIGIT = '0123456789'
    HEX_DIGIT = '0123456789ABCDEFabcdef'
    GEN_DELIMS = ':/?#[]@'
    SUB_DELIMS = "!$&'()*+,;="
    UNRESERVED = ALPHA + DIGIT + '-._~'
    RESERVED = GEN_DELIMS + SUB_DELIMS
    VAR_START = ALPHA + DIGIT + '_'
    VAR_CHAR = VAR_START + '.'
