"""Module for URI Template expansion."""

from __future__ import annotations

from .expansions import ExpansionFailedError
from .uritemplate import ExpansionInvalidError, ExpansionReservedError, URITemplate
from .variable import Variable, VariableInvalidError


__all__ = (
    'URITemplate',
    'Variable',
    'ExpansionInvalidError',
    'ExpansionReservedError',
    'VariableInvalidError',
    'ExpansionFailedError',
)


def expand(template: str, **kwargs) -> (str | None):
    try:
        templ = URITemplate(template)
        return templ.expand(**kwargs)
    except Exception:
        return None


def partial(template: str, **kwargs) -> (str | None):
    try:
        templ = URITemplate(template)
        return str(templ.partial(**kwargs))
    except Exception:
        return None


def validate(template: str) -> bool:
    try:
        URITemplate(template)
        return True
    except Exception:
        return False
