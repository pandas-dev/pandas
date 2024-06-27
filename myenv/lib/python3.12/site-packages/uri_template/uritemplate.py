"""Process URI templates per http://tools.ietf.org/html/rfc6570."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from .expansions import (CommaExpansion, Expansion,
                         FormStyleQueryContinuation, FormStyleQueryExpansion,
                         FragmentExpansion, LabelExpansion, Literal,
                         PathExpansion, PathStyleExpansion,
                         ReservedCommaExpansion, ReservedExpansion, SimpleExpansion)

if (TYPE_CHECKING):
    from collections.abc import Iterable
    from .variable import Variable


class ExpansionReservedError(Exception):
    """Exception thrown for reserved but unsupported expansions."""

    expansion: str

    def __init__(self, expansion: str) -> None:
        self.expansion = expansion

    def __str__(self) -> str:
        """Convert to string."""
        return 'Unsupported expansion: ' + self.expansion


class ExpansionInvalidError(Exception):
    """Exception thrown for unknown expansions."""

    expansion: str

    def __init__(self, expansion: str) -> None:
        self.expansion = expansion

    def __str__(self) -> str:
        """Convert to string."""
        return 'Bad expansion: ' + self.expansion


class URITemplate:
    """
    URI Template object.

    Constructor may raise ExpansionReservedError, ExpansionInvalidError, or VariableInvalidError.
    """

    expansions: list[Expansion]

    def __init__(self, template: str) -> None:
        self.expansions = []
        parts = re.split(r'(\{[^\}]*\})', template)
        for part in parts:
            if (part):
                if (('{' == part[0]) and ('}' == part[-1])):
                    expansion = part[1:-1]
                    if (re.match('^([a-zA-Z0-9_]|%[0-9a-fA-F][0-9a-fA-F]).*$', expansion)):
                        self.expansions.append(SimpleExpansion(expansion))
                    elif ('+' == part[1]):
                        self.expansions.append(ReservedExpansion(expansion))
                    elif ('#' == part[1]):
                        self.expansions.append(FragmentExpansion(expansion))
                    elif ('.' == part[1]):
                        self.expansions.append(LabelExpansion(expansion))
                    elif ('/' == part[1]):
                        self.expansions.append(PathExpansion(expansion))
                    elif (';' == part[1]):
                        self.expansions.append(PathStyleExpansion(expansion))
                    elif ('?' == part[1]):
                        self.expansions.append(FormStyleQueryExpansion(expansion))
                    elif ('&' == part[1]):
                        self.expansions.append(FormStyleQueryContinuation(expansion))
                    elif (',' == part[1]):
                        if ((1 < len(part)) and ('+' == part[2])):
                            self.expansions.append(ReservedCommaExpansion(expansion))
                        else:
                            self.expansions.append(CommaExpansion(expansion))
                    elif (part[1] in '=!@|'):
                        raise ExpansionReservedError(part)
                    else:
                        raise ExpansionInvalidError(part)
                else:
                    if (('{' not in part) and ('}' not in part)):
                        self.expansions.append(Literal(part))
                    else:
                        raise ExpansionInvalidError(part)

    @property
    def variables(self) -> Iterable[Variable]:
        """Get all variables in template."""
        vars: dict[str, Variable] = {}
        for expansion in self.expansions:
            for var in expansion.variables:
                vars[var.name] = var
        return vars.values()

    @property
    def variable_names(self) -> Iterable[str]:
        """Get names of all variables in template."""
        vars: dict[str, Variable] = {}
        for expansion in self.expansions:
            for var in expansion.variables:
                vars[var.name] = var
        return [var.name for var in vars.values()]

    def expand(self, **kwargs) -> str:
        """
        Expand the template.

        May raise ExpansionFailed if a composite value is passed to a variable with a prefix modifier.
        """
        expanded = [expansion.expand(kwargs) for expansion in self.expansions]
        return ''.join([expansion for expansion in expanded if (expansion is not None)])

    def partial(self, **kwargs) -> URITemplate:
        """
        Expand the template, preserving expansions for missing variables.

        May raise ExpansionFailed if a composite value is passed to a variable with a prefix modifier.
        """
        expanded = [expansion.partial(kwargs) for expansion in self.expansions]
        return URITemplate(''.join(expanded))

    @property
    def expanded(self) -> bool:
        """Determine if template is fully expanded."""
        return (str(self) == self.expand())

    def __str__(self) -> str:
        """Convert to string, returns original template."""
        return ''.join([str(expansion) for expansion in self.expansions])

    def __repr__(self) -> str:
        """Convert to string, returns original template."""
        return str(self)
