"""A pandoc filter used in converting notebooks to Latex.
Converts links between notebooks to Latex cross-references.
"""

import re

from pandocfilters import RawInline, applyJSONFilters, stringify  # type:ignore[import-untyped]


def resolve_references(source):
    """
    This applies the resolve_one_reference to the text passed in via the source argument.

    This expects content in the form of a string encoded JSON object as represented
    internally in ``pandoc``.
    """
    return applyJSONFilters([resolve_one_reference], source)


def resolve_one_reference(key, val, fmt, meta):
    """
    This takes a tuple of arguments that are compatible with ``pandocfilters.walk()`` that
    allows identifying hyperlinks in the document and transforms them into valid LaTeX
    \\hyperref{} calls so that linking to headers between cells is possible.

    See the documentation in ``pandocfilters.walk()`` for further information on the meaning
    and specification of ``key``, ``val``, ``fmt``, and ``meta``.
    """

    if key == "Link":
        text = stringify(val[1])
        target = val[2][0]
        m = re.match(r"#(.+)$", target)
        if m:
            # pandoc automatically makes labels for headings.
            label = m.group(1).lower()
            label = re.sub(r"[^\w-]+", "", label)  # Strip HTML entities
            text = re.sub(r"_", r"\_", text)  # Escape underscores in display text
            return RawInline("tex", rf"\hyperref[{label}]{{{text}}}")
    return None
    # Other elements will be returned unchanged.
