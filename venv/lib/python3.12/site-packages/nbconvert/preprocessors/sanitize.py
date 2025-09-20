"""
NBConvert Preprocessor for sanitizing HTML rendering of notebooks.
"""

import warnings

from bleach import ALLOWED_ATTRIBUTES, ALLOWED_TAGS, clean
from traitlets import Any, Bool, List, Set, Unicode

from .base import Preprocessor

_USE_BLEACH_CSS_SANITIZER = False
_USE_BLEACH_STYLES = False


try:
    # bleach[css] >=5.0
    from bleach.css_sanitizer import ALLOWED_CSS_PROPERTIES as ALLOWED_STYLES
    from bleach.css_sanitizer import CSSSanitizer

    _USE_BLEACH_CSS_SANITIZER = True
    _USE_BLEACH_STYLES = False
except ImportError:
    try:
        # bleach <5
        from bleach import ALLOWED_STYLES  # type:ignore[attr-defined, no-redef]

        _USE_BLEACH_CSS_SANITIZER = False
        _USE_BLEACH_STYLES = True
        warnings.warn(
            "Support for bleach <5 will be removed in a future version of nbconvert",
            DeprecationWarning,
            stacklevel=2,
        )

    except ImportError:
        warnings.warn(
            "The installed bleach/tinycss2 do not provide CSS sanitization, "
            "please upgrade to bleach >=5",
            UserWarning,
            stacklevel=2,
        )


__all__ = ["SanitizeHTML"]


class SanitizeHTML(Preprocessor):
    """A preprocessor to sanitize html."""

    # Bleach config.
    attributes = Any(
        config=True,
        default_value=ALLOWED_ATTRIBUTES,
        help="Allowed HTML tag attributes",
    )
    tags = List(
        Unicode(),
        config=True,
        default_value=ALLOWED_TAGS,  # type:ignore[arg-type]
        help="List of HTML tags to allow",
    )
    styles = List(
        Unicode(),
        config=True,
        default_value=ALLOWED_STYLES,  # type:ignore[arg-type]
        help="Allowed CSS styles if <style> tag is allowed",
    )
    strip = Bool(
        config=True,
        default_value=False,
        help="If True, remove unsafe markup entirely instead of escaping",
    )
    strip_comments = Bool(
        config=True,
        default_value=True,
        help="If True, strip comments from escaped HTML",
    )

    # Display data config.
    safe_output_keys = Set(
        config=True,
        default_value={
            "metadata",  # Not a mimetype per-se, but expected and safe.
            "text/plain",
            "text/latex",
            "application/json",
            "image/png",
            "image/jpeg",
        },
        help="Cell output mimetypes to render without modification",
    )
    sanitized_output_types = Set(
        config=True,
        default_value={
            "text/html",
            "text/markdown",
        },
        help="Cell output types to display after escaping with Bleach.",
    )

    def preprocess_cell(self, cell, resources, cell_index):
        """
        Sanitize potentially-dangerous contents of the cell.

        Cell Types:
          raw:
            Sanitize literal HTML
          markdown:
            Sanitize literal HTML
          code:
            Sanitize outputs that could result in code execution
        """
        if cell.cell_type == "raw":
            # Sanitize all raw cells anyway.
            # Only ones with the text/html mimetype should be emitted
            # but erring on the side of safety maybe.
            cell.source = self.sanitize_html_tags(cell.source)
            return cell, resources
        if cell.cell_type == "markdown":
            cell.source = self.sanitize_html_tags(cell.source)
            return cell, resources
        if cell.cell_type == "code":
            cell.outputs = self.sanitize_code_outputs(cell.outputs)
            return cell, resources
        return None

    def sanitize_code_outputs(self, outputs):
        """
        Sanitize code cell outputs.

        Removes 'text/javascript' fields from display_data outputs, and
        runs `sanitize_html_tags` over 'text/html'.
        """
        for output in outputs:
            # These are always ascii, so nothing to escape.
            if output["output_type"] in ("stream", "error"):
                continue
            data = output.data
            to_remove = []
            for key in data:
                if key in self.safe_output_keys:
                    continue
                if key in self.sanitized_output_types:
                    self.log.info("Sanitizing %s", key)
                    data[key] = self.sanitize_html_tags(data[key])
                else:
                    # Mark key for removal. (Python doesn't allow deletion of
                    # keys from a dict during iteration)
                    to_remove.append(key)
            for key in to_remove:
                self.log.info("Removing %s", key)
                del data[key]
        return outputs

    def sanitize_html_tags(self, html_str):
        """
        Sanitize a string containing raw HTML tags.
        """
        kwargs = {
            "tags": self.tags,
            "attributes": self.attributes,
            "strip": self.strip,
            "strip_comments": self.strip_comments,
        }

        if _USE_BLEACH_CSS_SANITIZER:
            css_sanitizer = CSSSanitizer(allowed_css_properties=self.styles)
            kwargs.update(css_sanitizer=css_sanitizer)
        elif _USE_BLEACH_STYLES:
            kwargs.update(styles=self.styles)

        return clean(html_str, **kwargs)


def _get_default_css_sanitizer():
    if _USE_BLEACH_CSS_SANITIZER:
        return CSSSanitizer(allowed_css_properties=ALLOWED_STYLES)
    return None
