# Class base Preprocessors
# Backwards compatibility for imported name
from nbclient.exceptions import CellExecutionError

from .base import Preprocessor
from .clearmetadata import ClearMetadataPreprocessor
from .clearoutput import ClearOutputPreprocessor
from .coalescestreams import CoalesceStreamsPreprocessor
from .convertfigures import ConvertFiguresPreprocessor
from .csshtmlheader import CSSHTMLHeaderPreprocessor
from .execute import ExecutePreprocessor
from .extractattachments import ExtractAttachmentsPreprocessor
from .extractoutput import ExtractOutputPreprocessor
from .highlightmagics import HighlightMagicsPreprocessor
from .latex import LatexPreprocessor
from .regexremove import RegexRemovePreprocessor
from .svg2pdf import SVG2PDFPreprocessor
from .tagremove import TagRemovePreprocessor

__all__ = [
    "CellExecutionError",
    "Preprocessor",
    "ClearMetadataPreprocessor",
    "ClearOutputPreprocessor",
    "CoalesceStreamsPreprocessor",
    "ConvertFiguresPreprocessor",
    "CSSHTMLHeaderPreprocessor",
    "ExecutePreprocessor",
    "ExtractAttachmentsPreprocessor",
    "ExtractOutputPreprocessor",
    "HighlightMagicsPreprocessor",
    "LatexPreprocessor",
    "RegexRemovePreprocessor",
    "SVG2PDFPreprocessor",
    "TagRemovePreprocessor",
]
