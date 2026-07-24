from collections.abc import Mapping, Sequence
from xml.etree.ElementTree import Element

from markdown import inlinepatterns, util
from markdown.core import Markdown
from markdown.extensions import Extension
from markdown.inlinepatterns import HtmlInlineProcessor

punctClass: str
endOfWordClass: str
closeClass: str
openingQuotesBase: str
substitutions: Mapping[str, str]
singleQuoteStartRe: str
doubleQuoteStartRe: str
doubleQuoteSetsRe: str
singleQuoteSetsRe: str
doubleQuoteSetsRe2: str
singleQuoteSetsRe2: str
decadeAbbrRe: str
openingDoubleQuotesRegex: str
closingDoubleQuotesRegex: str
closingDoubleQuotesRegex2: str
openingSingleQuotesRegex: str
closingSingleQuotesRegex: str
closingSingleQuotesRegex2: str
remainingSingleQuotesRegex: str
remainingDoubleQuotesRegex: str
HTML_STRICT_RE: str

class SubstituteTextPattern(HtmlInlineProcessor):
    replace: Sequence[int | str | Element]
    def __init__(self, pattern: str, replace: Sequence[int | str | Element], md: Markdown) -> None: ...

class SmartyExtension(Extension):
    substitutions: dict[str, str]
    def __init__(self, **kwargs) -> None: ...
    def educateDashes(self, md: Markdown) -> None: ...
    def educateEllipses(self, md: Markdown) -> None: ...
    def educateAngledQuotes(self, md: Markdown) -> None: ...
    def educateQuotes(self, md: Markdown) -> None: ...
    inlinePatterns: util.Registry[inlinepatterns.Pattern]

def makeExtension(**kwargs) -> SmartyExtension: ...
