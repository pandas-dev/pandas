# flake8: noqa
"""
Shim module between Bleach and html5lib. This makes it easier to upgrade the
html5lib library without having to change a lot of code.
"""

import re
import string
import warnings

# ignore html5lib deprecation warnings to use bleach; we are bleach
# apply before we import submodules that import html5lib
warnings.filterwarnings(
    "ignore",
    message="html5lib's sanitizer is deprecated",
    category=DeprecationWarning,
    module="bleach._vendor.html5lib",
)

from bleach._vendor.html5lib import (  # noqa: E402 module level import not at top of file
    HTMLParser,
    getTreeWalker,
)
from bleach._vendor.html5lib import (
    constants,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib.constants import (  # noqa: E402 module level import not at top of file
    namespaces,
    prefixes,
)
from bleach._vendor.html5lib.constants import (
    _ReparseException as ReparseException,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib.filters.base import (
    Filter,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib.filters.sanitizer import (
    allowed_protocols,
    allowed_css_properties,
    allowed_svg_properties,
    attr_val_is_uri,
    svg_attr_val_allows_ref,
    svg_allow_local_href,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib.filters.sanitizer import (
    Filter as SanitizerFilter,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib._inputstream import (
    HTMLInputStream,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib.serializer import (
    escape,
    HTMLSerializer,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib._tokenizer import (
    attributeMap,
    HTMLTokenizer,
)  # noqa: E402 module level import not at top of file
from bleach._vendor.html5lib._trie import (
    Trie,
)  # noqa: E402 module level import not at top of file


#: Map of entity name to expanded entity
ENTITIES = constants.entities

#: Trie of html entity string -> character representation
ENTITIES_TRIE = Trie(ENTITIES)

#: Token type constants--these never change
TAG_TOKEN_TYPES = {
    constants.tokenTypes["StartTag"],
    constants.tokenTypes["EndTag"],
    constants.tokenTypes["EmptyTag"],
}
TAG_TOKEN_TYPE_START = constants.tokenTypes["StartTag"]
TAG_TOKEN_TYPE_END = constants.tokenTypes["EndTag"]
TAG_TOKEN_TYPE_CHARACTERS = constants.tokenTypes["Characters"]
TAG_TOKEN_TYPE_PARSEERROR = constants.tokenTypes["ParseError"]


#: List of valid HTML tags, from WHATWG HTML Living Standard as of 2018-10-17
#: https://html.spec.whatwg.org/multipage/indices.html#elements-3
HTML_TAGS = frozenset(
    (
        "a",
        "abbr",
        "address",
        "area",
        "article",
        "aside",
        "audio",
        "b",
        "base",
        "bdi",
        "bdo",
        "blockquote",
        "body",
        "br",
        "button",
        "canvas",
        "caption",
        "cite",
        "code",
        "col",
        "colgroup",
        "data",
        "datalist",
        "dd",
        "del",
        "details",
        "dfn",
        "dialog",
        "div",
        "dl",
        "dt",
        "em",
        "embed",
        "fieldset",
        "figcaption",
        "figure",
        "footer",
        "form",
        "h1",
        "h2",
        "h3",
        "h4",
        "h5",
        "h6",
        "head",
        "header",
        "hgroup",
        "hr",
        "html",
        "i",
        "iframe",
        "img",
        "input",
        "ins",
        "kbd",
        "keygen",
        "label",
        "legend",
        "li",
        "link",
        "map",
        "mark",
        "menu",
        "meta",
        "meter",
        "nav",
        "noscript",
        "object",
        "ol",
        "optgroup",
        "option",
        "output",
        "p",
        "param",
        "picture",
        "pre",
        "progress",
        "q",
        "rp",
        "rt",
        "ruby",
        "s",
        "samp",
        "script",
        "section",
        "select",
        "slot",
        "small",
        "source",
        "span",
        "strong",
        "style",
        "sub",
        "summary",
        "sup",
        "table",
        "tbody",
        "td",
        "template",
        "textarea",
        "tfoot",
        "th",
        "thead",
        "time",
        "title",
        "tr",
        "track",
        "u",
        "ul",
        "var",
        "video",
        "wbr",
    )
)


#: List of block level HTML tags, as per https://github.com/mozilla/bleach/issues/369
#: from mozilla on 2019.07.11
#: https://developer.mozilla.org/en-US/docs/Web/HTML/Block-level_elements#Elements
HTML_TAGS_BLOCK_LEVEL = frozenset(
    (
        "address",
        "article",
        "aside",
        "blockquote",
        "details",
        "dialog",
        "dd",
        "div",
        "dl",
        "dt",
        "fieldset",
        "figcaption",
        "figure",
        "footer",
        "form",
        "h1",
        "h2",
        "h3",
        "h4",
        "h5",
        "h6",
        "header",
        "hgroup",
        "hr",
        "li",
        "main",
        "nav",
        "ol",
        "p",
        "pre",
        "section",
        "table",
        "ul",
    )
)


class InputStreamWithMemory:
    """Wraps an HTMLInputStream to remember characters since last <

    This wraps existing HTMLInputStream classes to keep track of the stream
    since the last < which marked an open tag state.

    """

    def __init__(self, inner_stream):
        self._inner_stream = inner_stream
        self.reset = self._inner_stream.reset
        self.position = self._inner_stream.position
        self._buffer = []

    @property
    def errors(self):
        return self._inner_stream.errors

    @property
    def charEncoding(self):
        return self._inner_stream.charEncoding

    @property
    def changeEncoding(self):
        return self._inner_stream.changeEncoding

    def char(self):
        c = self._inner_stream.char()
        # char() can return None if EOF, so ignore that
        if c:
            self._buffer.append(c)
        return c

    def charsUntil(self, characters, opposite=False):
        chars = self._inner_stream.charsUntil(characters, opposite=opposite)
        self._buffer.extend(list(chars))
        return chars

    def unget(self, char):
        if self._buffer:
            self._buffer.pop(-1)
        return self._inner_stream.unget(char)

    def get_tag(self):
        """Returns the stream history since last '<'

        Since the buffer starts at the last '<' as as seen by tagOpenState(),
        we know that everything from that point to when this method is called
        is the "tag" that is being tokenized.

        """
        return "".join(self._buffer)

    def start_tag(self):
        """Resets stream history to just '<'

        This gets called by tagOpenState() which marks a '<' that denotes an
        open tag. Any time we see that, we reset the buffer.

        """
        self._buffer = ["<"]


class BleachHTMLTokenizer(HTMLTokenizer):
    """Tokenizer that doesn't consume character entities"""

    def __init__(self, consume_entities=False, **kwargs):
        super().__init__(**kwargs)

        self.consume_entities = consume_entities

        # Wrap the stream with one that remembers the history
        self.stream = InputStreamWithMemory(self.stream)

        # Remember the last token emitted; needed for block element spacing
        self.emitted_last_token = None

    def __iter__(self):
        last_error_token = None

        for token in super().__iter__():
            if last_error_token is not None:
                if (
                    last_error_token["data"] == "invalid-character-in-attribute-name"
                    and token["type"] in TAG_TOKEN_TYPES
                    and token.get("data")
                ):
                    # token["data"] is an html5lib attributeMap
                    # (OrderedDict 3.7+ and dict otherwise)
                    # of attr name to attr value
                    #
                    # Remove attribute names that have ', " or < in them
                    # because those characters are invalid for attribute names.
                    token["data"] = attributeMap(
                        (attr_name, attr_value)
                        for attr_name, attr_value in token["data"].items()
                        if (
                            '"' not in attr_name
                            and "'" not in attr_name
                            and "<" not in attr_name
                        )
                    )
                    last_error_token = None
                    yield token

                elif (
                    last_error_token["data"] == "expected-closing-tag-but-got-char"
                    and self.parser.tags is not None
                    and token["data"].lower().strip() not in self.parser.tags
                ):
                    # We've got either a malformed tag or a pseudo-tag or
                    # something that html5lib wants to turn into a malformed
                    # comment which Bleach clean() will drop so we interfere
                    # with the token stream to handle it more correctly.
                    #
                    # If this is an allowed tag, it's malformed and we just let
                    # the html5lib parser deal with it--we don't enter into this
                    # block.
                    #
                    # If this is not an allowed tag, then we convert it to
                    # characters and it'll get escaped in the sanitizer.
                    token["data"] = self.stream.get_tag()
                    token["type"] = TAG_TOKEN_TYPE_CHARACTERS

                    last_error_token = None
                    yield token

                elif token["type"] == TAG_TOKEN_TYPE_PARSEERROR:
                    # If the token is a parse error, then let the last_error_token
                    # go, and make token the new last_error_token
                    yield last_error_token
                    last_error_token = token

                else:
                    yield last_error_token
                    yield token
                    last_error_token = None

                continue

            # If the token is a ParseError, we hold on to it so we can get the
            # next token and potentially fix it.
            if token["type"] == TAG_TOKEN_TYPE_PARSEERROR:
                last_error_token = token
                continue

            yield token

        if last_error_token:
            if last_error_token["data"] == "eof-in-tag-name":
                # Handle the case where the text being parsed ends with <
                # followed by a series of characters. It's treated as a tag
                # name that abruptly ends, but we should treat that like
                # character data
                yield {"type": TAG_TOKEN_TYPE_CHARACTERS, "data": self.stream.get_tag()}

            elif last_error_token["data"] in (
                "duplicate-attribute",
                "eof-in-attribute-name",
                "eof-in-attribute-value-no-quotes",
                "expected-end-of-tag-but-got-eof",
            ):
                # Handle the case where the text being parsed ends with <
                # followed by characters and then space and then:
                #
                # * more characters
                # * more characters repeated with a space between (e.g. "abc abc")
                # * more characters and then a space and then an EOF (e.g. "abc def ")
                #
                # These cases are treated as a tag name followed by an
                # attribute that abruptly ends, but we should treat that like
                # character data instead.
                yield {"type": TAG_TOKEN_TYPE_CHARACTERS, "data": self.stream.get_tag()}

            else:
                yield last_error_token

    def consumeEntity(self, allowedChar=None, fromAttribute=False):
        # If this tokenizer is set to consume entities, then we can let the
        # superclass do its thing.
        if self.consume_entities:
            return super().consumeEntity(allowedChar, fromAttribute)

        # If this tokenizer is set to not consume entities, then we don't want
        # to consume and convert them, so this overrides the html5lib tokenizer's
        # consumeEntity so that it's now a no-op.
        #
        # However, when that gets called, it's consumed an &, so we put that back in
        # the stream.
        if fromAttribute:
            self.currentToken["data"][-1][1] += "&"

        else:
            self.tokenQueue.append({"type": TAG_TOKEN_TYPE_CHARACTERS, "data": "&"})

    def tagOpenState(self):
        # This state marks a < that is either a StartTag, EndTag, EmptyTag,
        # or ParseError. In all cases, we want to drop any stream history
        # we've collected so far and we do that by calling start_tag() on
        # the input stream wrapper.
        self.stream.start_tag()
        return super().tagOpenState()

    def emitCurrentToken(self):
        token = self.currentToken

        if (
            self.parser.tags is not None
            and token["type"] in TAG_TOKEN_TYPES
            and token["name"].lower() not in self.parser.tags
        ):
            # If this is a start/end/empty tag for a tag that's not in our
            # allowed list, then it gets stripped or escaped. In both of these
            # cases it gets converted to a Characters token.
            if self.parser.strip:
                if (
                    self.emitted_last_token
                    and token["type"] == TAG_TOKEN_TYPE_START
                    and token["name"].lower() in HTML_TAGS_BLOCK_LEVEL
                ):
                    # If this is a block level tag we're stripping, we drop it
                    # for a newline because that's what a browser would parse
                    # it as
                    new_data = "\n"
                else:
                    # For all other things being stripped, we throw in an empty
                    # string token
                    new_data = ""

            else:
                # If we're escaping the token, we want to escape the exact
                # original string. Since tokenizing also normalizes data
                # and this is a tag-like thing, we've lost some information.
                # So we go back through the stream to get the original
                # string and use that.
                new_data = self.stream.get_tag()

            new_token = {"type": TAG_TOKEN_TYPE_CHARACTERS, "data": new_data}

            self.currentToken = self.emitted_last_token = new_token
            self.tokenQueue.append(new_token)
            self.state = self.dataState
            return

        self.emitted_last_token = self.currentToken
        super().emitCurrentToken()


class BleachHTMLParser(HTMLParser):
    """Parser that uses BleachHTMLTokenizer"""

    def __init__(self, tags, strip, consume_entities, **kwargs):
        """
        :arg tags: set of allowed tags--everything else is either stripped or
            escaped; if None, then this doesn't look at tags at all
        :arg strip: whether to strip disallowed tags (True) or escape them (False);
            if tags=None, then this doesn't have any effect
        :arg consume_entities: whether to consume entities (default behavior) or
            leave them as is when tokenizing (BleachHTMLTokenizer-added behavior)

        """
        self.tags = (
            frozenset((tag.lower() for tag in tags)) if tags is not None else None
        )
        self.strip = strip
        self.consume_entities = consume_entities
        super().__init__(**kwargs)

    def _parse(
        self, stream, innerHTML=False, container="div", scripting=True, **kwargs
    ):
        # set scripting=True to parse <noscript> as though JS is enabled to
        # match the expected context in browsers
        #
        # https://html.spec.whatwg.org/multipage/scripting.html#the-noscript-element
        #
        # Override HTMLParser so we can swap out the tokenizer for our own.
        self.innerHTMLMode = innerHTML
        self.container = container
        self.scripting = scripting
        self.tokenizer = BleachHTMLTokenizer(
            stream=stream, consume_entities=self.consume_entities, parser=self, **kwargs
        )
        self.reset()

        try:
            self.mainLoop()
        except ReparseException:
            self.reset()
            self.mainLoop()


def convert_entity(value):
    """Convert an entity (minus the & and ; part) into what it represents

    This handles numeric, hex, and text entities.

    :arg value: the string (minus the ``&`` and ``;`` part) to convert

    :returns: unicode character or None if it's an ambiguous ampersand that
        doesn't match a character entity

    """
    if value[0] == "#":
        if len(value) < 2:
            return None

        if value[1] in ("x", "X"):
            # hex-encoded code point
            int_as_string, base = value[2:], 16
        else:
            # decimal code point
            int_as_string, base = value[1:], 10

        if int_as_string == "":
            return None

        code_point = int(int_as_string, base)
        if 0 < code_point < 0x110000:
            return chr(code_point)
        else:
            return None

    return ENTITIES.get(value, None)


def convert_entities(text):
    """Converts all found entities in the text

    :arg text: the text to convert entities in

    :returns: unicode text with converted entities

    """
    if "&" not in text:
        return text

    new_text = []
    for part in next_possible_entity(text):
        if not part:
            continue

        if part.startswith("&"):
            entity = match_entity(part)
            if entity is not None:
                converted = convert_entity(entity)

                # If it's not an ambiguous ampersand, then replace with the
                # unicode character. Otherwise, we leave the entity in.
                if converted is not None:
                    new_text.append(converted)
                    remainder = part[len(entity) + 2 :]
                    if part:
                        new_text.append(remainder)
                    continue

        new_text.append(part)

    return "".join(new_text)


def match_entity(stream):
    """Returns first entity in stream or None if no entity exists

    Note: For Bleach purposes, entities must start with a "&" and end with a
    ";". This ignores ambiguous character entities that have no ";" at the end.

    :arg stream: the character stream

    :returns: the entity string without "&" or ";" if it's a valid character
        entity; ``None`` otherwise

    """
    # Nix the & at the beginning
    if stream[0] != "&":
        raise ValueError('Stream should begin with "&"')

    stream = stream[1:]

    stream = list(stream)
    possible_entity = ""
    end_characters = "<&=;" + string.whitespace

    # Handle number entities
    if stream and stream[0] == "#":
        possible_entity = "#"
        stream.pop(0)

        if stream and stream[0] in ("x", "X"):
            allowed = "0123456789abcdefABCDEF"
            possible_entity += stream.pop(0)
        else:
            allowed = "0123456789"

        # FIXME(willkg): Do we want to make sure these are valid number
        # entities? This doesn't do that currently.
        while stream and stream[0] not in end_characters:
            c = stream.pop(0)
            if c not in allowed:
                break
            possible_entity += c

        if possible_entity and stream and stream[0] == ";":
            return possible_entity
        return None

    # Handle character entities
    while stream and stream[0] not in end_characters:
        c = stream.pop(0)
        possible_entity += c
        if not ENTITIES_TRIE.has_keys_with_prefix(possible_entity):
            # If it's not a prefix, then it's not an entity and we're
            # out
            return None

    if possible_entity and stream and stream[0] == ";":
        return possible_entity

    return None


AMP_SPLIT_RE = re.compile("(&)")


def next_possible_entity(text):
    """Takes a text and generates a list of possible entities

    :arg text: the text to look at

    :returns: generator where each part (except the first) starts with an
        "&"

    """
    for i, part in enumerate(AMP_SPLIT_RE.split(text)):
        if i == 0:
            yield part
        elif i % 2 == 0:
            yield "&" + part


class BleachHTMLSerializer(HTMLSerializer):
    """HTMLSerializer that undoes & -> &amp; in attributes and sets
    escape_rcdata to True
    """

    # per the HTMLSerializer.__init__ docstring:
    #
    # Whether to escape characters that need to be
    # escaped within normal elements within rcdata elements such as
    # style.
    #
    escape_rcdata = True

    def escape_base_amp(self, stoken):
        """Escapes just bare & in HTML attribute values"""
        # First, undo escaping of &. We need to do this because html5lib's
        # HTMLSerializer expected the tokenizer to consume all the character
        # entities and convert them to their respective characters, but the
        # BleachHTMLTokenizer doesn't do that. For example, this fixes
        # &amp;entity; back to &entity; .
        stoken = stoken.replace("&amp;", "&")

        # However, we do want all bare & that are not marking character
        # entities to be changed to &amp;, so let's do that carefully here.
        for part in next_possible_entity(stoken):
            if not part:
                continue

            if part.startswith("&"):
                entity = match_entity(part)
                # Only leave entities in that are not ambiguous. If they're
                # ambiguous, then we escape the ampersand.
                if entity is not None and convert_entity(entity) is not None:
                    yield f"&{entity};"

                    # Length of the entity plus 2--one for & at the beginning
                    # and one for ; at the end
                    part = part[len(entity) + 2 :]
                    if part:
                        yield part
                    continue

            yield part.replace("&", "&amp;")

    def serialize(self, treewalker, encoding=None):
        """Wrap HTMLSerializer.serialize and conver & to &amp; in attribute values

        Note that this converts & to &amp; in attribute values where the & isn't
        already part of an unambiguous character entity.

        """
        in_tag = False
        after_equals = False

        for stoken in super().serialize(treewalker, encoding):
            if in_tag:
                if stoken == ">":
                    in_tag = False

                elif after_equals:
                    if stoken != '"':
                        yield from self.escape_base_amp(stoken)

                        after_equals = False
                        continue

                elif stoken == "=":
                    after_equals = True

                yield stoken
            else:
                if stoken.startswith("<"):
                    in_tag = True
                yield stoken
