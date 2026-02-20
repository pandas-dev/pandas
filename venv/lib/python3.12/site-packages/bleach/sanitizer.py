from itertools import chain
import re
import warnings

from xml.sax.saxutils import unescape

from bleach import html5lib_shim
from bleach import parse_shim


#: Set of allowed tags
ALLOWED_TAGS = frozenset(
    (
        "a",
        "abbr",
        "acronym",
        "b",
        "blockquote",
        "code",
        "em",
        "i",
        "li",
        "ol",
        "strong",
        "ul",
    )
)


#: Map of allowed attributes by tag
ALLOWED_ATTRIBUTES = {
    "a": ["href", "title"],
    "abbr": ["title"],
    "acronym": ["title"],
}

#: Set of allowed protocols
ALLOWED_PROTOCOLS = frozenset(("http", "https", "mailto"))

#: Invisible characters--0 to and including 31 except 9 (tab), 10 (lf), and 13 (cr)
INVISIBLE_CHARACTERS = "".join(
    [chr(c) for c in chain(range(0, 9), range(11, 13), range(14, 32))]
)

#: Regexp for characters that are invisible
INVISIBLE_CHARACTERS_RE = re.compile("[" + INVISIBLE_CHARACTERS + "]", re.UNICODE)

#: String to replace invisible characters with. This can be a character, a
#: string, or even a function that takes a Python re matchobj
INVISIBLE_REPLACEMENT_CHAR = "?"


class NoCssSanitizerWarning(UserWarning):
    pass


class Cleaner:
    """Cleaner for cleaning HTML fragments of malicious content

    This cleaner is a security-focused function whose sole purpose is to remove
    malicious content from a string such that it can be displayed as content in
    a web page.

    To use::

        from bleach.sanitizer import Cleaner

        cleaner = Cleaner()

        for text in all_the_yucky_things:
            sanitized = cleaner.clean(text)

    .. Note::

       This cleaner is not designed to use to transform content to be used in
       non-web-page contexts.

    .. Warning::

       This cleaner is not thread-safe--the html parser has internal state.
       Create a separate cleaner per thread!


    """

    def __init__(
        self,
        tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRIBUTES,
        protocols=ALLOWED_PROTOCOLS,
        strip=False,
        strip_comments=True,
        filters=None,
        css_sanitizer=None,
    ):
        """Initializes a Cleaner

        :arg set tags: set of allowed tags; defaults to
            ``bleach.sanitizer.ALLOWED_TAGS``

        :arg dict attributes: allowed attributes; can be a callable, list or dict;
            defaults to ``bleach.sanitizer.ALLOWED_ATTRIBUTES``

        :arg set protocols: set of allowed protocols for links; defaults
            to ``bleach.sanitizer.ALLOWED_PROTOCOLS``

        :arg bool strip: whether or not to strip disallowed elements

        :arg bool strip_comments: whether or not to strip HTML comments

        :arg list filters: list of html5lib Filter classes to pass streamed content through

            .. seealso:: http://html5lib.readthedocs.io/en/latest/movingparts.html#filters

            .. Warning::

               Using filters changes the output of ``bleach.Cleaner.clean``.
               Make sure the way the filters change the output are secure.

        :arg CSSSanitizer css_sanitizer: instance with a "sanitize_css" method for
            sanitizing style attribute values and style text; defaults to None

        """
        self.tags = tags
        self.attributes = attributes
        self.protocols = protocols
        self.strip = strip
        self.strip_comments = strip_comments
        self.filters = filters or []
        self.css_sanitizer = css_sanitizer

        self.parser = html5lib_shim.BleachHTMLParser(
            tags=self.tags,
            strip=self.strip,
            consume_entities=False,
            namespaceHTMLElements=False,
        )
        self.walker = html5lib_shim.getTreeWalker("etree")
        self.serializer = html5lib_shim.BleachHTMLSerializer(
            quote_attr_values="always",
            omit_optional_tags=False,
            escape_lt_in_attrs=True,
            # We want to leave entities as they are without escaping or
            # resolving or expanding
            resolve_entities=False,
            # Bleach has its own sanitizer, so don't use the html5lib one
            sanitize=False,
            # clean preserves attr order
            alphabetical_attributes=False,
        )

        if css_sanitizer is None:
            # FIXME(willkg): this doesn't handle when attributes or an
            # attributes value is a callable
            attributes_values = []
            if isinstance(attributes, list):
                attributes_values = attributes

            elif isinstance(attributes, dict):
                attributes_values = []
                for values in attributes.values():
                    if isinstance(values, (list, tuple)):
                        attributes_values.extend(values)

            if "style" in attributes_values:
                warnings.warn(
                    "'style' attribute specified, but css_sanitizer not set.",
                    category=NoCssSanitizerWarning,
                )

    def clean(self, text):
        """Cleans text and returns sanitized result as unicode

        :arg str text: text to be cleaned

        :returns: sanitized text as unicode

        :raises TypeError: if ``text`` is not a text type

        """
        if not isinstance(text, str):
            message = (
                f"argument cannot be of {text.__class__.__name__!r} type, "
                + "must be of text type"
            )
            raise TypeError(message)

        if not text:
            return ""

        dom = self.parser.parseFragment(text)
        filtered = BleachSanitizerFilter(
            source=self.walker(dom),
            allowed_tags=self.tags,
            attributes=self.attributes,
            strip_disallowed_tags=self.strip,
            strip_html_comments=self.strip_comments,
            css_sanitizer=self.css_sanitizer,
            allowed_protocols=self.protocols,
        )

        # Apply any filters after the BleachSanitizerFilter
        for filter_class in self.filters:
            filtered = filter_class(source=filtered)

        return self.serializer.render(filtered)


def attribute_filter_factory(attributes):
    """Generates attribute filter function for the given attributes value

    The attributes value can take one of several shapes. This returns a filter
    function appropriate to the attributes value. One nice thing about this is
    that there's less if/then shenanigans in the ``allow_token`` method.

    """
    if callable(attributes):
        return attributes

    if isinstance(attributes, dict):

        def _attr_filter(tag, attr, value):
            if tag in attributes:
                attr_val = attributes[tag]
                if callable(attr_val):
                    return attr_val(tag, attr, value)

                if attr in attr_val:
                    return True

            if "*" in attributes:
                attr_val = attributes["*"]
                if callable(attr_val):
                    return attr_val(tag, attr, value)

                return attr in attr_val

            return False

        return _attr_filter

    if isinstance(attributes, list):

        def _attr_filter(tag, attr, value):
            return attr in attributes

        return _attr_filter

    raise ValueError("attributes needs to be a callable, a list or a dict")


class BleachSanitizerFilter(html5lib_shim.SanitizerFilter):
    """html5lib Filter that sanitizes text

    This filter can be used anywhere html5lib filters can be used.

    """

    def __init__(
        self,
        source,
        allowed_tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRIBUTES,
        allowed_protocols=ALLOWED_PROTOCOLS,
        attr_val_is_uri=html5lib_shim.attr_val_is_uri,
        svg_attr_val_allows_ref=html5lib_shim.svg_attr_val_allows_ref,
        svg_allow_local_href=html5lib_shim.svg_allow_local_href,
        strip_disallowed_tags=False,
        strip_html_comments=True,
        css_sanitizer=None,
    ):
        """Creates a BleachSanitizerFilter instance

        :arg source: html5lib TreeWalker stream as an html5lib TreeWalker

        :arg set allowed_tags: set of allowed tags; defaults to
            ``bleach.sanitizer.ALLOWED_TAGS``

        :arg dict attributes: allowed attributes; can be a callable, list or dict;
            defaults to ``bleach.sanitizer.ALLOWED_ATTRIBUTES``

        :arg set allowed_protocols: set of allowed protocols for links; defaults
            to ``bleach.sanitizer.ALLOWED_PROTOCOLS``

        :arg attr_val_is_uri: set of attributes that have URI values

        :arg svg_attr_val_allows_ref: set of SVG attributes that can have
            references

        :arg svg_allow_local_href: set of SVG elements that can have local
            hrefs

        :arg bool strip_disallowed_tags: whether or not to strip disallowed
            tags

        :arg bool strip_html_comments: whether or not to strip HTML comments

        :arg CSSSanitizer css_sanitizer: instance with a "sanitize_css" method for
            sanitizing style attribute values and style text; defaults to None

        """
        # NOTE(willkg): This is the superclass of
        # html5lib.filters.sanitizer.Filter. We call this directly skipping the
        # __init__ for html5lib.filters.sanitizer.Filter because that does
        # things we don't need to do and kicks up the deprecation warning for
        # using Sanitizer.
        html5lib_shim.Filter.__init__(self, source)

        self.allowed_tags = frozenset(allowed_tags)
        self.allowed_protocols = frozenset(allowed_protocols)

        self.attr_filter = attribute_filter_factory(attributes)
        self.strip_disallowed_tags = strip_disallowed_tags
        self.strip_html_comments = strip_html_comments

        self.attr_val_is_uri = attr_val_is_uri
        self.svg_attr_val_allows_ref = svg_attr_val_allows_ref
        self.css_sanitizer = css_sanitizer
        self.svg_allow_local_href = svg_allow_local_href

    def sanitize_stream(self, token_iterator):
        for token in token_iterator:
            ret = self.sanitize_token(token)

            if not ret:
                continue

            if isinstance(ret, list):
                yield from ret
            else:
                yield ret

    def merge_characters(self, token_iterator):
        """Merge consecutive Characters tokens in a stream"""
        characters_buffer = []

        for token in token_iterator:
            if characters_buffer:
                if token["type"] == "Characters":
                    characters_buffer.append(token)
                    continue
                else:
                    # Merge all the characters tokens together into one and then
                    # operate on it.
                    new_token = {
                        "data": "".join(
                            [char_token["data"] for char_token in characters_buffer]
                        ),
                        "type": "Characters",
                    }
                    characters_buffer = []
                    yield new_token

            elif token["type"] == "Characters":
                characters_buffer.append(token)
                continue

            yield token

        new_token = {
            "data": "".join([char_token["data"] for char_token in characters_buffer]),
            "type": "Characters",
        }
        yield new_token

    def __iter__(self):
        return self.merge_characters(
            self.sanitize_stream(html5lib_shim.Filter.__iter__(self))
        )

    def sanitize_token(self, token):
        """Sanitize a token either by HTML-encoding or dropping.

        Unlike sanitizer.Filter, allowed_attributes can be a dict of {'tag':
        ['attribute', 'pairs'], 'tag': callable}.

        Here callable is a function with two arguments of attribute name and
        value. It should return true of false.

        Also gives the option to strip tags instead of encoding.

        :arg dict token: token to sanitize

        :returns: token or list of tokens

        """
        token_type = token["type"]
        if token_type in ["StartTag", "EndTag", "EmptyTag"]:
            if token["name"] in self.allowed_tags:
                return self.allow_token(token)

            elif self.strip_disallowed_tags:
                return None

            else:
                return self.disallowed_token(token)

        elif token_type == "Comment":
            if not self.strip_html_comments:
                # call lxml.sax.saxutils to escape &, <, and > in addition to " and '
                token["data"] = html5lib_shim.escape(
                    token["data"], entities={'"': "&quot;", "'": "&#x27;"}
                )
                return token
            else:
                return None

        elif token_type == "Characters":
            return self.sanitize_characters(token)

        else:
            return token

    def sanitize_characters(self, token):
        """Handles Characters tokens

        Our overridden tokenizer doesn't do anything with entities. However,
        that means that the serializer will convert all ``&`` in Characters
        tokens to ``&amp;``.

        Since we don't want that, we extract entities here and convert them to
        Entity tokens so the serializer will let them be.

        :arg token: the Characters token to work on

        :returns: a list of tokens

        """
        data = token.get("data", "")

        if not data:
            return token

        data = INVISIBLE_CHARACTERS_RE.sub(INVISIBLE_REPLACEMENT_CHAR, data)
        token["data"] = data

        # If there isn't a & in the data, we can return now
        if "&" not in data:
            return token

        new_tokens = []

        # For each possible entity that starts with a "&", we try to extract an
        # actual entity and re-tokenize accordingly
        for part in html5lib_shim.next_possible_entity(data):
            if not part:
                continue

            if part.startswith("&"):
                entity = html5lib_shim.match_entity(part)
                if entity is not None:
                    if entity == "amp":
                        # LinkifyFilter can't match urls across token boundaries
                        # which is problematic with &amp; since that shows up in
                        # querystrings all the time. This special-cases &amp;
                        # and converts it to a & and sticks it in as a
                        # Characters token. It'll get merged with surrounding
                        # tokens in the BleachSanitizerfilter.__iter__ and
                        # escaped in the serializer.
                        new_tokens.append({"type": "Characters", "data": "&"})
                    else:
                        new_tokens.append({"type": "Entity", "name": entity})

                    # Length of the entity plus 2--one for & at the beginning
                    # and one for ; at the end
                    remainder = part[len(entity) + 2 :]
                    if remainder:
                        new_tokens.append({"type": "Characters", "data": remainder})
                    continue

            new_tokens.append({"type": "Characters", "data": part})

        return new_tokens

    def sanitize_uri_value(self, value, allowed_protocols):
        """Checks a uri value to see if it's allowed

        :arg value: the uri value to sanitize
        :arg allowed_protocols: set of allowed protocols

        :returns: allowed value or None

        """
        # NOTE(willkg): This transforms the value into a normalized one that's
        # easier to match and verify, but shouldn't get returned since it's
        # vastly different than the original value.

        # Convert all character entities in the value
        normalized_uri = html5lib_shim.convert_entities(value)

        # Nix backtick, space characters, and control characters
        normalized_uri = re.sub(r"[`\000-\040\177-\240\s]+", "", normalized_uri)

        # Remove REPLACEMENT characters
        normalized_uri = normalized_uri.replace("\ufffd", "")

        # Lowercase it--this breaks the value, but makes it easier to match
        # against
        normalized_uri = normalized_uri.lower()

        try:
            # Drop attributes with uri values that have protocols that aren't
            # allowed
            parsed = parse_shim.urlparse(normalized_uri)
        except ValueError:
            # URI is impossible to parse, therefore it's not allowed
            return None

        if parsed.scheme:
            # If urlparse found a scheme, check that
            if parsed.scheme in allowed_protocols:
                return value

        else:
            # Allow uris that are just an anchor
            if normalized_uri.startswith("#"):
                return value

            # Handle protocols that urlparse doesn't recognize like "myprotocol"
            if (
                ":" in normalized_uri
                and normalized_uri.split(":")[0] in allowed_protocols
            ):
                return value

            # If there's no protocol/scheme specified, then assume it's "http" or
            # "https" and see if that's allowed
            if "http" in allowed_protocols or "https" in allowed_protocols:
                return value

        return None

    def allow_token(self, token):
        """Handles the case where we're allowing the tag"""
        if "data" in token:
            # Loop through all the attributes and drop the ones that are not
            # allowed, are unsafe or break other rules. Additionally, fix
            # attribute values that need fixing.
            #
            # At the end of this loop, we have the final set of attributes
            # we're keeping.
            attrs = {}
            for namespaced_name, val in token["data"].items():
                namespace, name = namespaced_name

                # Drop attributes that are not explicitly allowed
                #
                # NOTE(willkg): We pass in the attribute name--not a namespaced
                # name.
                if not self.attr_filter(token["name"], name, val):
                    continue

                # Drop attributes with uri values that use a disallowed protocol
                # Sanitize attributes with uri values
                if namespaced_name in self.attr_val_is_uri:
                    new_value = self.sanitize_uri_value(val, self.allowed_protocols)
                    if new_value is None:
                        continue
                    val = new_value

                # Drop values in svg attrs with non-local IRIs
                if namespaced_name in self.svg_attr_val_allows_ref:
                    new_val = re.sub(r"url\s*\(\s*[^#\s][^)]+?\)", " ", unescape(val))
                    new_val = new_val.strip()
                    if not new_val:
                        continue

                    else:
                        # Replace the val with the unescaped version because
                        # it's a iri
                        val = new_val

                # Drop href and xlink:href attr for svg elements with non-local IRIs
                if (None, token["name"]) in self.svg_allow_local_href:
                    if namespaced_name in [
                        (None, "href"),
                        (html5lib_shim.namespaces["xlink"], "href"),
                    ]:
                        if re.search(r"^\s*[^#\s]", val):
                            continue

                # If it's a style attribute, sanitize it
                if namespaced_name == (None, "style"):
                    if self.css_sanitizer:
                        val = self.css_sanitizer.sanitize_css(val)
                    else:
                        # FIXME(willkg): if style is allowed, but no
                        # css_sanitizer was set up, then this is probably a
                        # mistake and we should raise an error here
                        #
                        # For now, we're going to set the value to "" because
                        # there was no sanitizer set
                        val = ""

                # At this point, we want to keep the attribute, so add it in
                attrs[namespaced_name] = val

            token["data"] = attrs

        return token

    def disallowed_token(self, token):
        token_type = token["type"]
        if token_type == "EndTag":
            token["data"] = f"</{token['name']}>"

        elif token["data"]:
            assert token_type in ("StartTag", "EmptyTag")
            attrs = []
            for (ns, name), v in token["data"].items():
                # If we end up with a namespace, but no name, switch them so we
                # have a valid name to use.
                if ns and not name:
                    ns, name = name, ns

                # Figure out namespaced name if the namespace is appropriate
                # and exists; if the ns isn't in prefixes, then drop it.
                if ns is None or ns not in html5lib_shim.prefixes:
                    namespaced_name = name
                else:
                    namespaced_name = f"{html5lib_shim.prefixes[ns]}:{name}"

                # NOTE(willkg): HTMLSerializer escapes attribute values
                # already, so if we do it here (like HTMLSerializer does),
                # then we end up double-escaping.
                attrs.append(f' {namespaced_name}="{v}"')
            token["data"] = f"<{token['name']}{''.join(attrs)}>"

        else:
            token["data"] = f"<{token['name']}>"

        if token.get("selfClosing"):
            token["data"] = f"{token['data'][:-1]}/>"

        token["type"] = "Characters"

        del token["name"]
        return token
