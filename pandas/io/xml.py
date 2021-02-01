"""
:mod:`pandas.io.xml` is a module containing functionality for dealing with
XML IO.

"""

import io
from typing import Dict, List, Optional, Union
from urllib.error import HTTPError, URLError
from warnings import warn

from pandas._typing import FilePathOrBuffer
from pandas.errors import ParserError
from pandas.util._decorators import deprecate_nonkeyword_arguments

from pandas.core.dtypes.common import is_list_like

from pandas.core.frame import DataFrame

from pandas.io.common import is_url, stringify_path, urlopen
from pandas.io.parsers import TextParser


class _EtreeFrameParser:
    """
    Internal class to parse XML into DataFrames with the Python
    standard library XML modules: `xml.etree.ElementTree`.

    Parameters
    ----------
    io : str or file-like
        This can be either a string of raw XML, a valid URL,
        file or file-like object.

    xpath : str or regex
        The XPath expression to parse required set of nodes for
        migration to `Data Frame`. `etree` supports limited XPath.

    namespacess : dict
        The namespaces defined in XML document (`xmlns:namespace='URI')
        as dicts with key being namespace and value the URI.

    elems_only : bool
        Parse only the child elements at the specified `xpath`.

    attrs_only : bool
        Parse only the attributes at the specified `xpath`.

    names : list
        Column names for Data Frame of parsed XML data.

    encoding : str
        Encoding of xml object or document.

    stylesheet : str or file-like
        URL, file, file-like object, or a raw string containing XSLT,
        `etree` does not support XSLT but retained for consistency.

    See also
    --------
    pandas.io.xml._LxmlFrameParser

    Notes
    -----
    This class serves as fall back option if user does not have
    ``lxml`` installed or user specifically requests ``etree`` parser.
    """

    from xml.etree.ElementTree import Element, ElementTree

    def __init__(
        self,
        io,
        xpath,
        namespaces,
        elems_only,
        attrs_only,
        names,
        encoding,
        stylesheet,
    ):
        self.io = io
        self.xpath = xpath
        self.namespaces = namespaces
        self.elems_only = elems_only
        self.attrs_only = attrs_only
        self.names = names
        self.encoding = encoding
        self.stylesheet = stylesheet

    def parse_data(self) -> List[Dict[str, List[str]]]:
        """
        Parse xml data.

        This method will call the other internal methods to
        validate xpath, names, parse and return specific nodes.
        """

        if self.stylesheet:
            warn(
                "To use stylesheet, you need lxml installed. "
                "Nodes will be parsed on original XML at the xpath.",
                UserWarning,
            )

        self.xml_doc = self._parse_doc()

        self._validate_path()
        self._validate_names()

        return self._parse_nodes()

    def _parse_nodes(self) -> List[Dict[str, List[str]]]:
        """
        Parse xml nodes.

        This method will parse the children and attributes of elements
        in xpath, conditionally for only elements, only attributes
        or both while optionally renaming node names.

        Raises
        ------
        ValueError
            * If only elements and only attributes are specified.

        Notes
        -----
        Namespace URIs will be removed from return node values.Also,
        elements with missing children or attributes compared to siblings
        will have optional keys filled withi None values.
        """

        elems = self.xml_doc.findall(self.xpath, namespaces=self.namespaces)

        if self.elems_only and self.attrs_only:
            raise ValueError("Either element or attributes can be parsed not both.")
        elif self.elems_only:
            if self.names:
                dicts = [
                    {
                        **(
                            {el.tag: el.text.strip()}
                            if el.text and not el.text.isspace()
                            else {}
                        ),
                        **{
                            nm: ch.text.strip() if ch.text else None
                            for nm, ch in zip(self.names, el.findall("*"))
                        },
                    }
                    for el in elems
                ]
            else:
                dicts = [
                    {
                        ch.tag: ch.text.strip() if ch.text else None
                        for ch in el.findall("*")
                    }
                    for el in elems
                ]

        elif self.attrs_only:
            dicts = [el.attrib for el in elems]

        else:
            if self.names:
                dicts = [
                    {
                        **el.attrib,
                        **(
                            {el.tag: el.text.strip()}
                            if el.text and not el.text.isspace()
                            else {}
                        ),
                        **{
                            nm: ch.text.strip() if ch.text else None
                            for nm, ch in zip(self.names, el.findall("*"))
                        },
                    }
                    for el in elems
                ]

            else:
                dicts = [
                    {
                        **el.attrib,
                        **(
                            {el.tag: el.text.strip()}
                            if el.text and not el.text.isspace()
                            else {}
                        ),
                        **{
                            ch.tag: ch.text.strip() if ch.text else None
                            for ch in el.findall("*")
                        },
                    }
                    for el in elems
                ]

        if self.namespaces:
            dicts = [
                {k.split("}")[1] if "}" in k else k: v for k, v in d.items()}
                for d in dicts
            ]

        keys = list(dict.fromkeys([k for d in dicts for k in d.keys()]))
        dicts = [{k: d[k] if k in d.keys() else None for k in keys} for d in dicts]

        if self.names:
            dicts = [
                {nm: v for nm, (k, v) in zip(self.names, d.items())} for d in dicts
            ]

        return dicts

    def _validate_path(self) -> None:
        """
        Validate xpath.

        This method checks for sytnax, evaluation, or empty nodes return.

        Raises
        ------
        SyntaxError
            * If xpah is not supported or issues with namespaces.

        Notes
        -----
        `etree` supports limited XPath. If user attempts a more complex
        expression syntax error will raise.
        """

        msg = (
            "xpath does not return any nodes. "
            "If document uses namespaces denoted with "
            "xmlns, be sure to define namespaces and "
            "use them in xpath."
        )
        try:
            elems = self.xml_doc.find(self.xpath, namespaces=self.namespaces)
            if elems is None:
                raise ValueError(msg)

            if elems is not None and elems.find("*") is None and elems.attrib is None:
                raise ValueError(msg)

        except (KeyError, SyntaxError):
            raise SyntaxError(
                "You have used an incorrect or unsupported XPath "
                "expression for etree library or you used an "
                "undeclared namespace prefix."
            )

    def _validate_names(self) -> None:
        """
        Validate names.

        This method will check if names is a list-like and aligns
        with length of parse nodes.

        Raises
        ------
        ValueError
            * If value is not a list and less then length of nodes.
        """
        if self.names:
            children = self.xml_doc.find(
                self.xpath, namespaces=self.namespaces
            ).findall("*")

            if is_list_like(self.names):
                if len(self.names) < len(children):
                    raise ValueError(
                        "names does not match length of child elements in xpath."
                    )
            else:
                raise TypeError(
                    f"{type(self.names).__name__} is not a valid type for names"
                )

    def _convert_io(self) -> Union[None, str]:
        """
        Convert io object to string.

        This method will convert io object into a string or keep
        as string, depending on object type.
        """

        obj = None

        if isinstance(self.io, str):
            obj = self.io

        if isinstance(self.io, bytes):
            obj = self.io.decode(self.encoding)

        if isinstance(self.io, io.StringIO):
            obj = self.io.getvalue()

        if isinstance(self.io, io.BytesIO):
            obj = self.io.getvalue().decode(self.encoding)

        if isinstance(self.io, io.TextIOWrapper):
            obj = self.io.read()

        if isinstance(self.io, io.BufferedReader):
            obj = self.io.read().decode(self.encoding)

        return obj

    def _parse_doc(self) -> Union[Element, ElementTree]:
        """
        Build tree from io.

        This method will parse io object into tree for parsing
        conditionally by its specific object type.

        Raises
        ------
        HttpError
            * If URL cannot be reached.

        OSError
            * If file cannot be found.

        ParseError
            * If xml document conntains syntax issues.

        ValueError
            * If io object is not readable as string or file-like object.
        """

        from xml.etree.ElementTree import ParseError, fromstring, parse

        current_doc = self._convert_io()
        if current_doc:
            is_xml = current_doc.startswith(("<?xml", "<"))
        else:
            raise ValueError("io is not a url, file, or xml string")

        is_xml = (
            (current_doc.decode(self.encoding).startswith(("<?xml", "<")))
            if isinstance(current_doc, bytes)
            else current_doc.startswith(("<?xml", "<"))
        )

        try:
            if is_url(current_doc):
                with urlopen(current_doc) as f:
                    r = parse(f)
            elif is_xml:
                r = fromstring(current_doc)
            else:
                r = parse(current_doc)
        except (URLError, HTTPError, OSError, ParseError) as e:
            raise e

        return r


class _LxmlFrameParser:
    """
    Internal class to parse XML into DataFrames with third-party
    full-featured XML library, `lxml`, that supports
    XPath 1.0 and XSLT 1.0.

    Parameters
    ----------
    io : str or file-like
        This can be either a string of raw XML, a valid URL,
        file or file-like object.

    xpath : str or regex
        The XPath expression to parse required set of nodes for
        migration to `Data Frame`.

    namespacess : dict
        The namespaces defined in XML document (`xmlns:namespace='URI')
        as dicts with key being namespace and value the URI.

    elems_only : bool
        Parse only the child elements at the specified `xpath`.

    attrs_only : bool
        Parse only the attributes at the specified `xpath`.

    names : list
        Column names for Data Frame of parsed XML data.

    encoding : str
        Encoding of xml object or document.

    stylesheet : str or file-like
        URL, file, file-like object, or a raw string containing XSLT.

    See also
    --------
    pandas.io.xml._EtreeFrameParser

    Notes
    -----
    This is the default class called with `_EtreeFrameParser` serving
    as fall back option if user does not have ``lxml`` installed.
    With `lxml`, the user enjoys the full scope of funcationality and
    efficiency.
    """

    def __init__(
        self,
        io,
        xpath,
        namespaces,
        elems_only,
        attrs_only,
        names,
        encoding,
        stylesheet,
    ):
        self.io = io
        self.xpath = xpath
        self.namespaces = namespaces
        self.elems_only = elems_only
        self.attrs_only = attrs_only
        self.names = names
        self.encoding = encoding
        self.stylesheet = stylesheet
        self.is_style = False

        self.compression = "infer"

    def parse_data(self) -> List[Dict[str, List[str]]]:
        """
        Parse xml data.

        This method will call the other internal methods to
        validate xpath, names, optionally parse and run XSLT,
        and parse original or transformed XML and return specific nodes.
        """

        self.xml_doc = self._parse_doc()

        if self.stylesheet:
            self.is_style = True
            self.xsl_doc = self._parse_doc()
            self.xml_doc = self._transform_doc()

        self._validate_path()
        self._validate_names()

        return self._parse_nodes()

    def _parse_nodes(self) -> List[Dict[str, List[str]]]:
        """
        Parse xml nodes.

        This method will parse the children and attributes of elements
        in xpath, conditionally for only elements, only attributes
        or both while optionally renaming node names.

        Raises
        ------
        ValueError
            * If only elements and only attributes are specified.

        Notes
        -----
        Namespace URIs will be removed from return node values.Also,
        elements with missing children or attributes compared to siblings
        will have optional keys filled withi None values.
        """
        elems = self.xml_doc.xpath(self.xpath, namespaces=self.namespaces)

        if self.elems_only and self.attrs_only:
            raise ValueError("Either element or attributes can be parsed not both.")

        elif self.elems_only:
            if self.names:
                dicts = [
                    {
                        **(
                            {el.tag: el.text.strip()}
                            if el.text and not el.text.isspace()
                            else {}
                        ),
                        **{
                            nm: ch.text.strip() if ch.text else None
                            for nm, ch in zip(self.names, el.xpath("*"))
                        },
                    }
                    for el in elems
                ]
            else:
                dicts = [
                    {
                        ch.tag: ch.text.strip() if ch.text else None
                        for ch in el.xpath("*")
                    }
                    for el in elems
                ]

        elif self.attrs_only:
            dicts = [el.attrib for el in elems]

        else:
            if self.names:
                dicts = [
                    {
                        **el.attrib,
                        **(
                            {el.tag: el.text.strip()}
                            if el.text and not el.text.isspace()
                            else {}
                        ),
                        **{
                            nm: ch.text.strip() if ch.text else None
                            for nm, ch in zip(self.names, el.xpath("*"))
                        },
                    }
                    for el in elems
                ]
            else:
                dicts = [
                    {
                        **el.attrib,
                        **(
                            {el.tag: el.text.strip()}
                            if el.text and not el.text.isspace()
                            else {}
                        ),
                        **{
                            ch.tag: ch.text.strip() if ch.text else None
                            for ch in el.xpath("*")
                        },
                    }
                    for el in elems
                ]

        if self.namespaces or "}" in list(dicts[0].keys())[0]:
            dicts = [
                {k.split("}")[1] if "}" in k else k: v for k, v in d.items()}
                for d in dicts
            ]

        keys = list(dict.fromkeys([k for d in dicts for k in d.keys()]))
        dicts = [{k: d[k] if k in d.keys() else None for k in keys} for d in dicts]

        if self.names:
            dicts = [
                {nm: v for nm, (k, v) in zip(self.names, d.items())} for d in dicts
            ]

        return dicts

    def _transform_doc(self):
        """
        Transform original tree using stylesheet.

        This method will transform original xml using XSLT script into
        am ideally flatter xml document for easier parsing and migration
        to Data Frame.
        """
        from lxml.etree import XSLT, XSLTApplyError, XSLTParseError

        try:
            transformer = XSLT(self.xsl_doc)
            new_doc = transformer(self.xml_doc)
        except (XSLTApplyError, XSLTParseError) as e:
            raise e

        return new_doc

    def _validate_path(self) -> None:
        """
        Validate xpath.

        This method checks for sytnax, evaluation, or empty nodes return.

        Raises
        ------
        SyntaxError
            * If xpah is not supported or issues with namespaces.

        Notes
        -----
        `etree` supports limited XPath. If user attempts a more complex
        expression syntax error will raise.
        """
        from lxml.etree import XPathEvalError, XPathSyntaxError

        try:
            elems = self.xml_doc.xpath(self.xpath, namespaces=self.namespaces)
            children = self.xml_doc.xpath(self.xpath + "/*", namespaces=self.namespaces)
            attrs = self.xml_doc.xpath(self.xpath + "/@*", namespaces=self.namespaces)

            if (elems == [] and attrs == [] and children == []) or (
                elems != [] and attrs == [] and children == []
            ):
                raise ValueError(
                    "xpath does not return any nodes. "
                    "Be sure row level nodes are in xpath. "
                    "If document uses namespaces denoted with "
                    "xmlns, be sure to define namespaces and "
                    "use them in xpath."
                )
        except (XPathEvalError, XPathSyntaxError, TypeError) as e:
            raise e

    def _validate_names(self) -> None:
        """
        Validate names.

        This method will check if names is a list and aligns with
        length of parse nodes.

        Raises
        ------
        ValueError
            * If value is not a list and less then length of nodes.
        """
        if self.names:
            children = self.xml_doc.xpath(
                self.xpath + "[1]/*", namespaces=self.namespaces
            )

            if is_list_like(self.names):
                if len(self.names) < len(children):
                    raise ValueError(
                        "names does not match length of child elements in xpath."
                    )
            else:
                raise TypeError(
                    f"{type(self.names).__name__} is not a valid type for names"
                )

    def _convert_io(self) -> Union[None, str]:
        """
        Convert filepath_or_buffer object to string.

        This method will convert io object into a string or keep
        as string, depending on object type.
        """

        obj = None

        if isinstance(self.raw_doc, str):
            obj = self.raw_doc

        if isinstance(self.raw_doc, bytes):
            obj = self.raw_doc.decode(self.encoding)

        if isinstance(self.raw_doc, io.StringIO):
            obj = self.raw_doc.getvalue()

        if isinstance(self.raw_doc, io.BytesIO):
            obj = self.raw_doc.getvalue().decode(self.encoding)

        if isinstance(self.raw_doc, io.TextIOWrapper):
            obj = self.raw_doc.read()

        if isinstance(self.raw_doc, io.BufferedReader):
            obj = self.raw_doc.read().decode(self.encoding)

        return obj

    def _parse_doc(self):
        """
        Build tree from io.

        This method will parse io object into tree for parsing
        conditionally by its specific object type.

        Raises
        ------
        HttpError
            * If URL cannot be reached.

        LookupError
            * If xml document has incorrect or unknown encoding.

        OSError
            * If file cannot be found.

        XMLSyntaxError
            * If xml document conntains syntax issues.

        ValueError
            * If io object is not readable as string or file-like object.
        """

        from lxml.etree import XML, XMLParser, XMLSyntaxError, parse

        self.raw_doc = self.stylesheet if self.is_style else self.io

        current_doc = self._convert_io()
        if current_doc:
            is_xml = current_doc.startswith(("<?xml", "<"))
        else:
            raise ValueError("io is not a url, file, or xml string")

        try:
            curr_parser = XMLParser(encoding=self.encoding)

            if is_url(current_doc):
                with urlopen(current_doc) as f:
                    r = parse(f, parser=curr_parser)
            elif is_xml:
                r = XML(bytes(current_doc, encoding=self.encoding))
            else:
                r = parse(current_doc, parser=curr_parser)
        except (LookupError, URLError, HTTPError, OSError, XMLSyntaxError) as e:
            raise e

        return r


def _data_to_frame(data, **kwargs) -> DataFrame:
    """
    Convert parsed data to Data Frame.

    This method will bind xml dictionary data of keys and values
    into named columns of Data Frame using the built-in TextParser
    class that build Data Frame and infers specific dtypes.
    """

    tags = [list(d.keys()) for d in data]
    nodes = [list(d.values()) for d in data]

    try:
        with TextParser(nodes, names=tags[0], **kwargs) as tp:
            return tp.read()
    except ParserError:
        raise ParserError(
            "XML document may be too complex for import. "
            "Try to flatten document and use distinct "
            "element and attribute names."
        )


def _parse(
    io,
    xpath,
    namespaces,
    elems_only,
    attrs_only,
    names,
    encoding,
    parser,
    stylesheet,
    **kwargs,
) -> DataFrame:
    """
    Call internal parsers.

    This method will conditionally call internal parsers:
    LxmlFrameParser and/or EtreeParser.

    Raises
    ------
    ValueError
        * If parser is not lxml or etree.e.

    Notes
    -----
    This method will raise a warning instead of module not found or
    import error if user does not have 1xml and then reverts to
    fallback option with etree parser.
    """

    if parser == "lxml":
        try:
            p = _LxmlFrameParser(
                io,
                xpath,
                namespaces,
                elems_only,
                attrs_only,
                names,
                encoding,
                stylesheet,
            )
        except ImportError:
            warn(
                "You do not have lxml installed (default parser). "
                "Instead, etree will be used.",
                ImportWarning,
            )

            p = _EtreeFrameParser(
                io,
                xpath,
                namespaces,
                elems_only,
                attrs_only,
                names,
                encoding,
                stylesheet,
            )

    elif parser == "etree":
        p = _EtreeFrameParser(
            io,
            xpath,
            namespaces,
            elems_only,
            attrs_only,
            names,
            encoding,
            stylesheet,
        )
    else:
        raise ValueError("Values for parser can only be lxml or etree.")

    data_dicts = p.parse_data()

    return _data_to_frame(data=data_dicts, **kwargs)


@deprecate_nonkeyword_arguments(version="2.0")
def read_xml(
    io: FilePathOrBuffer,
    xpath: Optional[str] = "./*",
    namespaces: Optional[Union[dict, List[dict]]] = None,
    elems_only: Optional[bool] = False,
    attrs_only: Optional[bool] = False,
    names: Optional[List[str]] = None,
    encoding: Optional[str] = "utf-8",
    parser: Optional[str] = "lxml",
    stylesheet: Optional[FilePathOrBuffer[str]] = None,
) -> DataFrame:
    r"""
    Read XML docuemnts into a ``DataFrame`` object.

    .. versionadded:: 1.3.0

    Parameters
    ----------
    io : str, path object or file-like object
        A URL, file-like object, or raw string containing XML.

    xpath : str, optional
        The XPath to parse required set of nodes for migration to DataFrame.
        XPath should return a collection of elements and not a single
        element. Note: The ``etree`` parser supports limited XPath
        expressions. For more complex XPath, use ``lxml`` which requires
        installation.

    namespaces : dict, optional
        The namespaces defined in XML document as dicts with key being
        namespace prefix and value the URI. There is no need to include all
        namespaces in XML, only the ones used in ``xpath`` expression.
        Note: if XML document uses default namespace denoted as
        `xmlns='<URI>'` without a prefix, you must assign any temporary
        namespace, like 'doc', to URI in order to parse any underlying
        nodes. For example, ::

            namespaces = {"doc": "https://example.com"}

    elems_only : bool, optional, default = False
        Parse only the child elements at the specified ``xpath``. By default,
        all child elements and non-empty text nodes are returned.

    attrs_only :  bool, optional, default = False
        Parse only the attributes at the specified ``xpath``.
        By default, all attributes are returned.

    names :  list-like, optional
        Column names for DataFrame of parsed XML data. Use this parameter to
        rename original element names and distinguish same named elements.

    encoding : str, optional, default = 'utf-8'
        Encoding of XML document.

    parser : {'lxml','etree'}, default='lxml'
        Parser module to use for retrieval of data. Only 'lxml' and
        'etree' are supported. With 'lxml' more complex XPath searches
        and ability to use XSLT stylesheet are supported. Default parser
        uses 'lxml'. If module is not installed a warning will raise and
        process will continue with 'etree'.

    stylesheet : str, path object or file-like object
        A URL, file-like object, or a raw string containing an XSLT script.
        This stylesheet should flatten complex, deeply nested XML documents.
        To use this feature you must have ``lxml`` module installed and use
        'lxml' as ``parser``. The ``xpath`` must reference nodes of
        transformed XML document generated after XSLT transformation and not
        the original XML document. Only XSLT 1.0 scripts and not later
        versions is currently supported.

    Returns
    -------
    df
        A DataFrame.

    See Also
    --------
    read_json : Convert a JSON string to pandas object.
    read_html : Read HTML tables into a list of DataFrame objects.

    Notes
    -----
    This method is best designed to import shallow XML documents in
    following format which is the ideal fit for the two-dimensions of a
    ``DataFrame`` (row by column). ::

            <root>
                <row>
                  <column1>data</column1>
                  <column2>data</column2>
                  <column3>data</column3>
                  ...
               </row>
               <row>
                  ...
               </row>
               ...
            </root>

    As a file format, XML documents can be designed any way including
    layout of elements and attributes as long as it conforms to W3C
    specifications. Therefore, this method is a convenience handler for
    a specific flatter design and not all possible XML structures.

    However, for more complex XML documents, ``stylesheet`` allows you to
    temporarily redesign original document with XSLT (a special purpose
    language) for a flatter version for migration to a DataFrame.

    This function will *always* return a single :class:`DataFrame` or raise
    exceptions due to issues with XML document, ``xpath``, or other
    parameters.

    Examples
    --------
    >>> xml = '''<?xml version='1.0' encoding='utf-8'?>
    ... <data xmlns="http://example.com">
    ...  <row>
    ...    <shape>square</shape>
    ...    <degrees>360</degrees>
    ...    <sides>4.0</sides>
    ...  </row>
    ...  <row>
    ...    <shape>circle</shape>
    ...    <degrees>360</degrees>
    ...    <sides/>
    ...  </row>
    ...  <row>
    ...    <shape>triangle</shape>
    ...    <degrees>180</degrees>
    ...    <sides>3.0</sides>
    ...  </row>
    ... </data>'''

    >>> df = pd.read_xml(xml)

    >>> df
          shape  degrees  sides
    0    square      360    4.0
    1    circle      360    NaN
    2  triangle      180    3.0

    >>> xml = '''<?xml version='1.0' encoding='utf-8'?>
    ... <data>
    ...   <row shape="square" degrees="360" sides="4.0"/>
    ...   <row shape="circle" degrees="360"/>
    ...   <row shape="triangle" degrees="180" sides="3.0"/>
    ... </data>"'''

    >>> df = pd.read_xml(xml, xpath=".//row")

    >>> df
          shape  degrees  sides
    0    square      360    4.0
    1    circle      360    NaN
    2  triangle      180    3.0

    >>> xml = '''<?xml version='1.0' encoding='utf-8'?>
    ... <doc:data xmlns:doc="https://example.com">
    ...   <doc:row>
    ...     <doc:shape>square</doc:shape>
    ...     <doc:degrees>360</doc:degrees>
    ...     <doc:sides>4.0</doc:sides>
    ...   </doc:row>
    ...   <doc:row>
    ...     <doc:shape>circle</doc:shape>
    ...     <doc:degrees>360</doc:degrees>
    ...     <doc:sides/>
    ...   </doc:row>
    ...   <doc:row>
    ...     <doc:shape>triangle</doc:shape>
    ...     <doc:degrees>180</doc:degrees>
    ...     <doc:sides>3.0</doc:sides>
    ...   </doc:row>
    ... </doc:data>'''

    >>> df = pd.read(xml,
                     xpath="//doc:row",
                     namespaces = {'doc': 'https://example.com'})

    >>> df
          shape  degrees  sides
    0    square      360    4.0
    1    circle      360    NaN
    2  triangle      180    3.0
    """

    io = stringify_path(io)

    return _parse(
        io=io,
        xpath=xpath,
        namespaces=namespaces,
        elems_only=elems_only,
        attrs_only=attrs_only,
        names=names,
        encoding=encoding,
        parser=parser,
        stylesheet=stylesheet,
    )
