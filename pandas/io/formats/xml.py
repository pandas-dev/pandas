"""
Module for formatting output data in XML.
"""

import codecs
import io
from typing import Dict, List, Optional, Union
from urllib.error import HTTPError, URLError
from warnings import warn

from pandas._typing import FilePathOrBuffer

from pandas.core.dtypes.common import is_list_like

from pandas.io.common import is_url, urlopen
from pandas.io.formats.format import DataFrameFormatter


class EtreeXMLFormatter:
    """
    Class for formatting data in xml using Python standard library
    modules: `xml.etree.ElementTree` and `xml.dom.minidom`.

    Parameters
    ----------
    io : str or file-like
        This can be either a string of raw XML, a valid URL,
        file or file-like object.

    index : bool
        Whether to include index in xml document.

    row_name : str
        Name for root of xml document. Default is 'data'.

    root_name : str
        Name for row elemens of xml document. Default is 'row'.

    na_rep : str
        Missing data representation.

    attrs_cols : list
        List of columns to write as attributes in row element.

    elem_cols : list
        List of columns to write as children in row element.

    namespacess : dict
        The namespaces to define in XML document as dicts with key
        being namespace and value the URI.

    prefix : str
        The prefix for each element in XML document including root.

    encoding : str
        Encoding of xml object or document.

    xml_declaration : bool
        Whether to include xml declaration at top line item in xml.

    pretty_print : bool
        Whether to write xml document with line breaks and indentation.

    stylesheet : str or file-like
        A URL, file, file-like object, or a raw string containing XSLT,
        `etree` does not support XSLT but retained for consistency.

    See also
    --------
    pandas.io.formats.xml.LxmlXMLFormatter

    Notes
    -----
    This class serves as fall back option if user does not have
    ``lxml`` installed or user specifically requests ``etree`` parser.
    """

    def __init__(
        self,
        formatter: DataFrameFormatter,
        io: Optional[FilePathOrBuffer[str]] = None,
        index: Optional[bool] = True,
        root_name: Optional[str] = "data",
        row_name: Optional[str] = "row",
        na_rep: Optional[str] = None,
        attr_cols: Optional[Union[str, List[str]]] = None,
        elem_cols: Optional[Union[str, List[str]]] = None,
        namespaces: Optional[Dict[str, str]] = None,
        prefix: Optional[str] = None,
        encoding: Optional[str] = "utf-8",
        xml_declaration: Optional[bool] = True,
        pretty_print: Optional[bool] = True,
        stylesheet: Optional[FilePathOrBuffer[str]] = None,
    ) -> None:
        self.fmt = formatter
        self.io = io
        self.index = index
        self.root_name = root_name
        self.row_name = row_name
        self.na_rep = na_rep
        self.attr_cols = attr_cols
        self.elem_cols = elem_cols
        self.namespaces = namespaces
        self.prefix = prefix
        self.encoding = encoding
        self.xml_declaration = xml_declaration
        self.pretty_print = pretty_print
        self.stylesheet = stylesheet
        self.frame = self.fmt.frame

        self.validate_columns()
        self.validate_encoding()
        self.orig_cols = self.fmt.frame.columns.tolist()
        self.frame_dicts = self.process_dataframe()
        self.handle_indexes()
        self.prefix_uri = self.get_prefix_uri()

    def build_tree(self) -> bytes:
        """
        Build tree from  data.

        This method initializes the root and builds attributes and elements
        with optional namespaces.
        """
        from xml.etree.ElementTree import Element, SubElement, tostring

        self.root = Element(
            f"{self.prefix_uri}{self.root_name}", attrib=self.other_namespaces()
        )

        for k, d in self.frame_dicts.items():
            self.d = d
            self.elem_row = SubElement(self.root, f"{self.prefix_uri}{self.row_name}")

            if self.attr_cols:
                self.build_attribs()
            if self.elem_cols:
                self.build_elems()
            if not self.attr_cols and not self.elem_cols:
                self.elem_cols = list(self.frame_dicts[0].keys())
                self.build_elems()

        self.out_xml = tostring(self.root, method="xml", encoding=self.encoding)

        if self.pretty_print:
            self.out_xml = self.prettify_tree()

        if not self.xml_declaration:
            self.out_xml = self.remove_declaration()

        if self.stylesheet:
            warn(
                "To use stylesheet, you need lxml installed. "
                "The non-transformed, original XML is returned instead.",
                UserWarning,
            )

        return self.out_xml

    def validate_columns(self) -> None:
        """
        Validate elems_cols and attrs_cols.

        This method will check if columns is list-like.

        Raises
        ------
        ValueError
            * If value is not a list and less then length of nodes.
        """
        if self.attr_cols and not is_list_like(self.attr_cols):
            raise TypeError(
                f"{type(self.attr_cols).__name__} is not a valid type for attr_cols"
            )

        if self.elem_cols and not is_list_like(self.elem_cols):
            raise TypeError(
                f"{type(self.elem_cols).__name__} is not a valid type for elem_cols"
            )

    def validate_encoding(self) -> None:
        """
        Validate encoding.

        This method will check if encoding is among listed under codecs.

        Raises
        ------
        LookupError
            * If encoding is not available in codecs.
        """

        try:
            codecs.lookup(self.encoding)
        except LookupError as e:
            raise e

    def process_dataframe(self) -> None:
        """
        Adjust Data Frame to fit xml output.

        This method will adjust underlying data frame for xml output,
        including replacing missing entities and including indexes.
        """

        na_dict = {"None": self.na_rep, "NaN": self.na_rep, "nan": self.na_rep}

        df = (
            (self.fmt.frame.reset_index().applymap(str).replace(na_dict))
            if self.index
            else self.fmt.frame.applymap(str).replace(na_dict)
        )

        return df.to_dict(orient="index")

    def handle_indexes(self) -> None:
        """
        Handle indexes.

        This method will add indexes into attr_cols or elem_cols.
        """

        indexes = [x for x in self.frame_dicts[0].keys() if x not in self.orig_cols]

        if self.attr_cols and self.index:
            self.attr_cols = list(indexes) + self.attr_cols

        if self.elem_cols and self.index:
            self.elem_cols = list(indexes) + self.elem_cols

    def get_prefix_uri(self) -> str:
        """
        Get uri of namespace prefix.

        This method retrieves corresponding URI to prefix in namespaces.
        """

        from xml.etree.ElementTree import register_namespace

        uri = ""
        if self.namespaces:
            for p, n in self.namespaces.items():
                register_namespace(p, n)
            if self.prefix:
                try:
                    uri = f"{{{self.namespaces[self.prefix]}}}"
                except (KeyError):
                    raise KeyError("prefix is not included in namespaces")
            else:
                uri = f'{{{self.namespaces[""]}}}'

        return uri

    def other_namespaces(self) -> dict:
        """
        Define other namespaces.

        This method will build dictionary of namespaces attributes
        for root element, conditionally with optional namespaces and
        prefix.
        """

        nmsp_dict = {}
        if self.namespaces and self.prefix is None:
            nmsp_dict = {"xmlns": n for p, n in self.namespaces.items() if p != ""}

        if self.namespaces and self.prefix:
            nmsp_dict = {"xmlns": n for p, n in self.namespaces.items() if p == ""}

        return nmsp_dict

    def build_attribs(self) -> None:
        """
        Create attributes of row.

        This method adds attributes using attr_cols to row element and
        works with tuples for multindex or hierarchical columns.
        """

        for col in self.attr_cols:
            flat_col = col
            if isinstance(col, tuple):
                flat_col = (
                    "".join(str(c) for c in col).strip()
                    if "" in col
                    else "_".join(str(c) for c in col).strip()
                )

            attr_name = f"{self.prefix_uri}{flat_col}"
            try:
                if self.d[col] is not None:
                    self.elem_row.attrib[attr_name] = str(self.d[col])
            except KeyError:
                raise KeyError(f"no valid column, {col}")

    def build_elems(self) -> None:
        """
        Create child elements of row.

        This method adds child elements using elem_cols to row element and
        works with tuples for multindex or hierarchical columns.
        """

        from xml.etree.ElementTree import SubElement

        for col in self.elem_cols:
            flat_col = col
            if isinstance(col, tuple):
                flat_col = (
                    "".join(str(c) for c in col).strip()
                    if "" in col
                    else "_".join(str(c) for c in col).strip()
                )

            elem_name = f"{self.prefix_uri}{flat_col}"
            try:
                val = None if self.d[col] in [None, ""] else str(self.d[col])
                SubElement(self.elem_row, elem_name).text = val
            except KeyError:
                raise KeyError(f"no valid column, {col}")

    def prettify_tree(self) -> bytes:
        """
        Output tree for pretty print format.

        This method will pretty print xml with line breaks and indentation.
        """

        from xml.dom.minidom import parseString

        dom = parseString(self.out_xml)

        return dom.toprettyxml(indent="  ", encoding=self.encoding)

    def remove_declaration(self) -> None:
        """
        Remove xml declaration.

        This method will remove xml declaration of working tree. Currently,
        pretty_print is not supported in etree.
        """

        return self.out_xml.split(b"?>")[-1].strip()

    def write_output(self) -> Optional[str]:
        xml_doc = self.build_tree()

        try:
            if self.io:
                with open(self.io, "wb") as f:
                    f.write(xml_doc)
                xml_doc = None
            else:
                xml_doc = xml_doc.decode(self.encoding).rstrip()
        except (UnicodeDecodeError, OSError) as e:
            raise e

        return xml_doc


class LxmlXMLFormatter:
    """
    Class for formatting data in xml using Python standard library
    modules: `xml.etree.ElementTree` and `xml.dom.minidom`.

    Parameters
    ----------
    io : str or file-like
        This can be either a string of raw XML, a valid URL,
        file or file-like object.

    index : bool
        Whether to include index in xml document.

    row_name : str
        Name for root of xml document. Default is 'data'.

    root_name : str
        Name for row elemens of xml document. Default is 'row'.

    na_rep : str
        Missing data representation.

    attrs_cols : list
        List of columns to write as attributes in row element.

    elem_cols : list
        List of columns to write as children in row element.

    namespacess : dict
        The namespaces to define in XML document as dicts with key
        being namespace and value the URI.

    prefix : str
        The prefix for each element in XML document including root.

    encoding : str
        Encoding of xml object or document.

    xml_declaration : bool
        Whether to include xml declaration at top line item in xml.

    pretty_print : bool
        Whether to write xml document with line breaks and indentation.

    stylesheet : str or file-like
        A URL, file, file-like object, or a raw string containing XSLT.

    See also
    --------
    pandas.io.formats.xml.EtreeXMLFormatter

    Notes
    -----
    This class serves as default option. If user does not have `lxml`
    installed, `to_xml` will fall back with EtreeXMLFormatter.
    """

    def __init__(
        self,
        formatter: DataFrameFormatter,
        io: Optional[FilePathOrBuffer[str]] = None,
        index: Optional[bool] = True,
        root_name: Optional[str] = "data",
        row_name: Optional[str] = "row",
        na_rep: Optional[str] = None,
        attr_cols: Optional[Union[str, List[str]]] = None,
        elem_cols: Optional[Union[str, List[str]]] = None,
        namespaces: Optional[Dict[str, str]] = None,
        prefix: Optional[str] = None,
        encoding: Optional[str] = "utf-8",
        xml_declaration: Optional[bool] = True,
        pretty_print: Optional[bool] = True,
        stylesheet: Optional[FilePathOrBuffer[str]] = None,
    ) -> None:
        self.fmt = formatter
        self.io = io
        self.index = index
        self.root_name = root_name
        self.row_name = row_name
        self.na_rep = na_rep
        self.attr_cols = attr_cols
        self.elem_cols = elem_cols
        self.namespaces = namespaces
        self.prefix = prefix
        self.encoding = encoding
        self.xml_declaration = xml_declaration
        self.pretty_print = pretty_print
        self.stylesheet = stylesheet

        self.validate_columns()
        self.validate_encoding()
        self.orig_cols = self.fmt.frame.columns.tolist()
        self.frame_dicts = self.process_dataframe()
        self.prefix_uri = self.get_prefix_uri()

        self.convert_empty_str_key()
        self.handle_indexes()

    def build_tree(self) -> bytes:
        """
        Build tree from  data.

        This method initializes the root and builds attributes and elements
        with optional namespaces.
        """
        from lxml.etree import Element, SubElement, tostring

        self.root = Element(f"{self.prefix_uri}{self.root_name}", nsmap=self.namespaces)

        for k, d in self.frame_dicts.items():
            self.d = d
            self.elem_row = SubElement(self.root, f"{self.prefix_uri}{self.row_name}")

            if self.attr_cols:
                self.build_attribs()

            if self.elem_cols:
                self.build_elems()

            if not self.attr_cols and not self.elem_cols:
                self.elem_cols = list(self.frame_dicts[0].keys())
                self.build_elems()

        self.out_xml = tostring(
            self.root,
            pretty_print=self.pretty_print,
            method="xml",
            encoding=self.encoding,
            xml_declaration=self.xml_declaration,
        )

        if self.stylesheet:
            self.out_xml = self.transform_doc()

        return self.out_xml

    def validate_columns(self) -> None:
        """
        Validate elems_cols and attrs_cols.

        This method will check if columns is list-like.

        Raises
        ------
        ValueError
            * If value is not a list and less then length of nodes.
        """
        if self.attr_cols and not is_list_like(self.attr_cols):
            raise TypeError(
                f"{type(self.attr_cols).__name__} is not a valid type for attr_cols"
            )

        if self.elem_cols and not is_list_like(self.elem_cols):
            raise TypeError(
                f"{type(self.elem_cols).__name__} is not a valid type for elem_cols"
            )

    def validate_encoding(self) -> None:
        """
        Validate encoding.

        This method will check if encoding is among listed under codecs.

        Raises
        ------
        LookupError
            * If encoding is not available in codecs.
        """

        try:
            codecs.lookup(self.encoding)
        except LookupError as e:
            raise e

    def process_dataframe(self) -> dict:
        """
        Adjust Data Frame to fit xml output.

        This method will adjust underlying data frame for xml output,
        including replacing missing entities and including indexes.
        """

        na_dict = {"None": self.na_rep, "NaN": self.na_rep, "nan": self.na_rep}

        df = (
            (self.fmt.frame.reset_index().applymap(str).replace(na_dict))
            if self.index
            else self.fmt.frame.applymap(str).replace(na_dict)
        )

        return df.to_dict(orient="index")

    def convert_empty_str_key(self) -> None:
        """
        Replace zero-lengh string in `namespaces`.

        This method will replce '' with None to align to `lxml`
        requirement that empty string prefixes are not allowed.
        """

        if self.namespaces and "" in self.namespaces.keys():
            self.namespaces[None] = self.namespaces.pop("", "default")

    def handle_indexes(self) -> None:
        """
        Handle indexes.

        This method will add indexes into attr_cols or elem_cols.
        """
        indexes = [x for x in self.frame_dicts[0].keys() if x not in self.orig_cols]

        if self.attr_cols and self.index:
            self.attr_cols = list(indexes) + self.attr_cols

        if self.elem_cols and self.index:
            self.elem_cols = list(indexes) + self.elem_cols

    def get_prefix_uri(self) -> str:
        """
        Get uri of namespace prefix.

        This method retrieves corresponding URI to prefix in namespaces.

        Raises
        ------
        ValueError
            *If prefix is not included in namespace dict.
        """

        uri = ""
        if self.namespaces:
            if self.prefix:
                try:
                    uri = f"{{{self.namespaces[self.prefix]}}}"
                except (KeyError):
                    raise KeyError("prefix is not included in namespaces")
            else:
                uri = f'{{{self.namespaces[""]}}}'

        return uri

    def build_attribs(self) -> None:
        """
        Create attributes of row.

        This method adds attributes using attr_cols to row element and
        works with tuples for multindex or hierarchical columns.
        """
        for col in self.attr_cols:
            flat_col = col
            if isinstance(col, tuple):
                flat_col = (
                    "".join(str(c) for c in col).strip()
                    if "" in col
                    else "_".join(str(c) for c in col).strip()
                )

            attr_name = f"{self.prefix_uri}{flat_col}"
            try:
                if self.d[col] is not None:
                    self.elem_row.attrib[attr_name] = self.d[col]
            except KeyError:
                raise KeyError(f"no valid column, {col}")

    def build_elems(self) -> None:
        """
        Create child elements of row.

        This method adds child elements using elem_cols to row element and
        works with tuples for multindex or hierarchical columns.
        """
        from lxml.etree import SubElement

        for col in self.elem_cols:
            flat_col = col
            if isinstance(col, tuple):
                flat_col = (
                    "".join(str(c) for c in col).strip()
                    if "" in col
                    else "_".join(str(c) for c in col).strip()
                )

            elem_name = f"{self.prefix_uri}{flat_col}"
            try:
                val = None if self.d[col] in [None, ""] else str(self.d[col])
                SubElement(self.elem_row, elem_name).text = val
            except KeyError:
                raise KeyError(f"no valid column, {col}")

    def convert_io(self) -> Union[None, str]:
        """
        Convert stylesheet object to string.

        This method will convert stylesheet object into a string or keep
        as string, depending on object type.
        """

        obj = None

        if isinstance(self.stylesheet, str):
            obj = self.stylesheet

        if isinstance(self.stylesheet, bytes):
            obj = self.stylesheet.decode(self.encoding)

        if isinstance(self.stylesheet, io.StringIO):
            obj = self.stylesheet.getvalue()

        if isinstance(self.stylesheet, io.BytesIO):
            obj = self.stylesheet.getvalue().decode(self.encoding)

        if isinstance(self.stylesheet, io.TextIOWrapper):
            obj = self.stylesheet.read()

        if isinstance(self.stylesheet, io.BufferedReader):
            obj = self.stylesheet.read().decode(self.encoding)

        return obj

    def parse_doc(self):
        """
        Build tree from stylesheet.

        This method will parse stylesheet object into tree for parsing
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

        current_doc = self.convert_io()
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

    def transform_doc(self) -> bytes:
        """
        Transform original tree using stylesheet.

        This method will transform built tree with XSLT script.
        """

        from lxml.etree import XSLT, XSLTApplyError, XSLTParseError

        xsl_doc = self.parse_doc()

        try:
            transformer = XSLT(xsl_doc)
            new_doc = transformer(self.root)

        except (XSLTApplyError, XSLTParseError) as e:
            raise e

        return bytes(new_doc)

    def write_output(self) -> Optional[str]:
        xml_doc = self.build_tree()

        try:
            if self.io:
                with open(self.io, "wb") as f:
                    f.write(xml_doc)
                xml_doc = None
            else:
                xml_doc = xml_doc.decode(self.encoding).rstrip()

        except (UnicodeDecodeError, OSError) as e:
            raise e

        return xml_doc
