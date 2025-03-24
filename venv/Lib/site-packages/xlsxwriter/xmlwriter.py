###############################################################################
#
# XMLwriter - A base class for XlsxWriter classes.
#
# Used in conjunction with XlsxWriter.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

# pylint: disable=dangerous-default-value

# Standard packages.
import re
from io import StringIO

# Compile performance critical regular expressions.
re_control_chars_1 = re.compile("(_x[0-9a-fA-F]{4}_)")
re_control_chars_2 = re.compile(r"([\x00-\x08\x0b-\x1f])")
xml_escapes = re.compile('["&<>\n]')


class XMLwriter:
    """
    Simple XML writer class.

    """

    def __init__(self):
        self.fh = None
        self.internal_fh = False

    def _set_filehandle(self, filehandle):
        # Set the writer filehandle directly. Mainly for testing.
        self.fh = filehandle
        self.internal_fh = False

    def _set_xml_writer(self, filename):
        # Set the XML writer filehandle for the object.
        if isinstance(filename, StringIO):
            self.internal_fh = False
            self.fh = filename
        else:
            self.internal_fh = True
            # pylint: disable-next=consider-using-with
            self.fh = open(filename, "w", encoding="utf-8")

    def _xml_close(self):
        # Close the XML filehandle if we created it.
        if self.internal_fh:
            self.fh.close()

    def _xml_declaration(self):
        # Write the XML declaration.
        self.fh.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')

    def _xml_start_tag(self, tag, attributes=[]):
        # Write an XML start tag with optional attributes.
        for key, value in attributes:
            value = self._escape_attributes(value)
            tag += f' {key}="{value}"'

        self.fh.write(f"<{tag}>")

    def _xml_start_tag_unencoded(self, tag, attributes=[]):
        # Write an XML start tag with optional, unencoded, attributes.
        # This is a minor speed optimization for elements that don't
        # need encoding.
        for key, value in attributes:
            tag += f' {key}="{value}"'

        self.fh.write(f"<{tag}>")

    def _xml_end_tag(self, tag):
        # Write an XML end tag.
        self.fh.write(f"</{tag}>")

    def _xml_empty_tag(self, tag, attributes=[]):
        # Write an empty XML tag with optional attributes.
        for key, value in attributes:
            value = self._escape_attributes(value)
            tag += f' {key}="{value}"'

        self.fh.write(f"<{tag}/>")

    def _xml_empty_tag_unencoded(self, tag, attributes=[]):
        # Write an empty XML tag with optional, unencoded, attributes.
        # This is a minor speed optimization for elements that don't
        # need encoding.
        for key, value in attributes:
            tag += f' {key}="{value}"'

        self.fh.write(f"<{tag}/>")

    def _xml_data_element(self, tag, data, attributes=[]):
        # Write an XML element containing data with optional attributes.
        end_tag = tag

        for key, value in attributes:
            value = self._escape_attributes(value)
            tag += f' {key}="{value}"'

        data = self._escape_data(data)
        data = self._escape_control_characters(data)

        self.fh.write(f"<{tag}>{data}</{end_tag}>")

    def _xml_string_element(self, index, attributes=[]):
        # Optimized tag writer for <c> cell string elements in the inner loop.
        attr = ""

        for key, value in attributes:
            value = self._escape_attributes(value)
            attr += f' {key}="{value}"'

        self.fh.write(f'<c{attr} t="s"><v>{index}</v></c>')

    def _xml_si_element(self, string, attributes=[]):
        # Optimized tag writer for shared strings <si> elements.
        attr = ""

        for key, value in attributes:
            value = self._escape_attributes(value)
            attr += f' {key}="{value}"'

        string = self._escape_data(string)

        self.fh.write(f"<si><t{attr}>{string}</t></si>")

    def _xml_rich_si_element(self, string):
        # Optimized tag writer for shared strings <si> rich string elements.

        self.fh.write(f"<si>{string}</si>")

    def _xml_number_element(self, number, attributes=[]):
        # Optimized tag writer for <c> cell number elements in the inner loop.
        attr = ""

        for key, value in attributes:
            value = self._escape_attributes(value)
            attr += f' {key}="{value}"'

        self.fh.write(f"<c{attr}><v>{number:.16G}</v></c>")

    def _xml_formula_element(self, formula, result, attributes=[]):
        # Optimized tag writer for <c> cell formula elements in the inner loop.
        attr = ""

        for key, value in attributes:
            value = self._escape_attributes(value)
            attr += f' {key}="{value}"'

        formula = self._escape_data(formula)
        result = self._escape_data(result)
        self.fh.write(f"<c{attr}><f>{formula}</f><v>{result}</v></c>")

    def _xml_inline_string(self, string, preserve, attributes=[]):
        # Optimized tag writer for inlineStr cell elements in the inner loop.
        attr = ""
        t_attr = ""

        # Set the <t> attribute to preserve whitespace.
        if preserve:
            t_attr = ' xml:space="preserve"'

        for key, value in attributes:
            value = self._escape_attributes(value)
            attr += f' {key}="{value}"'

        string = self._escape_data(string)

        self.fh.write(f'<c{attr} t="inlineStr"><is><t{t_attr}>{string}</t></is></c>')

    def _xml_rich_inline_string(self, string, attributes=[]):
        # Optimized tag writer for rich inlineStr in the inner loop.
        attr = ""

        for key, value in attributes:
            value = self._escape_attributes(value)
            attr += f' {key}="{value}"'

        self.fh.write(f'<c{attr} t="inlineStr"><is>{string}</is></c>')

    def _escape_attributes(self, attribute):
        # Escape XML characters in attributes.
        try:
            if not xml_escapes.search(attribute):
                return attribute
        except TypeError:
            return attribute

        attribute = (
            attribute.replace("&", "&amp;")
            .replace('"', "&quot;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace("\n", "&#xA;")
        )
        return attribute

    def _escape_data(self, data):
        # Escape XML characters in data sections of tags.  Note, this
        # is different from _escape_attributes() in that double quotes
        # are not escaped by Excel.
        try:
            if not xml_escapes.search(data):
                return data
        except TypeError:
            return data

        data = data.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        return data

    @staticmethod
    def _escape_control_characters(data):
        # Excel escapes control characters with _xHHHH_ and also escapes any
        # literal strings of that type by encoding the leading underscore.
        # So "\0" -> _x0000_ and "_x0000_" -> _x005F_x0000_.
        # The following substitutions deal with those cases.
        try:
            # Escape the escape.
            data = re_control_chars_1.sub(r"_x005F\1", data)
        except TypeError:
            return data

        # Convert control character to the _xHHHH_ escape.
        data = re_control_chars_2.sub(
            lambda match: f"_x{ord(match.group(1)):04X}_", data
        )

        # Escapes non characters in strings.
        data = data.replace("\uFFFE", "_xFFFE_").replace("\uFFFF", "_xFFFF_")

        return data
