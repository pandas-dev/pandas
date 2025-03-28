###############################################################################
#
# RichValueTypes - A class for writing the Excel XLSX rdRichValueTypes.xml file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

# Package imports.
from . import xmlwriter


class RichValueTypes(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX rdRichValueTypes.xml file.


    """

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self):
        # Assemble and write the XML file.

        # Write the XML declaration.
        self._xml_declaration()

        # Write the rvTypesInfo element.
        self._write_rv_types_info()

        # Write the global element.
        self._write_global()

        self._xml_end_tag("rvTypesInfo")

        # Close the file.
        self._xml_close()

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################

    def _write_rv_types_info(self):
        # Write the <rvTypesInfo> element.
        xmlns = "http://schemas.microsoft.com/office/spreadsheetml/2017/richdata2"
        xmlns_x = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        xmlns_mc = "http://schemas.openxmlformats.org/markup-compatibility/2006"
        mc_ignorable = "x"

        attributes = [
            ("xmlns", xmlns),
            ("xmlns:mc", xmlns_mc),
            ("mc:Ignorable", mc_ignorable),
            ("xmlns:x", xmlns_x),
        ]

        self._xml_start_tag("rvTypesInfo", attributes)

    def _write_global(self):
        # Write the <global> element.
        key_flags = [
            ["_Self", ["ExcludeFromFile", "ExcludeFromCalcComparison"]],
            ["_DisplayString", ["ExcludeFromCalcComparison"]],
            ["_Flags", ["ExcludeFromCalcComparison"]],
            ["_Format", ["ExcludeFromCalcComparison"]],
            ["_SubLabel", ["ExcludeFromCalcComparison"]],
            ["_Attribution", ["ExcludeFromCalcComparison"]],
            ["_Icon", ["ExcludeFromCalcComparison"]],
            ["_Display", ["ExcludeFromCalcComparison"]],
            ["_CanonicalPropertyNames", ["ExcludeFromCalcComparison"]],
            ["_ClassificationId", ["ExcludeFromCalcComparison"]],
        ]

        self._xml_start_tag("global")
        self._xml_start_tag("keyFlags")

        for key_flag in key_flags:
            # Write the key element.
            self._write_key(key_flag)

        self._xml_end_tag("keyFlags")
        self._xml_end_tag("global")

    def _write_key(self, key_flag):
        # Write the <key> element.
        name = key_flag[0]
        attributes = [("name", name)]

        self._xml_start_tag("key", attributes)

        # Write the flag element.
        for name in key_flag[1]:
            self._write_flag(name)

        self._xml_end_tag("key")

    def _write_flag(self, name):
        # Write the <flag> element.
        attributes = [
            ("name", name),
            ("value", "1"),
        ]

        self._xml_empty_tag("flag", attributes)
