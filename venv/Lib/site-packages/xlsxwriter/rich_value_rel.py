###############################################################################
#
# RichValueRel - A class for writing the Excel XLSX richValueRel.xml file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

# Package imports.
from . import xmlwriter


class RichValueRel(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX richValueRel.xml file.


    """

    ###########################################################################
    #
    # Public API.
    #
    ###########################################################################

    def __init__(self):
        """
        Constructor.

        """

        super().__init__()
        self.num_embedded_images = 0

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self):
        # Assemble and write the XML file.

        # Write the XML declaration.
        self._xml_declaration()

        # Write the richValueRels element.
        self._write_rich_value_rels()

        self._xml_end_tag("richValueRels")

        # Close the file.
        self._xml_close()

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################
    def _write_rich_value_rels(self):
        # Write the <richValueRels> element.
        xmlns = "http://schemas.microsoft.com/office/spreadsheetml/2022/richvaluerel"
        xmlns_r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

        attributes = [
            ("xmlns", xmlns),
            ("xmlns:r", xmlns_r),
        ]

        self._xml_start_tag("richValueRels", attributes)

        # Write the rel elements.
        for index in range(self.num_embedded_images):
            self._write_rel(index + 1)

    def _write_rel(self, index):
        # Write the <rel> element.
        r_id = f"rId{index}"
        attributes = [("r:id", r_id)]

        self._xml_empty_tag("rel", attributes)
