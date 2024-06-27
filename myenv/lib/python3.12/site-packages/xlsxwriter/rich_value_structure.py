###############################################################################
#
# RichValueStructure - A class for writing the Excel XLSX rdrichvaluestructure.xml file.
#
# SPDX-License-Identifier: BSD-2-Clause
# Copyright 2013-2024, John McNamara, jmcnamara@cpan.org
#

# Package imports.
from . import xmlwriter


class RichValueStructure(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX rdrichvaluestructure.xml file.


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

        super(RichValueStructure, self).__init__()
        self.has_embedded_descriptions = False

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self):
        # Assemble and write the XML file.

        # Write the XML declaration.
        self._xml_declaration()

        # Write the rvStructures element.
        self._write_rv_structures()

        self._xml_end_tag("rvStructures")

        # Close the file.
        self._xml_close()

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################
    def _write_rv_structures(self):
        # Write the <rvStructures> element.
        xmlns = "http://schemas.microsoft.com/office/spreadsheetml/2017/richdata"
        count = "1"

        attributes = [
            ("xmlns", xmlns),
            ("count", count),
        ]

        self._xml_start_tag("rvStructures", attributes)

        # Write the s element.
        self._write_s()

    def _write_s(self):
        # Write the <s> element.
        t = "_localImage"
        attributes = [("t", t)]

        self._xml_start_tag("s", attributes)

        # Write the k elements.
        self._write_k("_rvRel:LocalImageIdentifier", "i")
        self._write_k("CalcOrigin", "i")

        if self.has_embedded_descriptions:
            self._write_k("Text", "s")

        self._xml_end_tag("s")

    def _write_k(self, name, type):
        # Write the <k> element.
        attributes = [
            ("n", name),
            ("t", type),
        ]

        self._xml_empty_tag("k", attributes)
