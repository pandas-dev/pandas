###############################################################################
#
# RichValue - A class for writing the Excel XLSX rdrichvalue.xml file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

# Package imports.
from . import xmlwriter


class RichValue(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX rdrichvalue.xml file.


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
        self.embedded_images = []

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self):
        # Assemble and write the XML file.

        # Write the XML declaration.
        self._xml_declaration()

        # Write the rvData element.
        self._write_rv_data()

        self._xml_end_tag("rvData")

        # Close the file.
        self._xml_close()

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################
    def _write_rv_data(self):
        # Write the <rvData> element.
        xmlns = "http://schemas.microsoft.com/office/spreadsheetml/2017/richdata"

        attributes = [
            ("xmlns", xmlns),
            ("count", len(self.embedded_images)),
        ]

        self._xml_start_tag("rvData", attributes)

        for index, image_data in enumerate(self.embedded_images):
            # Write the rv element.
            self._write_rv(index, image_data[3], image_data[4])

    def _write_rv(self, index, description, decorative):
        # Write the <rv> element.
        attributes = [("s", 0)]
        value = 5

        if decorative:
            value = 6

        self._xml_start_tag("rv", attributes)

        # Write the v elements.
        self._write_v(index)
        self._write_v(value)

        if description:
            self._write_v(description)

        self._xml_end_tag("rv")

    def _write_v(self, data):
        # Write the <v> element.
        self._xml_data_element("v", data)
