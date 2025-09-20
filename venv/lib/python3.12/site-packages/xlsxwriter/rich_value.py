###############################################################################
#
# RichValue - A class for writing the Excel XLSX rdrichvalue.xml file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

from xlsxwriter.image import Image

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

    def __init__(self) -> None:
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

    def _assemble_xml_file(self) -> None:
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
    def _write_rv_data(self) -> None:
        # Write the <rvData> element.
        xmlns = "http://schemas.microsoft.com/office/spreadsheetml/2017/richdata"

        attributes = [
            ("xmlns", xmlns),
            ("count", len(self.embedded_images)),
        ]

        self._xml_start_tag("rvData", attributes)

        for index, image in enumerate(self.embedded_images):
            # Write the rv element.
            self._write_rv(index, image)

    def _write_rv(self, index, image: Image) -> None:
        # Write the <rv> element.
        attributes = [("s", 0)]
        value = 5

        if image.decorative:
            value = 6

        self._xml_start_tag("rv", attributes)

        # Write the v elements.
        self._write_v(index)
        self._write_v(value)

        if image.description:
            self._write_v(image.description)

        self._xml_end_tag("rv")

    def _write_v(self, data) -> None:
        # Write the <v> element.
        self._xml_data_element("v", data)
