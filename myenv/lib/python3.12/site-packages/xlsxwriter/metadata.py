###############################################################################
#
# Metadata - A class for writing the Excel XLSX Metadata file.
#
# SPDX-License-Identifier: BSD-2-Clause
# Copyright 2013-2024, John McNamara, jmcnamara@cpan.org
#

from . import xmlwriter


class Metadata(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX Metadata file.


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

        super(Metadata, self).__init__()
        self.has_dynamic_functions = False
        self.has_embedded_images = False
        self.num_embedded_images = 0

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self):
        # Assemble and write the XML file.

        if self.num_embedded_images > 0:
            self.has_embedded_images = True

        # Write the XML declaration.
        self._xml_declaration()

        # Write the metadata element.
        self._write_metadata()

        # Write the metadataTypes element.
        self._write_metadata_types()

        # Write the futureMetadata elements.
        if self.has_dynamic_functions:
            self._write_cell_future_metadata()
        if self.has_embedded_images:
            self._write_value_future_metadata()

        # Write the cellMetadata element.
        if self.has_dynamic_functions:
            self._write_cell_metadata()
        if self.has_embedded_images:
            self._write_value_metadata()

        self._xml_end_tag("metadata")

        # Close the file.
        self._xml_close()

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################

    def _write_metadata(self):
        # Write the <metadata> element.
        xmlns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
        schema = "http://schemas.microsoft.com/office/spreadsheetml"

        attributes = [("xmlns", xmlns)]

        if self.has_embedded_images:
            attributes.append(("xmlns:xlrd", schema + "/2017/richdata"))

        if self.has_dynamic_functions:
            attributes.append(("xmlns:xda", schema + "/2017/dynamicarray"))

        self._xml_start_tag("metadata", attributes)

    def _write_metadata_types(self):
        # Write the <metadataTypes> element.
        count = 0

        if self.has_dynamic_functions:
            count += 1
        if self.has_embedded_images:
            count += 1

        attributes = [("count", count)]

        self._xml_start_tag("metadataTypes", attributes)

        # Write the metadataType element.
        if self.has_dynamic_functions:
            self._write_cell_metadata_type()
        if self.has_embedded_images:
            self._write_value_metadata_type()

        self._xml_end_tag("metadataTypes")

    def _write_cell_metadata_type(self):
        # Write the <metadataType> element.
        attributes = [
            ("name", "XLDAPR"),
            ("minSupportedVersion", 120000),
            ("copy", 1),
            ("pasteAll", 1),
            ("pasteValues", 1),
            ("merge", 1),
            ("splitFirst", 1),
            ("rowColShift", 1),
            ("clearFormats", 1),
            ("clearComments", 1),
            ("assign", 1),
            ("coerce", 1),
            ("cellMeta", 1),
        ]

        self._xml_empty_tag("metadataType", attributes)

    def _write_value_metadata_type(self):
        # Write the <metadataType> element.
        attributes = [
            ("name", "XLRICHVALUE"),
            ("minSupportedVersion", 120000),
            ("copy", 1),
            ("pasteAll", 1),
            ("pasteValues", 1),
            ("merge", 1),
            ("splitFirst", 1),
            ("rowColShift", 1),
            ("clearFormats", 1),
            ("clearComments", 1),
            ("assign", 1),
            ("coerce", 1),
        ]

        self._xml_empty_tag("metadataType", attributes)

    def _write_cell_future_metadata(self):
        # Write the <futureMetadata> element.
        attributes = [
            ("name", "XLDAPR"),
            ("count", 1),
        ]

        self._xml_start_tag("futureMetadata", attributes)
        self._xml_start_tag("bk")
        self._xml_start_tag("extLst")
        self._write_cell_ext()
        self._xml_end_tag("extLst")
        self._xml_end_tag("bk")
        self._xml_end_tag("futureMetadata")

    def _write_value_future_metadata(self):
        # Write the <futureMetadata> element.
        attributes = [
            ("name", "XLRICHVALUE"),
            ("count", self.num_embedded_images),
        ]

        self._xml_start_tag("futureMetadata", attributes)

        for index in range(self.num_embedded_images):
            self._xml_start_tag("bk")
            self._xml_start_tag("extLst")
            self._write_value_ext(index)
            self._xml_end_tag("extLst")
            self._xml_end_tag("bk")

        self._xml_end_tag("futureMetadata")

    def _write_cell_ext(self):
        # Write the <ext> element.
        attributes = [("uri", "{bdbb8cdc-fa1e-496e-a857-3c3f30c029c3}")]

        self._xml_start_tag("ext", attributes)

        # Write the xda:dynamicArrayProperties element.
        self._write_xda_dynamic_array_properties()

        self._xml_end_tag("ext")

    def _write_xda_dynamic_array_properties(self):
        # Write the <xda:dynamicArrayProperties> element.
        attributes = [
            ("fDynamic", 1),
            ("fCollapsed", 0),
        ]

        self._xml_empty_tag("xda:dynamicArrayProperties", attributes)

    def _write_value_ext(self, index):
        # Write the <ext> element.
        attributes = [("uri", "{3e2802c4-a4d2-4d8b-9148-e3be6c30e623}")]

        self._xml_start_tag("ext", attributes)

        # Write the xlrd:rvb element.
        self._write_xlrd_rvb(index)

        self._xml_end_tag("ext")

    def _write_xlrd_rvb(self, index):
        # Write the <xlrd:rvb> element.
        attributes = [("i", index)]

        self._xml_empty_tag("xlrd:rvb", attributes)

    def _write_cell_metadata(self):
        # Write the <cellMetadata> element.
        attributes = [("count", 1)]

        self._xml_start_tag("cellMetadata", attributes)
        self._xml_start_tag("bk")

        # Write the rc element.
        self._write_rc(1, 0)

        self._xml_end_tag("bk")
        self._xml_end_tag("cellMetadata")

    def _write_value_metadata(self):
        # Write the <valueMetadata> element.
        count = self.num_embedded_images
        type = 1

        if self.has_dynamic_functions:
            type = 2

        attributes = [("count", count)]

        self._xml_start_tag("valueMetadata", attributes)

        # Write the rc elements.
        for index in range(self.num_embedded_images):
            self._xml_start_tag("bk")
            self._write_rc(type, index)
            self._xml_end_tag("bk")

        self._xml_end_tag("valueMetadata")

    def _write_rc(self, type, index):
        # Write the <rc> element.
        attributes = [
            ("t", type),
            ("v", index),
        ]

        self._xml_empty_tag("rc", attributes)
