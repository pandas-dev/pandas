###############################################################################
#
# Core - A class for writing the Excel XLSX Worksheet file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

from datetime import datetime, timezone
from typing import Dict, Union

from . import xmlwriter


class Core(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX Core file.


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

        self.properties = {}
        self.iso_date = ""

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self) -> None:
        # Assemble and write the XML file.

        # Set the creation date for the file.
        date = self.properties.get("created")
        if not isinstance(date, datetime):
            date = datetime.now(timezone.utc)

        self.iso_date = date.strftime("%Y-%m-%dT%H:%M:%SZ")

        # Write the XML declaration.
        self._xml_declaration()

        self._write_cp_core_properties()
        self._write_dc_title()
        self._write_dc_subject()
        self._write_dc_creator()
        self._write_cp_keywords()
        self._write_dc_description()
        self._write_cp_last_modified_by()
        self._write_dcterms_created()
        self._write_dcterms_modified()
        self._write_cp_category()
        self._write_cp_content_status()

        self._xml_end_tag("cp:coreProperties")

        # Close the file.
        self._xml_close()

    def _set_properties(self, properties: Dict[str, Union[str, datetime]]) -> None:
        # Set the document properties.
        self.properties = properties

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################

    def _write_cp_core_properties(self) -> None:
        # Write the <cp:coreProperties> element.

        xmlns_cp = (
            "http://schemas.openxmlformats.org/package/2006/"
            + "metadata/core-properties"
        )
        xmlns_dc = "http://purl.org/dc/elements/1.1/"
        xmlns_dcterms = "http://purl.org/dc/terms/"
        xmlns_dcmitype = "http://purl.org/dc/dcmitype/"
        xmlns_xsi = "http://www.w3.org/2001/XMLSchema-instance"

        attributes = [
            ("xmlns:cp", xmlns_cp),
            ("xmlns:dc", xmlns_dc),
            ("xmlns:dcterms", xmlns_dcterms),
            ("xmlns:dcmitype", xmlns_dcmitype),
            ("xmlns:xsi", xmlns_xsi),
        ]

        self._xml_start_tag("cp:coreProperties", attributes)

    def _write_dc_creator(self) -> None:
        # Write the <dc:creator> element.
        data = self.properties.get("author", "")

        self._xml_data_element("dc:creator", data)

    def _write_cp_last_modified_by(self) -> None:
        # Write the <cp:lastModifiedBy> element.
        data = self.properties.get("author", "")

        self._xml_data_element("cp:lastModifiedBy", data)

    def _write_dcterms_created(self) -> None:
        # Write the <dcterms:created> element.
        attributes = [("xsi:type", "dcterms:W3CDTF")]
        self._xml_data_element("dcterms:created", self.iso_date, attributes)

    def _write_dcterms_modified(self) -> None:
        # Write the <dcterms:modified> element.
        attributes = [("xsi:type", "dcterms:W3CDTF")]
        self._xml_data_element("dcterms:modified", self.iso_date, attributes)

    def _write_dc_title(self) -> None:
        # Write the <dc:title> element.
        if "title" in self.properties:
            data = self.properties["title"]
        else:
            return

        self._xml_data_element("dc:title", data)

    def _write_dc_subject(self) -> None:
        # Write the <dc:subject> element.
        if "subject" in self.properties:
            data = self.properties["subject"]
        else:
            return

        self._xml_data_element("dc:subject", data)

    def _write_cp_keywords(self) -> None:
        # Write the <cp:keywords> element.
        if "keywords" in self.properties:
            data = self.properties["keywords"]
        else:
            return

        self._xml_data_element("cp:keywords", data)

    def _write_dc_description(self) -> None:
        # Write the <dc:description> element.
        if "comments" in self.properties:
            data = self.properties["comments"]
        else:
            return

        self._xml_data_element("dc:description", data)

    def _write_cp_category(self) -> None:
        # Write the <cp:category> element.
        if "category" in self.properties:
            data = self.properties["category"]
        else:
            return

        self._xml_data_element("cp:category", data)

    def _write_cp_content_status(self) -> None:
        # Write the <cp:contentStatus> element.
        if "status" in self.properties:
            data = self.properties["status"]
        else:
            return

        self._xml_data_element("cp:contentStatus", data)
