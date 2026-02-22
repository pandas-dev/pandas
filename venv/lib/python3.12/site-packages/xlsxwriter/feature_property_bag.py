###############################################################################
#
# FeaturePropertyBag - A class for writing the Excel XLSX featurePropertyBag.xml
#                      file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

# Package imports.
from . import xmlwriter


class FeaturePropertyBag(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX FeaturePropertyBag file.


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

        self.feature_property_bags = set()

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self) -> None:
        # Assemble and write the XML file.

        # Write the XML declaration.
        self._xml_declaration()

        # Write the FeaturePropertyBags element.
        self._write_feature_property_bags()

        # Write the Checkbox bag element.
        self._write_checkbox_bag()

        # Write the XFControls bag element.
        self._write_xf_control_bag()

        # Write the XFComplement bag element.
        self._write_xf_compliment_bag()

        # Write the XFComplements bag element.
        self._write_xf_compliments_bag()

        # Write the DXFComplements bag element.
        if "DXFComplements" in self.feature_property_bags:
            self._write_dxf_compliments_bag()

        self._xml_end_tag("FeaturePropertyBags")

        # Close the file.
        self._xml_close()

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################

    def _write_feature_property_bags(self) -> None:
        # Write the <FeaturePropertyBags> element.

        xmlns = (
            "http://schemas.microsoft.com/office/spreadsheetml/2022/featurepropertybag"
        )

        attributes = [("xmlns", xmlns)]

        self._xml_start_tag("FeaturePropertyBags", attributes)

    def _write_checkbox_bag(self) -> None:
        # Write the Checkbox <bag> element.
        attributes = [("type", "Checkbox")]

        self._xml_empty_tag("bag", attributes)

    def _write_xf_control_bag(self) -> None:
        # Write the XFControls<bag> element.
        attributes = [("type", "XFControls")]

        self._xml_start_tag("bag", attributes)

        # Write the bagId element.
        self._write_bag_id("CellControl", 0)

        self._xml_end_tag("bag")

    def _write_xf_compliment_bag(self) -> None:
        # Write the XFComplement <bag> element.
        attributes = [("type", "XFComplement")]

        self._xml_start_tag("bag", attributes)

        # Write the bagId element.
        self._write_bag_id("XFControls", 1)

        self._xml_end_tag("bag")

    def _write_xf_compliments_bag(self) -> None:
        # Write the XFComplements <bag> element.
        attributes = [
            ("type", "XFComplements"),
            ("extRef", "XFComplementsMapperExtRef"),
        ]

        self._xml_start_tag("bag", attributes)
        self._xml_start_tag("a", [("k", "MappedFeaturePropertyBags")])

        self._write_bag_id("", 2)

        self._xml_end_tag("a")
        self._xml_end_tag("bag")

    def _write_dxf_compliments_bag(self) -> None:
        # Write the DXFComplements <bag> element.
        attributes = [
            ("type", "DXFComplements"),
            ("extRef", "DXFComplementsMapperExtRef"),
        ]

        self._xml_start_tag("bag", attributes)
        self._xml_start_tag("a", [("k", "MappedFeaturePropertyBags")])

        self._write_bag_id("", 2)

        self._xml_end_tag("a")
        self._xml_end_tag("bag")

    def _write_bag_id(self, key, bag_id) -> None:
        # Write the <bagId> element.
        attributes = []

        if key:
            attributes = [("k", key)]

        self._xml_data_element("bagId", bag_id, attributes)
