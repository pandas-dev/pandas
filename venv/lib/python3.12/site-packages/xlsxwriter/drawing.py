###############################################################################
#
# Drawing - A class for writing the Excel XLSX Drawing file.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

from enum import Enum

from xlsxwriter.color import Color
from xlsxwriter.url import Url

from . import xmlwriter
from .shape import Shape


class DrawingTypes(Enum):
    """
    Enum to represent different types of drawings in a worksheet.
    """

    NONE = 0
    CHART = 1
    IMAGE = 2
    SHAPE = 3


class DrawingInfo:
    """
    An internal class to represent a drawing object in an Excel worksheet.

    """

    def __init__(self) -> None:
        """
        Initialize a DrawingType instance with default values.
        """
        self._drawing_type = DrawingTypes.NONE
        self._anchor_type = None
        self._dimensions = []
        self._width = 0
        self._height = 0
        self._shape = None
        self._anchor = None
        self._url = None
        self._rel_index = 0
        self._name = None
        self._description = None
        self._decorative = False


class Drawing(xmlwriter.XMLwriter):
    """
    A class for writing the Excel XLSX Drawing file.


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

        self.drawings = []
        self.embedded = 0
        self.orientation = 0

    ###########################################################################
    #
    # Private API.
    #
    ###########################################################################

    def _assemble_xml_file(self) -> None:
        # Assemble and write the XML file.

        # Write the XML declaration.
        self._xml_declaration()

        # Write the xdr:wsDr element.
        self._write_drawing_workspace()

        if self.embedded:
            index = 0
            for drawing in self.drawings:
                # Write the xdr:twoCellAnchor element.
                index += 1
                self._write_two_cell_anchor(index, drawing)

        else:
            # Write the xdr:absoluteAnchor element.
            drawing = DrawingInfo()
            drawing._rel_index = 1
            self._write_absolute_anchor(1, drawing)

        self._xml_end_tag("xdr:wsDr")

        # Close the file.
        self._xml_close()

    def _add_drawing_object(self, drawing_object: DrawingInfo) -> None:
        # Add a chart, image or shape sub object to the drawing.
        self.drawings.append(drawing_object)

    ###########################################################################
    #
    # XML methods.
    #
    ###########################################################################

    def _write_drawing_workspace(self) -> None:
        # Write the <xdr:wsDr> element.
        schema = "http://schemas.openxmlformats.org/drawingml/"
        xmlns_xdr = schema + "2006/spreadsheetDrawing"
        xmlns_a = schema + "2006/main"

        attributes = [
            ("xmlns:xdr", xmlns_xdr),
            ("xmlns:a", xmlns_a),
        ]

        self._xml_start_tag("xdr:wsDr", attributes)

    def _write_two_cell_anchor(self, index: int, drawing: DrawingInfo) -> None:
        # Write the <xdr:twoCellAnchor> element.
        dimensions = drawing._dimensions
        col_from = dimensions[0]
        row_from = dimensions[1]
        col_from_offset = dimensions[2]
        row_from_offset = dimensions[3]
        col_to = dimensions[4]
        row_to = dimensions[5]
        col_to_offset = dimensions[6]
        row_to_offset = dimensions[7]
        col_absolute = dimensions[8]
        row_absolute = dimensions[9]

        attributes = []

        # Add attribute for positioning.
        if drawing._anchor == 2:
            attributes.append(("editAs", "oneCell"))
        elif drawing._anchor == 3:
            attributes.append(("editAs", "absolute"))

        # Add editAs attribute for shapes.
        if drawing._shape and drawing._shape.edit_as:
            attributes.append(("editAs", drawing._shape.edit_as))

        self._xml_start_tag("xdr:twoCellAnchor", attributes)

        # Write the xdr:from element.
        self._write_from(col_from, row_from, col_from_offset, row_from_offset)

        # Write the xdr:from element.
        self._write_to(col_to, row_to, col_to_offset, row_to_offset)

        if drawing._drawing_type == DrawingTypes.CHART:
            # Graphic frame.
            # Write the xdr:graphicFrame element for charts.
            self._write_graphic_frame(index, drawing)
        elif drawing._drawing_type == DrawingTypes.IMAGE:
            # Write the xdr:pic element.
            self._write_pic(index, col_absolute, row_absolute, drawing)
        else:
            # Write the xdr:sp element for shapes.
            self._write_sp(index, col_absolute, row_absolute, drawing)

        # Write the xdr:clientData element.
        self._write_client_data()

        self._xml_end_tag("xdr:twoCellAnchor")

    def _write_absolute_anchor(self, index: int, drawing: DrawingInfo) -> None:
        self._xml_start_tag("xdr:absoluteAnchor")
        # Write the <xdr:absoluteAnchor> element.

        # Different coordinates for horizontal (= 0) and vertical (= 1).
        if self.orientation == 0:
            # Write the xdr:pos element.
            self._write_pos(0, 0)

            # Write the xdr:ext element.
            self._write_xdr_ext(9308969, 6078325)

        else:
            # Write the xdr:pos element.
            self._write_pos(0, -47625)

            # Write the xdr:ext element.
            self._write_xdr_ext(6162675, 6124575)

        # Write the xdr:graphicFrame element.
        self._write_graphic_frame(index, drawing)

        # Write the xdr:clientData element.
        self._write_client_data()

        self._xml_end_tag("xdr:absoluteAnchor")

    def _write_from(self, col: int, row: int, col_offset, row_offset) -> None:
        # Write the <xdr:from> element.
        self._xml_start_tag("xdr:from")

        # Write the xdr:col element.
        self._write_col(col)

        # Write the xdr:colOff element.
        self._write_col_off(col_offset)

        # Write the xdr:row element.
        self._write_row(row)

        # Write the xdr:rowOff element.
        self._write_row_off(row_offset)

        self._xml_end_tag("xdr:from")

    def _write_to(self, col: int, row: int, col_offset, row_offset) -> None:
        # Write the <xdr:to> element.
        self._xml_start_tag("xdr:to")

        # Write the xdr:col element.
        self._write_col(col)

        # Write the xdr:colOff element.
        self._write_col_off(col_offset)

        # Write the xdr:row element.
        self._write_row(row)

        # Write the xdr:rowOff element.
        self._write_row_off(row_offset)

        self._xml_end_tag("xdr:to")

    def _write_col(self, data) -> None:
        # Write the <xdr:col> element.
        self._xml_data_element("xdr:col", data)

    def _write_col_off(self, data) -> None:
        # Write the <xdr:colOff> element.
        self._xml_data_element("xdr:colOff", data)

    def _write_row(self, data) -> None:
        # Write the <xdr:row> element.
        self._xml_data_element("xdr:row", data)

    def _write_row_off(self, data) -> None:
        # Write the <xdr:rowOff> element.
        self._xml_data_element("xdr:rowOff", data)

    def _write_pos(self, x, y) -> None:
        # Write the <xdr:pos> element.

        attributes = [("x", x), ("y", y)]

        self._xml_empty_tag("xdr:pos", attributes)

    def _write_xdr_ext(self, cx, cy) -> None:
        # Write the <xdr:ext> element.

        attributes = [("cx", cx), ("cy", cy)]

        self._xml_empty_tag("xdr:ext", attributes)

    def _write_graphic_frame(self, index: int, drawing: DrawingInfo) -> None:
        # Write the <xdr:graphicFrame> element.
        attributes = [("macro", "")]

        self._xml_start_tag("xdr:graphicFrame", attributes)

        # Write the xdr:nvGraphicFramePr element.
        self._write_nv_graphic_frame_pr(index, drawing)

        # Write the xdr:xfrm element.
        self._write_xfrm()

        # Write the a:graphic element.
        self._write_atag_graphic(drawing._rel_index)

        self._xml_end_tag("xdr:graphicFrame")

    def _write_nv_graphic_frame_pr(self, index: int, drawing: DrawingInfo) -> None:
        # Write the <xdr:nvGraphicFramePr> element.

        name = drawing._name
        if not name:
            name = "Chart " + str(index)

        self._xml_start_tag("xdr:nvGraphicFramePr")

        # Write the xdr:cNvPr element.
        self._write_c_nv_pr(index + 1, drawing, name)

        # Write the xdr:cNvGraphicFramePr element.
        self._write_c_nv_graphic_frame_pr()

        self._xml_end_tag("xdr:nvGraphicFramePr")

    def _write_c_nv_pr(self, index: int, drawing: DrawingInfo, name: str) -> None:
        # Write the <xdr:cNvPr> element.
        attributes = [("id", index), ("name", name)]

        # Add description attribute for images.
        if drawing._description and not drawing._decorative:
            attributes.append(("descr", drawing._description))

        if drawing._url or drawing._decorative:
            self._xml_start_tag("xdr:cNvPr", attributes)

            if drawing._url:
                self._write_a_hlink_click(drawing._url)

            if drawing._decorative:
                self._write_decorative()

            self._xml_end_tag("xdr:cNvPr")
        else:
            self._xml_empty_tag("xdr:cNvPr", attributes)

    def _write_decorative(self) -> None:
        self._xml_start_tag("a:extLst")

        self._write_uri_ext("{FF2B5EF4-FFF2-40B4-BE49-F238E27FC236}")
        self._write_a16_creation_id()
        self._xml_end_tag("a:ext")

        self._write_uri_ext("{C183D7F6-B498-43B3-948B-1728B52AA6E4}")
        self._write_adec_decorative()
        self._xml_end_tag("a:ext")

        self._xml_end_tag("a:extLst")

    def _write_uri_ext(self, uri) -> None:
        # Write the <a:ext> element.
        attributes = [("uri", uri)]

        self._xml_start_tag("a:ext", attributes)

    def _write_adec_decorative(self) -> None:
        # Write the <adec:decorative> element.
        xmlns = "http://schemas.microsoft.com/office/drawing/2017/decorative"
        val = "1"

        attributes = [
            ("xmlns:adec", xmlns),
            ("val", val),
        ]

        self._xml_empty_tag("adec:decorative", attributes)

    def _write_a16_creation_id(self) -> None:
        # Write the <a16:creationId> element.

        xmlns_a_16 = "http://schemas.microsoft.com/office/drawing/2014/main"
        creation_id = "{00000000-0008-0000-0000-000002000000}"

        attributes = [
            ("xmlns:a16", xmlns_a_16),
            ("id", creation_id),
        ]

        self._xml_empty_tag("a16:creationId", attributes)

    def _write_a_hlink_click(self, url: Url) -> None:
        # Write the <a:hlinkClick> element.
        schema = "http://schemas.openxmlformats.org/officeDocument/"
        xmlns_r = schema + "2006/relationships"

        attributes = [
            ("xmlns:r", xmlns_r),
            ("r:id", "rId" + str(url._rel_index)),
        ]

        if url._tip:
            attributes.append(("tooltip", url._tip))

        self._xml_empty_tag("a:hlinkClick", attributes)

    def _write_c_nv_graphic_frame_pr(self) -> None:
        # Write the <xdr:cNvGraphicFramePr> element.
        if self.embedded:
            self._xml_empty_tag("xdr:cNvGraphicFramePr")
        else:
            self._xml_start_tag("xdr:cNvGraphicFramePr")

            # Write the a:graphicFrameLocks element.
            self._write_a_graphic_frame_locks()

            self._xml_end_tag("xdr:cNvGraphicFramePr")

    def _write_a_graphic_frame_locks(self) -> None:
        # Write the <a:graphicFrameLocks> element.
        attributes = [("noGrp", 1)]

        self._xml_empty_tag("a:graphicFrameLocks", attributes)

    def _write_xfrm(self) -> None:
        # Write the <xdr:xfrm> element.
        self._xml_start_tag("xdr:xfrm")

        # Write the xfrmOffset element.
        self._write_xfrm_offset()

        # Write the xfrmOffset element.
        self._write_xfrm_extension()

        self._xml_end_tag("xdr:xfrm")

    def _write_xfrm_offset(self) -> None:
        # Write the <a:off> xfrm sub-element.

        attributes = [
            ("x", 0),
            ("y", 0),
        ]

        self._xml_empty_tag("a:off", attributes)

    def _write_xfrm_extension(self) -> None:
        # Write the <a:ext> xfrm sub-element.

        attributes = [
            ("cx", 0),
            ("cy", 0),
        ]

        self._xml_empty_tag("a:ext", attributes)

    def _write_atag_graphic(self, index: int) -> None:
        # Write the <a:graphic> element.
        self._xml_start_tag("a:graphic")

        # Write the a:graphicData element.
        self._write_atag_graphic_data(index)

        self._xml_end_tag("a:graphic")

    def _write_atag_graphic_data(self, index: int) -> None:
        # Write the <a:graphicData> element.
        uri = "http://schemas.openxmlformats.org/drawingml/2006/chart"

        attributes = [
            (
                "uri",
                uri,
            )
        ]

        self._xml_start_tag("a:graphicData", attributes)

        # Write the c:chart element.
        self._write_c_chart("rId" + str(index))

        self._xml_end_tag("a:graphicData")

    def _write_c_chart(self, r_id) -> None:
        # Write the <c:chart> element.

        schema = "http://schemas.openxmlformats.org/"
        xmlns_c = schema + "drawingml/2006/chart"
        xmlns_r = schema + "officeDocument/2006/relationships"

        attributes = [
            ("xmlns:c", xmlns_c),
            ("xmlns:r", xmlns_r),
            ("r:id", r_id),
        ]

        self._xml_empty_tag("c:chart", attributes)

    def _write_client_data(self) -> None:
        # Write the <xdr:clientData> element.
        self._xml_empty_tag("xdr:clientData")

    def _write_sp(
        self,
        index,
        col_absolute,
        row_absolute,
        drawing: DrawingInfo,
    ) -> None:
        # Write the <xdr:sp> element.

        if drawing._shape and drawing._shape.connect:
            attributes = [("macro", "")]
            self._xml_start_tag("xdr:cxnSp", attributes)

            # Write the xdr:nvCxnSpPr element.
            self._write_nv_cxn_sp_pr(drawing._shape)

            # Write the xdr:spPr element.
            self._write_xdr_sp_pr(col_absolute, row_absolute, drawing)

            self._xml_end_tag("xdr:cxnSp")
        else:
            # Add attribute for shapes.
            attributes = [("macro", ""), ("textlink", drawing._shape.textlink)]

            self._xml_start_tag("xdr:sp", attributes)

            # Write the xdr:nvSpPr element.
            self._write_nv_sp_pr(index, drawing)

            # Write the xdr:spPr element.
            self._write_xdr_sp_pr(col_absolute, row_absolute, drawing)

            # Write the xdr:style element.
            self._write_style()

            # Write the xdr:txBody element.
            if drawing._shape.text is not None:
                self._write_tx_body(drawing._shape)

            self._xml_end_tag("xdr:sp")

    def _write_nv_cxn_sp_pr(self, shape) -> None:
        # Write the <xdr:nvCxnSpPr> element.
        self._xml_start_tag("xdr:nvCxnSpPr")

        self._xml_start_tag("xdr:cNvCxnSpPr")

        attributes = [("noChangeShapeType", "1")]
        self._xml_empty_tag("a:cxnSpLocks", attributes)

        if shape.start:
            attributes = [("id", shape.start), ("idx", shape.start_index)]
            self._xml_empty_tag("a:stCxn", attributes)

        if shape.end:
            attributes = [("id", shape.end), ("idx", shape.end_index)]
            self._xml_empty_tag("a:endCxn", attributes)

        self._xml_end_tag("xdr:cNvCxnSpPr")
        self._xml_end_tag("xdr:nvCxnSpPr")

    def _write_nv_sp_pr(self, index: int, drawing: DrawingInfo) -> None:
        # Write the <xdr:NvSpPr> element.
        attributes = []

        self._xml_start_tag("xdr:nvSpPr")

        name = drawing._shape.name + " " + str(index)

        self._write_c_nv_pr(index + 1, drawing, name)

        if drawing._shape.name == "TextBox":
            attributes = [("txBox", 1)]

        self._xml_empty_tag("xdr:cNvSpPr", attributes)

        self._xml_end_tag("xdr:nvSpPr")

    def _write_pic(
        self,
        index: int,
        col_absolute: int,
        row_absolute: int,
        drawing: DrawingInfo,
    ) -> None:
        # Write the <xdr:pic> element.
        self._xml_start_tag("xdr:pic")

        # Write the xdr:nvPicPr element.
        self._write_nv_pic_pr(index, drawing)
        # Write the xdr:blipFill element.
        self._write_blip_fill(drawing._rel_index)

        # Write the xdr:spPr element.
        self._write_sp_pr(col_absolute, row_absolute, drawing)

        self._xml_end_tag("xdr:pic")

    def _write_nv_pic_pr(self, index: int, drawing: DrawingInfo) -> None:
        # Write the <xdr:nvPicPr> element.
        self._xml_start_tag("xdr:nvPicPr")

        name = "Picture " + str(index)

        # Write the xdr:cNvPr element.
        self._write_c_nv_pr(index + 1, drawing, name)

        # Write the xdr:cNvPicPr element.
        self._write_c_nv_pic_pr()

        self._xml_end_tag("xdr:nvPicPr")

    def _write_c_nv_pic_pr(self) -> None:
        # Write the <xdr:cNvPicPr> element.
        self._xml_start_tag("xdr:cNvPicPr")

        # Write the a:picLocks element.
        self._write_a_pic_locks()

        self._xml_end_tag("xdr:cNvPicPr")

    def _write_a_pic_locks(self) -> None:
        # Write the <a:picLocks> element.
        attributes = [("noChangeAspect", 1)]

        self._xml_empty_tag("a:picLocks", attributes)

    def _write_blip_fill(self, index: int) -> None:
        # Write the <xdr:blipFill> element.
        self._xml_start_tag("xdr:blipFill")

        # Write the a:blip element.
        self._write_a_blip(index)

        # Write the a:stretch element.
        self._write_a_stretch()

        self._xml_end_tag("xdr:blipFill")

    def _write_a_blip(self, index: int) -> None:
        # Write the <a:blip> element.
        schema = "http://schemas.openxmlformats.org/officeDocument/"
        xmlns_r = schema + "2006/relationships"
        r_embed = "rId" + str(index)

        attributes = [("xmlns:r", xmlns_r), ("r:embed", r_embed)]

        self._xml_empty_tag("a:blip", attributes)

    def _write_a_stretch(self) -> None:
        # Write the <a:stretch> element.
        self._xml_start_tag("a:stretch")

        # Write the a:fillRect element.
        self._write_a_fill_rect()

        self._xml_end_tag("a:stretch")

    def _write_a_fill_rect(self) -> None:
        # Write the <a:fillRect> element.
        self._xml_empty_tag("a:fillRect")

    def _write_sp_pr(self, col_absolute, row_absolute, drawing: DrawingInfo) -> None:
        # Write the <xdr:spPr> element, for charts.

        self._xml_start_tag("xdr:spPr")

        # Write the a:xfrm element.
        self._write_a_xfrm(col_absolute, row_absolute, drawing._width, drawing._height)

        # Write the a:prstGeom element.
        self._write_a_prst_geom(drawing._shape)

        self._xml_end_tag("xdr:spPr")

    def _write_xdr_sp_pr(
        self, col_absolute: int, row_absolute: int, drawing: DrawingInfo
    ) -> None:
        # Write the <xdr:spPr> element for shapes.
        self._xml_start_tag("xdr:spPr")

        # Write the a:xfrm element.
        self._write_a_xfrm(
            col_absolute, row_absolute, drawing._width, drawing._height, drawing._shape
        )

        # Write the a:prstGeom element.
        shape = drawing._shape
        self._write_a_prst_geom(shape)

        if shape.fill:
            if not shape.fill["defined"]:
                # Write the a:solidFill element.
                self._write_a_solid_fill_scheme("lt1")
            elif "none" in shape.fill:
                # Write the a:noFill element.
                self._xml_empty_tag("a:noFill")
            elif "color" in shape.fill:
                # Write the a:solidFill element.
                self._write_a_solid_fill(shape.fill["color"])

        if shape.gradient:
            # Write the a:gradFill element.
            self._write_a_grad_fill(shape.gradient)

        # Write the a:ln element.
        self._write_a_ln(shape.line)

        self._xml_end_tag("xdr:spPr")

    def _write_a_xfrm(
        self, col_absolute, row_absolute, width, height, shape=None
    ) -> None:
        # Write the <a:xfrm> element.
        attributes = []

        if shape:
            if shape.rotation:
                rotation = shape.rotation
                rotation *= 60000
                attributes.append(("rot", rotation))

            if shape.flip_h:
                attributes.append(("flipH", 1))
            if shape.flip_v:
                attributes.append(("flipV", 1))

        self._xml_start_tag("a:xfrm", attributes)

        # Write the a:off element.
        self._write_a_off(col_absolute, row_absolute)

        # Write the a:ext element.
        self._write_a_ext(width, height)

        self._xml_end_tag("a:xfrm")

    def _write_a_off(self, x, y) -> None:
        # Write the <a:off> element.
        attributes = [
            ("x", x),
            ("y", y),
        ]

        self._xml_empty_tag("a:off", attributes)

    def _write_a_ext(self, cx, cy) -> None:
        # Write the <a:ext> element.
        attributes = [
            ("cx", cx),
            ("cy", cy),
        ]

        self._xml_empty_tag("a:ext", attributes)

    def _write_a_prst_geom(self, shape=None) -> None:
        # Write the <a:prstGeom> element.
        attributes = [("prst", "rect")]

        self._xml_start_tag("a:prstGeom", attributes)

        # Write the a:avLst element.
        self._write_a_av_lst(shape)

        self._xml_end_tag("a:prstGeom")

    def _write_a_av_lst(self, shape=None) -> None:
        # Write the <a:avLst> element.
        adjustments = []

        if shape and shape.adjustments:
            adjustments = shape.adjustments

        if adjustments:
            self._xml_start_tag("a:avLst")

            i = 0
            for adj in adjustments:
                i += 1
                # Only connectors have multiple adjustments.
                if shape.connect:
                    suffix = i
                else:
                    suffix = ""

                # Scale Adjustments: 100,000 = 100%.
                adj_int = str(int(adj * 1000))

                attributes = [("name", "adj" + suffix), ("fmla", "val" + adj_int)]

                self._xml_empty_tag("a:gd", attributes)

            self._xml_end_tag("a:avLst")
        else:
            self._xml_empty_tag("a:avLst")

    def _write_a_solid_fill(self, color: Color) -> None:
        # Write the <a:solidFill> element.
        self._xml_start_tag("a:solidFill")

        # Write the a:srgbClr element.
        self._write_a_srgb_clr(color)

        self._xml_end_tag("a:solidFill")

    def _write_a_solid_fill_scheme(self, named_color, shade=None) -> None:
        attributes = [("val", named_color)]

        self._xml_start_tag("a:solidFill")

        if shade:
            self._xml_start_tag("a:schemeClr", attributes)
            self._write_a_shade(shade)
            self._xml_end_tag("a:schemeClr")
        else:
            self._xml_empty_tag("a:schemeClr", attributes)

        self._xml_end_tag("a:solidFill")

    def _write_a_ln(self, line) -> None:
        # Write the <a:ln> element.
        width = line.get("width", 0.75)

        # Round width to nearest 0.25, like Excel.
        width = int((width + 0.125) * 4) / 4.0

        # Convert to internal units.
        width = int(0.5 + (12700 * width))

        attributes = [("w", width), ("cmpd", "sng")]

        self._xml_start_tag("a:ln", attributes)

        if "none" in line:
            # Write the a:noFill element.
            self._xml_empty_tag("a:noFill")

        elif "color" in line:
            # Write the a:solidFill element.
            self._write_a_solid_fill(line["color"])

        else:
            # Write the a:solidFill element.
            self._write_a_solid_fill_scheme("lt1", "50000")

        # Write the line/dash type.
        line_type = line.get("dash_type")
        if line_type:
            # Write the a:prstDash element.
            self._write_a_prst_dash(line_type)

        self._xml_end_tag("a:ln")

    def _write_tx_body(self, shape) -> None:
        # Write the <xdr:txBody> element.
        attributes = []

        if shape.text_rotation != 0:
            if shape.text_rotation == 90:
                attributes.append(("vert", "vert270"))
            if shape.text_rotation == -90:
                attributes.append(("vert", "vert"))
            if shape.text_rotation == 270:
                attributes.append(("vert", "wordArtVert"))
            if shape.text_rotation == 271:
                attributes.append(("vert", "eaVert"))

        attributes.append(("wrap", "square"))
        attributes.append(("rtlCol", "0"))

        if not shape.align["defined"]:
            attributes.append(("anchor", "t"))
        else:
            if "vertical" in shape.align:
                align = shape.align["vertical"]
                if align == "top":
                    attributes.append(("anchor", "t"))
                elif align == "middle":
                    attributes.append(("anchor", "ctr"))
                elif align == "bottom":
                    attributes.append(("anchor", "b"))
            else:
                attributes.append(("anchor", "t"))

            if "horizontal" in shape.align:
                align = shape.align["horizontal"]
                if align == "center":
                    attributes.append(("anchorCtr", "1"))
            else:
                attributes.append(("anchorCtr", "0"))

        self._xml_start_tag("xdr:txBody")
        self._xml_empty_tag("a:bodyPr", attributes)
        self._xml_empty_tag("a:lstStyle")

        lines = shape.text.split("\n")

        # Set the font attributes.
        font = shape.font
        # pylint: disable=protected-access
        style_attrs = Shape._get_font_style_attributes(font)
        latin_attrs = Shape._get_font_latin_attributes(font)
        style_attrs.insert(0, ("lang", font["lang"]))

        if shape.textlink != "":
            attributes = [
                ("id", "{B8ADDEFE-BF52-4FD4-8C5D-6B85EF6FF707}"),
                ("type", "TxLink"),
            ]

            self._xml_start_tag("a:p")
            self._xml_start_tag("a:fld", attributes)

            self._write_font_run(font, style_attrs, latin_attrs, "a:rPr")

            self._xml_data_element("a:t", shape.text)
            self._xml_end_tag("a:fld")

            self._write_font_run(font, style_attrs, latin_attrs, "a:endParaRPr")

            self._xml_end_tag("a:p")
        else:
            for line in lines:
                self._xml_start_tag("a:p")

                if line == "":
                    self._write_font_run(font, style_attrs, latin_attrs, "a:endParaRPr")
                    self._xml_end_tag("a:p")
                    continue

                if "text" in shape.align:
                    if shape.align["text"] == "left":
                        self._xml_empty_tag("a:pPr", [("algn", "l")])
                    if shape.align["text"] == "center":
                        self._xml_empty_tag("a:pPr", [("algn", "ctr")])
                    if shape.align["text"] == "right":
                        self._xml_empty_tag("a:pPr", [("algn", "r")])

                self._xml_start_tag("a:r")

                self._write_font_run(font, style_attrs, latin_attrs, "a:rPr")

                self._xml_data_element("a:t", line)

                self._xml_end_tag("a:r")
                self._xml_end_tag("a:p")

        self._xml_end_tag("xdr:txBody")

    def _write_font_run(self, font, style_attrs, latin_attrs, run_type) -> None:
        # Write a:rPr or a:endParaRPr.
        has_color = font.get("color") is not None

        if latin_attrs or has_color:
            self._xml_start_tag(run_type, style_attrs)

            if has_color:
                self._write_a_solid_fill(font["color"])

            if latin_attrs:
                self._write_a_latin(latin_attrs)
                self._write_a_cs(latin_attrs)

            self._xml_end_tag(run_type)
        else:
            self._xml_empty_tag(run_type, style_attrs)

    def _write_style(self) -> None:
        # Write the <xdr:style> element.
        self._xml_start_tag("xdr:style")

        # Write the a:lnRef element.
        self._write_a_ln_ref()

        # Write the a:fillRef element.
        self._write_a_fill_ref()

        # Write the a:effectRef element.
        self._write_a_effect_ref()

        # Write the a:fontRef element.
        self._write_a_font_ref()

        self._xml_end_tag("xdr:style")

    def _write_a_ln_ref(self) -> None:
        # Write the <a:lnRef> element.
        attributes = [("idx", "0")]

        self._xml_start_tag("a:lnRef", attributes)

        # Write the a:scrgbClr element.
        self._write_a_scrgb_clr()

        self._xml_end_tag("a:lnRef")

    def _write_a_fill_ref(self) -> None:
        # Write the <a:fillRef> element.
        attributes = [("idx", "0")]

        self._xml_start_tag("a:fillRef", attributes)

        # Write the a:scrgbClr element.
        self._write_a_scrgb_clr()

        self._xml_end_tag("a:fillRef")

    def _write_a_effect_ref(self) -> None:
        # Write the <a:effectRef> element.
        attributes = [("idx", "0")]

        self._xml_start_tag("a:effectRef", attributes)

        # Write the a:scrgbClr element.
        self._write_a_scrgb_clr()

        self._xml_end_tag("a:effectRef")

    def _write_a_scrgb_clr(self) -> None:
        # Write the <a:scrgbClr> element.

        attributes = [
            ("r", "0"),
            ("g", "0"),
            ("b", "0"),
        ]

        self._xml_empty_tag("a:scrgbClr", attributes)

    def _write_a_font_ref(self) -> None:
        # Write the <a:fontRef> element.
        attributes = [("idx", "minor")]

        self._xml_start_tag("a:fontRef", attributes)

        # Write the a:schemeClr element.
        self._write_a_scheme_clr("dk1")

        self._xml_end_tag("a:fontRef")

    def _write_a_scheme_clr(self, val) -> None:
        # Write the <a:schemeClr> element.
        attributes = [("val", val)]

        self._xml_empty_tag("a:schemeClr", attributes)

    def _write_a_shade(self, shade) -> None:
        # Write the <a:shade> element.
        attributes = [("val", shade)]

        self._xml_empty_tag("a:shade", attributes)

    def _write_a_prst_dash(self, val) -> None:
        # Write the <a:prstDash> element.

        attributes = [("val", val)]

        self._xml_empty_tag("a:prstDash", attributes)

    def _write_a_grad_fill(self, gradient) -> None:
        # Write the <a:gradFill> element.

        attributes = [("flip", "none"), ("rotWithShape", "1")]

        if gradient["type"] == "linear":
            attributes = []

        self._xml_start_tag("a:gradFill", attributes)

        # Write the a:gsLst element.
        self._write_a_gs_lst(gradient)

        if gradient["type"] == "linear":
            # Write the a:lin element.
            self._write_a_lin(gradient["angle"])
        else:
            # Write the a:path element.
            self._write_a_path(gradient["type"])

            # Write the a:tileRect element.
            self._write_a_tile_rect(gradient["type"])

        self._xml_end_tag("a:gradFill")

    def _write_a_gs_lst(self, gradient) -> None:
        # Write the <a:gsLst> element.
        positions = gradient["positions"]
        colors = gradient["colors"]

        self._xml_start_tag("a:gsLst")

        for i, color in enumerate(colors):
            pos = int(positions[i] * 1000)
            attributes = [("pos", pos)]
            self._xml_start_tag("a:gs", attributes)

            # Write the a:srgbClr element.
            self._write_a_srgb_clr(color)

            self._xml_end_tag("a:gs")

        self._xml_end_tag("a:gsLst")

    def _write_a_lin(self, angle) -> None:
        # Write the <a:lin> element.

        angle = int(60000 * angle)

        attributes = [
            ("ang", angle),
            ("scaled", "0"),
        ]

        self._xml_empty_tag("a:lin", attributes)

    def _write_a_path(self, gradient_type) -> None:
        # Write the <a:path> element.

        attributes = [("path", gradient_type)]

        self._xml_start_tag("a:path", attributes)

        # Write the a:fillToRect element.
        self._write_a_fill_to_rect(gradient_type)

        self._xml_end_tag("a:path")

    def _write_a_fill_to_rect(self, gradient_type) -> None:
        # Write the <a:fillToRect> element.

        if gradient_type == "shape":
            attributes = [
                ("l", "50000"),
                ("t", "50000"),
                ("r", "50000"),
                ("b", "50000"),
            ]
        else:
            attributes = [
                ("l", "100000"),
                ("t", "100000"),
            ]

        self._xml_empty_tag("a:fillToRect", attributes)

    def _write_a_tile_rect(self, gradient_type) -> None:
        # Write the <a:tileRect> element.

        if gradient_type == "shape":
            attributes = []
        else:
            attributes = [
                ("r", "-100000"),
                ("b", "-100000"),
            ]

        self._xml_empty_tag("a:tileRect", attributes)

    def _write_a_srgb_clr(self, color: Color) -> None:
        # Write the <a:srgbClr> element.
        attributes = [("val", color._rgb_hex_value())]

        self._xml_empty_tag("a:srgbClr", attributes)

    def _write_a_latin(self, attributes) -> None:
        # Write the <a:latin> element.
        self._xml_empty_tag("a:latin", attributes)

    def _write_a_cs(self, attributes) -> None:
        # Write the <a:latin> element.
        self._xml_empty_tag("a:cs", attributes)
