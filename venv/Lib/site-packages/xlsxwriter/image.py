###############################################################################
#
# Image - A class for representing image objects in Excel.
#
# SPDX-License-Identifier: BSD-2-Clause
#
# Copyright (c) 2013-2025, John McNamara, jmcnamara@cpan.org
#

import hashlib
import os
from io import BytesIO
from pathlib import Path
from struct import unpack
from typing import Tuple, Union

from xlsxwriter.url import Url

from .exceptions import UndefinedImageSize, UnsupportedImageFormat

DEFAULT_DPI = 96.0


class Image:
    """
    A class to represent an image in an Excel worksheet.

    """

    def __init__(self, source: Union[str, Path, BytesIO]) -> None:
        """
        Initialize an Image instance.

        Args:
            source (Union[str, Path, BytesIO]): The filename, Path or BytesIO
            object of the image.
        """
        if isinstance(source, (str, Path)):
            self.filename = source
            self.image_data = None
            self.image_name = os.path.basename(source)
        elif isinstance(source, BytesIO):
            self.filename = ""
            self.image_data = source
            self.image_name = ""
        else:
            raise ValueError("Source must be a filename (str) or a BytesIO object.")

        self._row: int = 0
        self._col: int = 0
        self._x_offset: int = 0
        self._y_offset: int = 0
        self._x_scale: float = 1.0
        self._y_scale: float = 1.0
        self._url: Union[Url, None] = None
        self._anchor: int = 2
        self._description: Union[str, None] = None
        self._decorative: bool = False
        self._header_position: Union[str, None] = None
        self._ref_id: Union[str, None] = None

        # Derived properties.
        self._image_extension: str = ""
        self._width: float = 0.0
        self._height: float = 0.0
        self._x_dpi: float = DEFAULT_DPI
        self._y_dpi: float = DEFAULT_DPI
        self._digest: Union[str, None] = None

        self._get_image_properties()

    def __repr__(self) -> str:
        """
        Return a string representation of the main properties of the Image
        instance.
        """
        return (
            f"Image:\n"
            f"    filename   = {self.filename!r}\n"
            f"    image_name = {self.image_name!r}\n"
            f"    image_type = {self.image_type!r}\n"
            f"    width      = {self._width}\n"
            f"    height     = {self._height}\n"
            f"    x_dpi      = {self._x_dpi}\n"
            f"    y_dpi      = {self._y_dpi}\n"
        )

    @property
    def image_type(self) -> str:
        """Get the image type (e.g., 'PNG', 'JPEG')."""
        return self._image_extension.upper()

    @property
    def width(self) -> float:
        """Get the width of the image."""
        return self._width

    @property
    def height(self) -> float:
        """Get the height of the image."""
        return self._height

    @property
    def x_dpi(self) -> float:
        """Get the horizontal DPI of the image."""
        return self._x_dpi

    @property
    def y_dpi(self) -> float:
        """Get the vertical DPI of the image."""
        return self._y_dpi

    @property
    def description(self) -> Union[str, None]:
        """Get the description/alt-text of the image."""
        return self._description

    @description.setter
    def description(self, value: str) -> None:
        """Set the description/alt-text of the image."""
        if value:
            self._description = value

    @property
    def decorative(self) -> bool:
        """Get whether the image is decorative."""
        return self._decorative

    @decorative.setter
    def decorative(self, value: bool) -> None:
        """Set whether the image is decorative."""
        self._decorative = value

    @property
    def url(self) -> Union[Url, None]:
        """Get the image url."""
        return self._url

    @url.setter
    def url(self, value: Url) -> None:
        """Set the image url."""
        if value:
            self._url = value

    def _set_user_options(self, options=None) -> None:
        """
        This handles the additional optional parameters to ``insert_button()``.
        """
        if options is None:
            return

        if not self._url:
            self._url = Url.from_options(options)
            if self._url:
                self._url._set_object_link()

        self._anchor = options.get("object_position", self._anchor)
        self._x_scale = options.get("x_scale", self._x_scale)
        self._y_scale = options.get("y_scale", self._y_scale)
        self._x_offset = options.get("x_offset", self._x_offset)
        self._y_offset = options.get("y_offset", self._y_offset)
        self._decorative = options.get("decorative", self._decorative)
        self.image_data = options.get("image_data", self.image_data)
        self._description = options.get("description", self._description)

        # For backward compatibility with older parameter name.
        self._anchor = options.get("positioning", self._anchor)

    def _get_image_properties(self) -> None:
        # Extract dimension information from the image file.
        height = 0.0
        width = 0.0
        x_dpi = DEFAULT_DPI
        y_dpi = DEFAULT_DPI

        if self.image_data:
            # Read the image data from the user supplied byte stream.
            data = self.image_data.getvalue()
        else:
            # Open the image file and read in the data.
            with open(self.filename, "rb") as fh:
                data = fh.read()

        # Get the image digest to check for duplicates.
        digest = hashlib.sha256(data).hexdigest()

        # Look for some common image file markers.
        png_marker = unpack("3s", data[1:4])[0]
        jpg_marker = unpack(">H", data[:2])[0]
        bmp_marker = unpack("2s", data[:2])[0]
        gif_marker = unpack("4s", data[:4])[0]
        emf_marker = (unpack("4s", data[40:44]))[0]
        emf_marker1 = unpack("<L", data[:4])[0]

        if png_marker == b"PNG":
            (image_type, width, height, x_dpi, y_dpi) = self._process_png(data)

        elif jpg_marker == 0xFFD8:
            (image_type, width, height, x_dpi, y_dpi) = self._process_jpg(data)

        elif bmp_marker == b"BM":
            (image_type, width, height) = self._process_bmp(data)

        elif emf_marker1 == 0x9AC6CDD7:
            (image_type, width, height, x_dpi, y_dpi) = self._process_wmf(data)

        elif emf_marker1 == 1 and emf_marker == b" EMF":
            (image_type, width, height, x_dpi, y_dpi) = self._process_emf(data)

        elif gif_marker == b"GIF8":
            (image_type, width, height, x_dpi, y_dpi) = self._process_gif(data)

        else:
            raise UnsupportedImageFormat(
                f"{self.filename}: Unknown or unsupported image file format."
            )

        # Check that we found the required data.
        if not height or not width:
            raise UndefinedImageSize(
                f"{self.filename}: no size data found in image file."
            )

        # Set a default dpi for images with 0 dpi.
        if x_dpi == 0:
            x_dpi = DEFAULT_DPI
        if y_dpi == 0:
            y_dpi = DEFAULT_DPI

        self._image_extension = image_type
        self._width = width
        self._height = height
        self._x_dpi = x_dpi
        self._y_dpi = y_dpi
        self._digest = digest

    def _process_png(
        self,
        data: bytes,
    ) -> Tuple[str, float, float, float, float]:
        # Extract width and height information from a PNG file.
        offset = 8
        data_length = len(data)
        end_marker = False
        width = 0.0
        height = 0.0
        x_dpi = DEFAULT_DPI
        y_dpi = DEFAULT_DPI

        # Search through the image data to read the height and width in the
        # IHDR element. Also read the DPI in the pHYs element.
        while not end_marker and offset < data_length:
            length = unpack(">I", data[offset + 0 : offset + 4])[0]
            marker = unpack("4s", data[offset + 4 : offset + 8])[0]

            # Read the image dimensions.
            if marker == b"IHDR":
                width = unpack(">I", data[offset + 8 : offset + 12])[0]
                height = unpack(">I", data[offset + 12 : offset + 16])[0]

            # Read the image DPI.
            if marker == b"pHYs":
                x_density = unpack(">I", data[offset + 8 : offset + 12])[0]
                y_density = unpack(">I", data[offset + 12 : offset + 16])[0]
                units = unpack("b", data[offset + 16 : offset + 17])[0]

                if units == 1 and x_density > 0 and y_density > 0:
                    x_dpi = x_density * 0.0254
                    y_dpi = y_density * 0.0254

            if marker == b"IEND":
                end_marker = True
                continue

            offset = offset + length + 12

        return "png", width, height, x_dpi, y_dpi

    def _process_jpg(self, data: bytes) -> Tuple[str, float, float, float, float]:
        # Extract width and height information from a JPEG file.
        offset = 2
        data_length = len(data)
        end_marker = False
        width = 0.0
        height = 0.0
        x_dpi = DEFAULT_DPI
        y_dpi = DEFAULT_DPI

        # Search through the image data to read the JPEG markers.
        while not end_marker and offset < data_length:
            marker = unpack(">H", data[offset + 0 : offset + 2])[0]
            length = unpack(">H", data[offset + 2 : offset + 4])[0]

            # Read the height and width in the 0xFFCn elements (except C4, C8
            # and CC which aren't SOF markers).
            if (
                (marker & 0xFFF0) == 0xFFC0
                and marker != 0xFFC4
                and marker != 0xFFC8
                and marker != 0xFFCC
            ):
                height = unpack(">H", data[offset + 5 : offset + 7])[0]
                width = unpack(">H", data[offset + 7 : offset + 9])[0]

            # Read the DPI in the 0xFFE0 element.
            if marker == 0xFFE0:
                units = unpack("b", data[offset + 11 : offset + 12])[0]
                x_density = unpack(">H", data[offset + 12 : offset + 14])[0]
                y_density = unpack(">H", data[offset + 14 : offset + 16])[0]

                if units == 1:
                    x_dpi = x_density
                    y_dpi = y_density

                if units == 2:
                    x_dpi = x_density * 2.54
                    y_dpi = y_density * 2.54

                # Workaround for incorrect dpi.
                if x_dpi == 1:
                    x_dpi = DEFAULT_DPI
                if y_dpi == 1:
                    y_dpi = DEFAULT_DPI

            if marker == 0xFFDA:
                end_marker = True
                continue

            offset = offset + length + 2

        return "jpeg", width, height, x_dpi, y_dpi

    def _process_gif(self, data: bytes) -> Tuple[str, float, float, float, float]:
        # Extract width and height information from a GIF file.
        x_dpi = DEFAULT_DPI
        y_dpi = DEFAULT_DPI

        width = unpack("<h", data[6:8])[0]
        height = unpack("<h", data[8:10])[0]

        return "gif", width, height, x_dpi, y_dpi

    def _process_bmp(self, data: bytes) -> Tuple[str, float, float]:
        # Extract width and height information from a BMP file.
        width = unpack("<L", data[18:22])[0]
        height = unpack("<L", data[22:26])[0]
        return "bmp", width, height

    def _process_wmf(self, data: bytes) -> Tuple[str, float, float, float, float]:
        # Extract width and height information from a WMF file.
        x_dpi = DEFAULT_DPI
        y_dpi = DEFAULT_DPI

        # Read the bounding box, measured in logical units.
        x1 = unpack("<h", data[6:8])[0]
        y1 = unpack("<h", data[8:10])[0]
        x2 = unpack("<h", data[10:12])[0]
        y2 = unpack("<h", data[12:14])[0]

        # Read the number of logical units per inch. Used to scale the image.
        inch = unpack("<H", data[14:16])[0]

        # Convert to rendered height and width.
        width = float((x2 - x1) * x_dpi) / inch
        height = float((y2 - y1) * y_dpi) / inch

        return "wmf", width, height, x_dpi, y_dpi

    def _process_emf(self, data: bytes) -> Tuple[str, float, float, float, float]:
        # Extract width and height information from a EMF file.

        # Read the bounding box, measured in logical units.
        bound_x1 = unpack("<l", data[8:12])[0]
        bound_y1 = unpack("<l", data[12:16])[0]
        bound_x2 = unpack("<l", data[16:20])[0]
        bound_y2 = unpack("<l", data[20:24])[0]

        # Convert the bounds to width and height.
        width = bound_x2 - bound_x1
        height = bound_y2 - bound_y1

        # Read the rectangular frame in units of 0.01mm.
        frame_x1 = unpack("<l", data[24:28])[0]
        frame_y1 = unpack("<l", data[28:32])[0]
        frame_x2 = unpack("<l", data[32:36])[0]
        frame_y2 = unpack("<l", data[36:40])[0]

        # Convert the frame bounds to mm width and height.
        width_mm = 0.01 * (frame_x2 - frame_x1)
        height_mm = 0.01 * (frame_y2 - frame_y1)

        # Get the dpi based on the logical size.
        x_dpi = width * 25.4 / width_mm
        y_dpi = height * 25.4 / height_mm

        # This is to match Excel's calculation. It is probably to account for
        # the fact that the bounding box is inclusive-inclusive. Or a bug.
        width += 1
        height += 1

        return "emf", width, height, x_dpi, y_dpi
