import io
import os
import re
import struct
from decimal import Decimal
from typing import BinaryIO, NamedTuple, Protocol, Tuple, Union, runtime_checkable
from urllib.parse import urlparse
from urllib.request import urlopen

from xml.etree import ElementTree

_UNIT_KM = -3
_UNIT_100M = -2
_UNIT_10M = -1
_UNIT_1M = 0
_UNIT_10CM = 1
_UNIT_CM = 2
_UNIT_MM = 3
_UNIT_0_1MM = 4
_UNIT_0_01MM = 5
_UNIT_UM = 6
_UNIT_INCH = 6

_TIFF_TYPE_SIZES = {
  1: 1,
  2: 1,
  3: 2,
  4: 4,
  5: 8,
  6: 1,
  7: 1,
  8: 2,
  9: 4,
  10: 8,
  11: 4,
  12: 8,
}

_HEIF_BRANDS = {
    b'avif', b'avis',
    b'heic', b'heix', b'hevc', b'hevx',
    b'mif1', b'msf1',
}

_HEIF_IROT_TO_EXIF = {
    0: 1,
    1: 6,
    2: 3,
    3: 8,
}

_JPEG_NO_SOF_MARKERS = {0xc4, 0xc8, 0xcc}


@runtime_checkable
class ReadSeekBinary(Protocol):
    def read(self, size: int = -1) -> bytes:
        ...

    def seek(self, offset: int, whence: int = 0) -> int:
        ...


PathInput = Union[str, bytes, os.PathLike]
FileInput = Union[PathInput, BinaryIO, ReadSeekBinary]


class ImageInfo(NamedTuple):
    width: int = -1
    height: int = -1
    rotation: int = -1
    xdpi: int = -1
    ydpi: int = -1
    colors: int = -1
    channels: int = -1


def _open_file(filepath: FileInput):
    if isinstance(filepath, ReadSeekBinary):
        return filepath, False
    if isinstance(filepath, str):
        parsed = urlparse(filepath)
        if parsed.scheme in ("http", "https"):
            with urlopen(filepath) as response:
                return io.BytesIO(response.read()), True
    return open(filepath, 'rb'), True


def _convertToDPI(density, unit):
    if unit == _UNIT_KM:
        return int(density * 0.0000254 + 0.5)
    elif unit == _UNIT_100M:
        return int(density * 0.000254 + 0.5)
    elif unit == _UNIT_10M:
        return int(density * 0.00254 + 0.5)
    elif unit == _UNIT_1M:
        return int(density * 0.0254 + 0.5)
    elif unit == _UNIT_10CM:
        return int(density * 0.254 + 0.5)
    elif unit == _UNIT_CM:
        return int(density * 2.54 + 0.5)
    elif unit == _UNIT_MM:
        return int(density * 25.4 + 0.5)
    elif unit == _UNIT_0_1MM:
        return density * 254
    elif unit == _UNIT_0_01MM:
        return density * 2540
    elif unit == _UNIT_UM:
        return density * 25400
    return density


def _convertToPx(value):
    matched = re.match(r"(\d+(?:\.\d+)?)?([a-z]*)$", value)
    if not matched:
        raise ValueError("unknown length value: %s" % value)

    length, unit = matched.groups()
    length = Decimal(length)
    if unit == "":
        return float(length)
    elif unit == "cm":
        return float(length * Decimal("96") / Decimal("2.54"))
    elif unit == "mm":
        return float(length * Decimal("96") / Decimal("25.4"))
    elif unit == "in":
        return float(length * Decimal("96"))
    elif unit == "pc":
        return float(length * Decimal("96") / Decimal("6"))
    elif unit == "pt":
        return float(length * Decimal("96") / Decimal("72"))
    elif unit == "px":
        return float(length)

    raise ValueError("unknown unit type: %s" % unit)


def _get_size(fhandle):
    height = -1
    width = -1
    fhandle.seek(0)
    head = fhandle.read(64)
    size = len(head)
    # handle GIFs
    if size >= 10 and head[:6] in (b'GIF87a', b'GIF89a'):
        # Check to see if content_type is correct
        try:
            width, height = struct.unpack("<hh", head[6:10])
        except struct.error:
            raise ValueError("Invalid GIF file")
    # see png edition spec bytes are below chunk length then and finally the
    elif size >= 24 and head.startswith(b'\211PNG\r\n\032\n') and head[12:16] == b'IHDR':
        try:
            width, height = struct.unpack(">LL", head[16:24])
        except struct.error:
            raise ValueError("Invalid PNG file")
    # Maybe this is for an older PNG version.
    elif size >= 16 and head.startswith(b'\211PNG\r\n\032\n'):
        # Check to see if we have the right content type
        try:
            width, height = struct.unpack(">LL", head[8:16])
        except struct.error:
            raise ValueError("Invalid PNG file")
    # handle JPEGs
    elif size >= 2 and head.startswith(b'\377\330'):
        try:
            fhandle.seek(0)
            _seek_to_jpeg_sof(fhandle)
            # We are at a SOFn block
            fhandle.seek(1, 1)  # Skip `precision' byte.
            height, width = struct.unpack('>HH', fhandle.read(4))
        except (struct.error, ValueError):
            raise ValueError("Invalid JPEG file")
    # handle JPEG2000s
    elif size >= 12 and head.startswith(b'\x00\x00\x00\x0cjP  \r\n\x87\n'):
        fhandle.seek(48)
        try:
            height, width = struct.unpack('>LL', fhandle.read(8))
        except struct.error:
            raise ValueError("Invalid JPEG2000 file")
    # handle AVIF/HEIF
    elif size >= 16 and head[4:8] == b'ftyp':
        ftyp_size = struct.unpack('>L', head[:4])[0]
        if ftyp_size < 8:
            raise ValueError("Invalid HEIF file")
        fhandle.seek(8)
        ftyp_payload = fhandle.read(ftyp_size - 8)
        if any(brand in ftyp_payload for brand in _HEIF_BRANDS):
            width, height, _, _ = _read_heif_metadata(fhandle)
            if width != -1 and height != -1:
                return width, height
            raise ValueError("Invalid HEIF file")
    # handle big endian TIFF
    elif size >= 8 and head.startswith(b"\x4d\x4d\x00\x2a"):
        offset = struct.unpack('>L', head[4:8])[0]
        fhandle.seek(offset)
        ifdsize = struct.unpack(">H", fhandle.read(2))[0]
        for i in range(ifdsize):
            tag, datatype, count, data = struct.unpack(">HHLL", fhandle.read(12))
            if tag == 256:
                if datatype == 3:
                    width = int(data / 65536)
                elif datatype == 4:
                    width = data
                else:
                    raise ValueError("Invalid TIFF file: width column data type should be SHORT/LONG.")
            elif tag == 257:
                if datatype == 3:
                    height = int(data / 65536)
                elif datatype == 4:
                    height = data
                else:
                    raise ValueError("Invalid TIFF file: height column data type should be SHORT/LONG.")
            if width != -1 and height != -1:
                break
        if width == -1 or height == -1:
            raise ValueError("Invalid TIFF file: width and/or height IDS entries are missing.")
    elif size >= 8 and head.startswith(b"\x49\x49\x2a\x00"):
        offset = struct.unpack('<L', head[4:8])[0]
        fhandle.seek(offset)
        ifdsize = struct.unpack("<H", fhandle.read(2))[0]
        for i in range(ifdsize):
            tag, datatype, count, data = struct.unpack("<HHLL", fhandle.read(12))
            if tag == 256:
                width = data
            elif tag == 257:
                height = data
            if width != -1 and height != -1:
                break
        if width == -1 or height == -1:
            raise ValueError("Invalid TIFF file: width and/or height IDS entries are missing.")
    # handle little endian BigTiff
    elif size >= 8 and head.startswith(b"\x49\x49\x2b\x00"):
        bytesize_offset = struct.unpack('<L', head[4:8])[0]
        if bytesize_offset != 8:
            raise ValueError('Invalid BigTIFF file: Expected offset to be 8, found {} instead.'.format(offset))
        offset = struct.unpack('<Q', head[8:16])[0]
        fhandle.seek(offset)
        ifdsize = struct.unpack("<Q", fhandle.read(8))[0]
        for i in range(ifdsize):
            tag, datatype, count, data = struct.unpack("<HHQQ", fhandle.read(20))
            if tag == 256:
                width = data
            elif tag == 257:
                height = data
            if width != -1 and height != -1:
                break
        if width == -1 or height == -1:
            raise ValueError("Invalid BigTIFF file: width and/or height IDS entries are missing.")

    # handle SVGs
    elif size >= 5 and (head.startswith(b'<?xml') or head.startswith(b'<svg')):
        fhandle.seek(0)
        data = fhandle.read(1024)
        try:
            data = data.decode('utf-8')
            width = re.search(r'[^-]width="(.*?)"', data).group(1)
            height = re.search(r'[^-]height="(.*?)"', data).group(1)
        except Exception:
            raise ValueError("Invalid SVG file")
        width = _convertToPx(width)
        height = _convertToPx(height)

    # handle Netpbm
    elif head[:1] == b"P" and head[1:2] in b"123456":
        fhandle.seek(2)
        sizes = []

        while True:
            next_chr = fhandle.read(1)

            if next_chr.isspace():
                continue

            if next_chr == b"":
                raise ValueError("Invalid Netpbm file")

            if next_chr == b"#":
                fhandle.readline()
                continue

            if not next_chr.isdigit():
                raise ValueError("Invalid character found on Netpbm file")

            size = next_chr
            next_chr = fhandle.read(1)

            while next_chr.isdigit():
                size += next_chr
                next_chr = fhandle.read(1)

            sizes.append(int(size))

            if len(sizes) == 2:
                break

            fhandle.seek(-1, os.SEEK_CUR)
        width, height = sizes
    elif head.startswith(b"RIFF") and head[8:12] == b"WEBP":
        if head[12:16] == b"VP8 ":
            width, height = struct.unpack("<HH", head[26:30])
        elif head[12:16] == b"VP8X":
            width = struct.unpack("<I", head[24:27] + b"\0")[0] + 1
            height = struct.unpack("<I", head[27:30] + b"\0")[0] + 1
        elif head[12:16] == b"VP8L":
            b = head[21:25]
            width = (((b[1] & 63) << 8) | b[0]) + 1
            height = (((b[3] & 15) << 10) | (b[2] << 2) | ((b[1] & 192) >> 6)) + 1
        else:
            raise ValueError("Unsupported WebP file")
    elif head.startswith(b'BM'):
        width, height = struct.unpack("<ll", head[18:26])

    return width, height


def _read_jpeg_exif_rotation(fhandle):
    fhandle.seek(0)
    head = fhandle.read(2)
    if not head.startswith(b'\377\330'):
        return -1

    while True:
        marker_start = fhandle.read(1)
        if not marker_start:
            break
        while marker_start == b'\xff':
            marker_code = fhandle.read(1)
            if marker_code != b'\xff':
                break
        else:
            continue

        if not marker_code or marker_code in (b'\xd9', b'\xda'):
            break

        try:
            segment_size = struct.unpack('>H', fhandle.read(2))[0]
        except struct.error:
            break
        if segment_size < 2:
            break

        payload = fhandle.read(segment_size - 2)
        if marker_code != b'\xe1' or not payload.startswith(b'Exif\x00\x00'):
            continue

        return _read_orientation_from_exif_payload(payload[6:])

    return -1


def _read_jpeg_segment_header(fhandle):
    marker_byte = fhandle.read(1)
    while marker_byte == b'\xff':
        marker_byte = fhandle.read(1)
    if not marker_byte:
        raise ValueError("Unexpected end of JPEG file")
    marker = marker_byte[0]
    segment_size = struct.unpack('>H', fhandle.read(2))[0] - 2
    if segment_size < 0:
        raise ValueError("Invalid JPEG segment size")
    return marker, segment_size


def _seek_to_jpeg_sof(fhandle):
    block_size = 2
    marker = 0
    while not (0xc0 <= marker <= 0xcf and marker not in _JPEG_NO_SOF_MARKERS):
        fhandle.seek(block_size, 1)
        marker, block_size = _read_jpeg_segment_header(fhandle)
    return marker, block_size


def _read_orientation_from_exif_payload(exif_data):
    if len(exif_data) < 8:
        return -1
    endian_token = exif_data[:2]
    if endian_token == b'II':
        endian = '<'
    elif endian_token == b'MM':
        endian = '>'
    else:
        return -1

    try:
        first_ifd_offset = struct.unpack(endian + 'L', exif_data[4:8])[0]
    except struct.error:
        return -1
    if first_ifd_offset + 2 > len(exif_data):
        return -1

    try:
        ifd_count = struct.unpack(endian + 'H', exif_data[first_ifd_offset:first_ifd_offset + 2])[0]
    except struct.error:
        return -1
    cursor = first_ifd_offset + 2

    for _ in range(ifd_count):
        if cursor + 12 > len(exif_data):
            return -1
        try:
            tag, datatype, count, value = struct.unpack(endian + 'HHLL', exif_data[cursor:cursor + 12])
        except struct.error:
            return -1
        if tag == 0x0112 and datatype == 3 and count == 1:
            return int(value / 65536) if endian == '>' else value & 0xFFFF
        cursor += 12
    return -1


def _read_heif_exif_rotation(fhandle):
    _, _, property_rotation, exif_rotation = _read_heif_metadata(fhandle)
    if property_rotation != -1:
        return property_rotation
    if exif_rotation != -1:
        return exif_rotation

    fhandle.seek(0)
    data = fhandle.read()
    marker = b'Exif\x00\x00'
    start = data.find(marker)
    if start == -1:
        return -1
    return _read_orientation_from_exif_payload(data[start + len(marker):])


def _iter_iso_boxes(data, start, end):
    offset = start
    while offset + 8 <= end:
        size = struct.unpack('>L', data[offset:offset + 4])[0]
        box_type = data[offset + 4:offset + 8]
        header_size = 8
        if size == 1:
            if offset + 16 > end:
                return
            size = struct.unpack('>Q', data[offset + 8:offset + 16])[0]
            header_size = 16
        elif size == 0:
            size = end - offset
        if size < header_size or offset + size > end:
            return
        yield offset, size, box_type, header_size
        offset += size


def _read_heif_metadata(fhandle):
    fhandle.seek(0)
    data = fhandle.read()

    meta_box = None
    for offset, size, box_type, header_size in _iter_iso_boxes(data, 0, len(data)):
        if box_type == b'meta':
            meta_box = (offset, size, header_size)
            break
    if meta_box is None:
        return -1, -1, -1, -1

    meta_offset, meta_size, meta_header = meta_box
    meta_start = meta_offset + meta_header + 4
    meta_end = meta_offset + meta_size

    primary_item_id = None
    properties = []
    associations = {}
    item_types = {}
    item_extents = {}

    for offset, size, box_type, header_size in _iter_iso_boxes(data, meta_start, meta_end):
        payload_start = offset + header_size
        payload_end = offset + size
        if box_type == b'pitm':
            version = data[payload_start]
            if version == 0 and payload_start + 6 <= payload_end:
                primary_item_id = struct.unpack('>H', data[payload_start + 4:payload_start + 6])[0]
            elif version > 0 and payload_start + 8 <= payload_end:
                primary_item_id = struct.unpack('>L', data[payload_start + 4:payload_start + 8])[0]
        elif box_type == b'iinf' and payload_start + 6 <= payload_end:
            version = data[payload_start]
            if version == 0:
                entry_count = struct.unpack('>H', data[payload_start + 4:payload_start + 6])[0]
                cursor = payload_start + 6
            else:
                if payload_start + 8 > payload_end:
                    continue
                entry_count = struct.unpack('>L', data[payload_start + 4:payload_start + 8])[0]
                cursor = payload_start + 8

            for _ in range(entry_count):
                if cursor + 8 > payload_end:
                    break
                entry_size = struct.unpack('>L', data[cursor:cursor + 4])[0]
                entry_type = data[cursor + 4:cursor + 8]
                entry_end = cursor + entry_size
                if entry_size < 8 or entry_end > payload_end:
                    break
                if entry_type == b'infe' and cursor + 13 <= payload_end:
                    infe_payload = cursor + 8
                    infe_version = data[infe_payload]
                    if infe_version == 2 and infe_payload + 12 <= entry_end:
                        item_id = struct.unpack('>H', data[infe_payload + 4:infe_payload + 6])[0]
                        item_type = data[infe_payload + 8:infe_payload + 12]
                        item_types[item_id] = item_type
                    elif infe_version >= 3 and infe_payload + 16 <= entry_end:
                        item_id = struct.unpack('>L', data[infe_payload + 4:infe_payload + 8])[0]
                        item_type = data[infe_payload + 12:infe_payload + 16]
                        item_types[item_id] = item_type
                cursor = entry_end
        elif box_type == b'iloc' and payload_start + 8 <= payload_end:
            version = data[payload_start]
            cursor = payload_start + 4

            if cursor + 2 > payload_end:
                continue
            offset_size = data[cursor] >> 4
            length_size = data[cursor] & 0x0F
            cursor += 1

            base_offset_size = data[cursor] >> 4
            index_size = (data[cursor] & 0x0F) if version in (1, 2) else 0
            cursor += 1

            if version < 2:
                if cursor + 2 > payload_end:
                    continue
                item_count = struct.unpack('>H', data[cursor:cursor + 2])[0]
                cursor += 2
            else:
                if cursor + 4 > payload_end:
                    continue
                item_count = struct.unpack('>L', data[cursor:cursor + 4])[0]
                cursor += 4

            for _ in range(item_count):
                if version < 2:
                    if cursor + 2 > payload_end:
                        break
                    item_id = struct.unpack('>H', data[cursor:cursor + 2])[0]
                    cursor += 2
                else:
                    if cursor + 4 > payload_end:
                        break
                    item_id = struct.unpack('>L', data[cursor:cursor + 4])[0]
                    cursor += 4

                if version in (1, 2):
                    if cursor + 2 > payload_end:
                        break
                    cursor += 2

                if cursor + 2 > payload_end:
                    break
                cursor += 2

                if cursor + base_offset_size > payload_end:
                    break
                base_offset = int.from_bytes(data[cursor:cursor + base_offset_size], 'big') if base_offset_size else 0
                cursor += base_offset_size

                if cursor + 2 > payload_end:
                    break
                extent_count = struct.unpack('>H', data[cursor:cursor + 2])[0]
                cursor += 2

                extents = []
                for _ in range(extent_count):
                    if version in (1, 2) and index_size:
                        if cursor + index_size > payload_end:
                            break
                        cursor += index_size
                    if cursor + offset_size + length_size > payload_end:
                        break
                    extent_offset = int.from_bytes(data[cursor:cursor + offset_size], 'big') if offset_size else 0
                    cursor += offset_size
                    extent_length = int.from_bytes(data[cursor:cursor + length_size], 'big') if length_size else 0
                    cursor += length_size
                    extents.append((base_offset + extent_offset, extent_length))

                if extents:
                    item_extents[item_id] = extents
        elif box_type == b'iprp':
            for p_offset, p_size, p_type, p_header in _iter_iso_boxes(data, payload_start, payload_end):
                p_payload_start = p_offset + p_header
                p_payload_end = p_offset + p_size
                if p_type == b'ipco':
                    properties = list(_iter_iso_boxes(data, p_payload_start, p_payload_end))
                elif p_type == b'ipma' and p_payload_start + 8 <= p_payload_end:
                    flags = int.from_bytes(data[p_payload_start + 1:p_payload_start + 4], 'big')
                    is_large_index = bool(flags & 1)
                    cursor = p_payload_start + 4
                    if cursor + 4 > p_payload_end:
                        continue
                    entry_count = struct.unpack('>L', data[cursor:cursor + 4])[0]
                    cursor += 4
                    for _ in range(entry_count):
                        if cursor + 3 > p_payload_end:
                            break
                        item_id = struct.unpack('>H', data[cursor:cursor + 2])[0]
                        cursor += 2
                        assoc_count = data[cursor]
                        cursor += 1
                        item_props = []
                        for _ in range(assoc_count):
                            if is_large_index:
                                if cursor + 2 > p_payload_end:
                                    break
                                value = struct.unpack('>H', data[cursor:cursor + 2])[0]
                                cursor += 2
                                item_props.append(value & 0x7FFF)
                            else:
                                if cursor + 1 > p_payload_end:
                                    break
                                value = data[cursor]
                                cursor += 1
                                item_props.append(value & 0x7F)
                        associations[item_id] = item_props

    if not properties:
        return -1, -1, -1, -1

    target_indexes = associations.get(primary_item_id, list(range(1, len(properties) + 1)))
    width = height = rotation = exif_rotation = -1
    for index in target_indexes:
        if not (1 <= index <= len(properties)):
            continue
        p_offset, p_size, p_type, p_header = properties[index - 1]
        p_payload_start = p_offset + p_header
        if p_type == b'ispe' and p_payload_start + 12 <= p_offset + p_size:
            width, height = struct.unpack('>LL', data[p_payload_start + 4:p_payload_start + 12])
        elif p_type == b'irot' and p_payload_start + 5 <= p_offset + p_size:
            rotation = _HEIF_IROT_TO_EXIF.get(data[p_payload_start + 4] & 0x03, -1)

    for item_id, item_type in item_types.items():
        if item_type != b'Exif':
            continue
        for extent_offset, extent_length in item_extents.get(item_id, []):
            if extent_length < 8:
                continue
            extent_end = extent_offset + extent_length
            if extent_offset < 0 or extent_end > len(data):
                continue
            exif_item = data[extent_offset:extent_end]
            tiff_offset = 4 + struct.unpack('>L', exif_item[:4])[0]
            if tiff_offset + 8 > len(exif_item):
                continue
            exif_rotation = _read_orientation_from_exif_payload(exif_item[tiff_offset:])
            if exif_rotation != -1:
                break
        if exif_rotation != -1:
            break

    return width, height, rotation, exif_rotation


def _read_tiff_rotation(fhandle):
    fhandle.seek(0)
    head = fhandle.read(16)
    if len(head) < 8:
        return -1

    if head.startswith(b"MM\x00*"):
        endian = '>'
        is_bigtiff = False
    elif head.startswith(b"II*\x00"):
        endian = '<'
        is_bigtiff = False
    elif head.startswith(b"II+\x00"):
        endian = '<'
        is_bigtiff = True
    else:
        return -1

    try:
        if is_bigtiff:
            if len(head) < 16:
                return -1
            bytesize = struct.unpack(endian + 'H', head[4:6])[0]
            if bytesize != 8:
                return -1
            ifd_offset = struct.unpack(endian + 'Q', head[8:16])[0]
            fhandle.seek(ifd_offset)
            entry_count = struct.unpack(endian + 'Q', fhandle.read(8))[0]
            for _ in range(entry_count):
                entry = fhandle.read(20)
                if len(entry) < 20:
                    return -1
                tag, datatype = struct.unpack(endian + 'HH', entry[:4])
                count = struct.unpack(endian + 'Q', entry[4:12])[0]
                value_field = entry[12:20]
                if tag == 274 and count == 1:
                    if datatype == 3:
                        return struct.unpack(endian + 'H', value_field[:2])[0]
                    if datatype == 4:
                        return struct.unpack(endian + 'L', value_field[:4])[0]
            return -1

        ifd_offset = struct.unpack(endian + 'L', head[4:8])[0]
        fhandle.seek(ifd_offset)
        entry_count = struct.unpack(endian + 'H', fhandle.read(2))[0]
        for _ in range(entry_count):
            entry = fhandle.read(12)
            if len(entry) < 12:
                return -1
            tag, datatype, count = struct.unpack(endian + 'HHL', entry[:8])
            value_field = entry[8:12]
            if tag == 274 and count == 1:
                if datatype == 3:
                    return struct.unpack(endian + 'H', value_field[:2])[0]
                if datatype == 4:
                    return struct.unpack(endian + 'L', value_field)[0]
    except struct.error:
        return -1

    return -1


def _get_rotation(fhandle):
    rotation = _read_jpeg_exif_rotation(fhandle)
    if rotation != -1:
        return rotation
    rotation = _read_heif_exif_rotation(fhandle)
    if rotation != -1:
        return rotation
    return _read_tiff_rotation(fhandle)


def _is_rotation_swapped(rotation):
    return rotation in {5, 6, 7, 8}


def _get_dpi(fhandle):
    xDPI = -1
    yDPI = -1

    fhandle.seek(0)
    head = fhandle.read(24)
    size = len(head)
    # handle GIFs
    # GIFs doesn't have density
    if size >= 10 and head[:6] in (b'GIF87a', b'GIF89a'):
        pass
    # see png edition spec bytes are below chunk length then and finally the
    elif size >= 24 and head.startswith(b'\211PNG\r\n\032\n'):
        chunkOffset = 8
        chunk = head[8:]
        while True:
            chunkType = chunk[4:8]
            if chunkType == b'pHYs':
                try:
                    xDensity, yDensity, unit = struct.unpack(">LLB", chunk[8:])
                except struct.error:
                    raise ValueError("Invalid PNG file")
                if unit:
                    xDPI = _convertToDPI(xDensity, _UNIT_1M)
                    yDPI = _convertToDPI(yDensity, _UNIT_1M)
                else:  # no unit
                    xDPI = xDensity
                    yDPI = yDensity
                break
            elif chunkType == b'IDAT':
                break
            else:
                try:
                    dataSize, = struct.unpack(">L", chunk[0:4])
                except struct.error:
                    raise ValueError("Invalid PNG file")
                chunkOffset += dataSize + 12
                fhandle.seek(chunkOffset)
                chunk = fhandle.read(17)
    # handle JPEGs
    elif size >= 2 and head.startswith(b'\377\330'):
        try:
            fhandle.seek(0)
            block_size = 2
            marker = 0
            while not 0xc0 <= marker <= 0xcf:
                fhandle.seek(block_size, 1)
                marker, block_size = _read_jpeg_segment_header(fhandle)
                if marker == 0xe0:  # APP0 marker
                    fhandle.seek(7, 1)
                    unit, xDensity, yDensity = struct.unpack(">BHH", fhandle.read(5))
                    if unit == 1 or unit == 0:
                        xDPI = xDensity
                        yDPI = yDensity
                    elif unit == 2:
                        xDPI = _convertToDPI(xDensity, _UNIT_CM)
                        yDPI = _convertToDPI(yDensity, _UNIT_CM)
                    break
        except (struct.error, ValueError):
            raise ValueError("Invalid JPEG file")
    # handle JPEG2000s
    elif size >= 12 and head.startswith(b'\x00\x00\x00\x0cjP  \r\n\x87\n'):
        fhandle.seek(32)
        # skip JP2 image header box
        headerSize = struct.unpack('>L', fhandle.read(4))[0] - 8
        fhandle.seek(4, 1)
        foundResBox = False
        try:
            while headerSize > 0:
                boxHeader = fhandle.read(8)
                boxType = boxHeader[4:]
                if boxType == b'res ':  # find resolution super box
                    foundResBox = True
                    headerSize -= 8
                    break
                boxSize, = struct.unpack('>L', boxHeader[:4])
                fhandle.seek(boxSize - 8, 1)
                headerSize -= boxSize
            if foundResBox:
                while headerSize > 0:
                    boxHeader = fhandle.read(8)
                    boxType = boxHeader[4:]
                    if boxType == b'resd':  # Display resolution box
                        yDensity, xDensity, yUnit, xUnit = struct.unpack(">HHBB", fhandle.read(10))
                        xDPI = _convertToDPI(xDensity, xUnit)
                        yDPI = _convertToDPI(yDensity, yUnit)
                        break
                    boxSize, = struct.unpack('>L', boxHeader[:4])
                    fhandle.seek(boxSize - 8, 1)
                    headerSize -= boxSize
        except struct.error:
            raise ValueError("Invalid JPEG2000 file")

    return xDPI, yDPI


def _get_colors(fhandle):
    colors = -1
    fhandle.seek(0)
    head = fhandle.read(32)
    if len(head) >= 11 and head[:6] in (b'GIF87a', b'GIF89a'):
        packed = head[10]
        if packed & 0x80:
            colors = 2 ** ((packed & 0x07) + 1)
    elif len(head) >= 26 and head.startswith(b'\211PNG\r\n\032\n') and head[12:16] == b'IHDR':
        bit_depth = head[24]
        color_type = head[25]
        channels = {
            0: 1,
            2: 3,
            3: 1,
            4: 1,
            6: 3,
        }.get(color_type)
        if channels:
            colors = 2 ** (bit_depth * channels)
    return colors


def _get_channels(fhandle):
    channels = -1

    fhandle.seek(0)
    head = fhandle.read(32)
    size = len(head)

    if size >= 26 and head.startswith(b'\211PNG\r\n\032\n') and head[12:16] == b'IHDR':
        color_type = head[25]
        channels = {
            0: 1,
            2: 3,
            3: 1,
            4: 2,
            6: 4,
        }.get(color_type, -1)
    elif size >= 2 and head.startswith(b'\377\330'):
        try:
            fhandle.seek(0)
            _seek_to_jpeg_sof(fhandle)
            fhandle.seek(5, 1)
            channels = struct.unpack('>B', fhandle.read(1))[0]
        except (struct.error, ValueError):
            raise ValueError("Invalid JPEG file")
    elif size >= 11 and head[:6] in (b'GIF87a', b'GIF89a'):
        channels = 3
    elif size >= 26 and head.startswith(b'BM'):
        bit_depth = struct.unpack('<H', head[28:30])[0]
        if bit_depth <= 8:
            channels = 1
        elif bit_depth == 24:
            channels = 3
        elif bit_depth == 32:
            channels = 4

    return channels

  
def get_info(filepath: FileInput, *, size: bool = True, dpi: bool = True, colors: bool = True,
             exif_rotation: bool = True, channels: bool = True) -> ImageInfo:
    fhandle, should_close = _open_file(filepath)
    try:
        width = height = rotation = xdpi = ydpi = color_count = channel_count = -1
        if size:
            width, height = _get_size(fhandle)
            rotation = _get_rotation(fhandle)
            if exif_rotation and _is_rotation_swapped(rotation):
                width, height = height, width
        if dpi:
            xdpi, ydpi = _get_dpi(fhandle)
        if colors:
            color_count = _get_colors(fhandle)
        if channels:
            channel_count = _get_channels(fhandle)
        return ImageInfo(width=width, height=height, rotation=rotation, xdpi=xdpi, ydpi=ydpi, colors=color_count, channels=channel_count)
    finally:
        if should_close:
            fhandle.close()


def get(filepath: FileInput, *, exif_rotation: bool = True) -> Tuple[int, int]:
    """
    Return (width, height) for a given img file content.
    Set exif_rotation=False to return stored dimensions as-is.
    :type filepath: Union[bytes, str, pathlib.Path]
    :rtype Tuple[int, int]
    """
    try:
        info = get_info(filepath, size=True, dpi=False, colors=False, exif_rotation=exif_rotation, channels=False)
    except Exception:
        return -1, -1
    return info.width, info.height


def getDPI(filepath: FileInput) -> Tuple[int, int]:
    """
    Return (x DPI, y DPI) for a given img file content
    no requirements
    :type filepath: Union[bytes, str, pathlib.Path]
    :rtype Tuple[int, int]
    """
    try:
        info = get_info(filepath, size=False, dpi=True, colors=False)
    except Exception:
        return -1, -1
    return info.xdpi, info.ydpi
