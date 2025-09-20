# Python implementation of low level MySQL client-server protocol
# http://dev.mysql.com/doc/internals/en/client-server-protocol.html

from .charset import MBLENGTH
from .constants import FIELD_TYPE, SERVER_STATUS
from . import err

import struct
import sys


DEBUG = False

NULL_COLUMN = 251
UNSIGNED_CHAR_COLUMN = 251
UNSIGNED_SHORT_COLUMN = 252
UNSIGNED_INT24_COLUMN = 253
UNSIGNED_INT64_COLUMN = 254


def dump_packet(data):  # pragma: no cover
    def printable(data):
        if 32 <= data < 127:
            return chr(data)
        return "."

    try:
        print("packet length:", len(data))
        for i in range(1, 7):
            f = sys._getframe(i)
            print("call[%d]: %s (line %d)" % (i, f.f_code.co_name, f.f_lineno))
        print("-" * 66)
    except ValueError:
        pass
    dump_data = [data[i : i + 16] for i in range(0, min(len(data), 256), 16)]
    for d in dump_data:
        print(
            " ".join(f"{x:02X}" for x in d)
            + "   " * (16 - len(d))
            + " " * 2
            + "".join(printable(x) for x in d)
        )
    print("-" * 66)
    print()


class MysqlPacket:
    """Representation of a MySQL response packet.

    Provides an interface for reading/parsing the packet results.
    """

    __slots__ = ("_position", "_data")

    def __init__(self, data, encoding):
        self._position = 0
        self._data = data

    def get_all_data(self):
        return self._data

    def read(self, size):
        """Read the first 'size' bytes in packet and advance cursor past them."""
        result = self._data[self._position : (self._position + size)]
        if len(result) != size:
            error = (
                "Result length not requested length:\n"
                f"Expected={size}.  Actual={len(result)}.  Position: {self._position}.  Data Length: {len(self._data)}"
            )
            if DEBUG:
                print(error)
                self.dump()
            raise AssertionError(error)
        self._position += size
        return result

    def read_all(self):
        """Read all remaining data in the packet.

        (Subsequent read() will return errors.)
        """
        result = self._data[self._position :]
        self._position = None  # ensure no subsequent read()
        return result

    def advance(self, length):
        """Advance the cursor in data buffer 'length' bytes."""
        new_position = self._position + length
        if new_position < 0 or new_position > len(self._data):
            raise Exception(
                f"Invalid advance amount ({length}) for cursor.  Position={new_position}"
            )
        self._position = new_position

    def rewind(self, position=0):
        """Set the position of the data buffer cursor to 'position'."""
        if position < 0 or position > len(self._data):
            raise Exception("Invalid position to rewind cursor to: %s." % position)
        self._position = position

    def get_bytes(self, position, length=1):
        """Get 'length' bytes starting at 'position'.

        Position is start of payload (first four packet header bytes are not
        included) starting at index '0'.

        No error checking is done.  If requesting outside end of buffer
        an empty string (or string shorter than 'length') may be returned!
        """
        return self._data[position : (position + length)]

    def read_uint8(self):
        result = self._data[self._position]
        self._position += 1
        return result

    def read_uint16(self):
        result = struct.unpack_from("<H", self._data, self._position)[0]
        self._position += 2
        return result

    def read_uint24(self):
        low, high = struct.unpack_from("<HB", self._data, self._position)
        self._position += 3
        return low + (high << 16)

    def read_uint32(self):
        result = struct.unpack_from("<I", self._data, self._position)[0]
        self._position += 4
        return result

    def read_uint64(self):
        result = struct.unpack_from("<Q", self._data, self._position)[0]
        self._position += 8
        return result

    def read_string(self):
        end_pos = self._data.find(b"\0", self._position)
        if end_pos < 0:
            return None
        result = self._data[self._position : end_pos]
        self._position = end_pos + 1
        return result

    def read_length_encoded_integer(self):
        """Read a 'Length Coded Binary' number from the data buffer.

        Length coded numbers can be anywhere from 1 to 9 bytes depending
        on the value of the first byte.
        """
        c = self.read_uint8()
        if c == NULL_COLUMN:
            return None
        if c < UNSIGNED_CHAR_COLUMN:
            return c
        elif c == UNSIGNED_SHORT_COLUMN:
            return self.read_uint16()
        elif c == UNSIGNED_INT24_COLUMN:
            return self.read_uint24()
        elif c == UNSIGNED_INT64_COLUMN:
            return self.read_uint64()

    def read_length_coded_string(self):
        """Read a 'Length Coded String' from the data buffer.

        A 'Length Coded String' consists first of a length coded
        (unsigned, positive) integer represented in 1-9 bytes followed by
        that many bytes of binary data.  (For example "cat" would be "3cat".)
        """
        length = self.read_length_encoded_integer()
        if length is None:
            return None
        return self.read(length)

    def read_struct(self, fmt):
        s = struct.Struct(fmt)
        result = s.unpack_from(self._data, self._position)
        self._position += s.size
        return result

    def is_ok_packet(self):
        # https://dev.mysql.com/doc/internals/en/packet-OK_Packet.html
        return self._data[0] == 0 and len(self._data) >= 7

    def is_eof_packet(self):
        # http://dev.mysql.com/doc/internals/en/generic-response-packets.html#packet-EOF_Packet
        # Caution: \xFE may be LengthEncodedInteger.
        # If \xFE is LengthEncodedInteger header, 8bytes followed.
        return self._data[0] == 0xFE and len(self._data) < 9

    def is_auth_switch_request(self):
        # http://dev.mysql.com/doc/internals/en/connection-phase-packets.html#packet-Protocol::AuthSwitchRequest
        return self._data[0] == 0xFE

    def is_extra_auth_data(self):
        # https://dev.mysql.com/doc/internals/en/successful-authentication.html
        return self._data[0] == 1

    def is_resultset_packet(self):
        field_count = self._data[0]
        return 1 <= field_count <= 250

    def is_load_local_packet(self):
        return self._data[0] == 0xFB

    def is_error_packet(self):
        return self._data[0] == 0xFF

    def check_error(self):
        if self.is_error_packet():
            self.raise_for_error()

    def raise_for_error(self):
        self.rewind()
        self.advance(1)  # field_count == error (we already know that)
        errno = self.read_uint16()
        if DEBUG:
            print("errno =", errno)
        err.raise_mysql_exception(self._data)

    def dump(self):
        dump_packet(self._data)


class FieldDescriptorPacket(MysqlPacket):
    """A MysqlPacket that represents a specific column's metadata in the result.

    Parsing is automatically done and the results are exported via public
    attributes on the class such as: db, table_name, name, length, type_code.
    """

    def __init__(self, data, encoding):
        MysqlPacket.__init__(self, data, encoding)
        self._parse_field_descriptor(encoding)

    def _parse_field_descriptor(self, encoding):
        """Parse the 'Field Descriptor' (Metadata) packet.

        This is compatible with MySQL 4.1+ (not compatible with MySQL 4.0).
        """
        self.catalog = self.read_length_coded_string()
        self.db = self.read_length_coded_string()
        self.table_name = self.read_length_coded_string().decode(encoding)
        self.org_table = self.read_length_coded_string().decode(encoding)
        self.name = self.read_length_coded_string().decode(encoding)
        self.org_name = self.read_length_coded_string().decode(encoding)
        (
            self.charsetnr,
            self.length,
            self.type_code,
            self.flags,
            self.scale,
        ) = self.read_struct("<xHIBHBxx")
        # 'default' is a length coded binary and is still in the buffer?
        # not used for normal result sets...

    def description(self):
        """Provides a 7-item tuple compatible with the Python PEP249 DB Spec."""
        return (
            self.name,
            self.type_code,
            None,  # TODO: display_length; should this be self.length?
            self.get_column_length(),  # 'internal_size'
            self.get_column_length(),  # 'precision'  # TODO: why!?!?
            self.scale,
            self.flags % 2 == 0,
        )

    def get_column_length(self):
        if self.type_code == FIELD_TYPE.VAR_STRING:
            mblen = MBLENGTH.get(self.charsetnr, 1)
            return self.length // mblen
        return self.length

    def __str__(self):
        return "{} {!r}.{!r}.{!r}, type={}, flags={:x}".format(
            self.__class__,
            self.db,
            self.table_name,
            self.name,
            self.type_code,
            self.flags,
        )


class OKPacketWrapper:
    """
    OK Packet Wrapper. It uses an existing packet object, and wraps
    around it, exposing useful variables while still providing access
    to the original packet objects variables and methods.
    """

    def __init__(self, from_packet):
        if not from_packet.is_ok_packet():
            raise ValueError(
                "Cannot create "
                + str(self.__class__.__name__)
                + " object from invalid packet type"
            )

        self.packet = from_packet
        self.packet.advance(1)

        self.affected_rows = self.packet.read_length_encoded_integer()
        self.insert_id = self.packet.read_length_encoded_integer()
        self.server_status, self.warning_count = self.read_struct("<HH")
        self.message = self.packet.read_all()
        self.has_next = self.server_status & SERVER_STATUS.SERVER_MORE_RESULTS_EXISTS

    def __getattr__(self, key):
        return getattr(self.packet, key)


class EOFPacketWrapper:
    """
    EOF Packet Wrapper. It uses an existing packet object, and wraps
    around it, exposing useful variables while still providing access
    to the original packet objects variables and methods.
    """

    def __init__(self, from_packet):
        if not from_packet.is_eof_packet():
            raise ValueError(
                f"Cannot create '{self.__class__}' object from invalid packet type"
            )

        self.packet = from_packet
        self.warning_count, self.server_status = self.packet.read_struct("<xhh")
        if DEBUG:
            print("server_status=", self.server_status)
        self.has_next = self.server_status & SERVER_STATUS.SERVER_MORE_RESULTS_EXISTS

    def __getattr__(self, key):
        return getattr(self.packet, key)


class LoadLocalPacketWrapper:
    """
    Load Local Packet Wrapper. It uses an existing packet object, and wraps
    around it, exposing useful variables while still providing access
    to the original packet objects variables and methods.
    """

    def __init__(self, from_packet):
        if not from_packet.is_load_local_packet():
            raise ValueError(
                f"Cannot create '{self.__class__}' object from invalid packet type"
            )

        self.packet = from_packet
        self.filename = self.packet.get_all_data()[1:]
        if DEBUG:
            print("filename=", self.filename)

    def __getattr__(self, key):
        return getattr(self.packet, key)
