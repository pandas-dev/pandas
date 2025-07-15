class Type(object):
    """
    Types supported by Parquet.  These types are intended to be used in combination
    with the encodings to control the on disk storage format.
    For example INT16 is not included as a type since a good encoding of INT32
    would handle this.

    """
    BOOLEAN = 0
    INT32 = 1
    INT64 = 2
    INT96 = 3
    FLOAT = 4
    DOUBLE = 5
    BYTE_ARRAY = 6
    FIXED_LEN_BYTE_ARRAY = 7

    _VALUES_TO_NAMES = {
        0: "BOOLEAN",
        1: "INT32",
        2: "INT64",
        3: "INT96",
        4: "FLOAT",
        5: "DOUBLE",
        6: "BYTE_ARRAY",
        7: "FIXED_LEN_BYTE_ARRAY",
    }

    _NAMES_TO_VALUES = {
        "BOOLEAN": 0,
        "INT32": 1,
        "INT64": 2,
        "INT96": 3,
        "FLOAT": 4,
        "DOUBLE": 5,
        "BYTE_ARRAY": 6,
        "FIXED_LEN_BYTE_ARRAY": 7,
    }


class ConvertedType(object):
    """
    DEPRECATED: Common types used by frameworks(e.g. hive, pig) using parquet.
    ConvertedType is superseded by LogicalType.  This enum should not be extended.

    See LogicalTypes.md for conversion between ConvertedType and LogicalType.

    """
    UTF8 = 0
    MAP = 1
    MAP_KEY_VALUE = 2
    LIST = 3
    ENUM = 4
    DECIMAL = 5
    DATE = 6
    TIME_MILLIS = 7
    TIME_MICROS = 8
    TIMESTAMP_MILLIS = 9
    TIMESTAMP_MICROS = 10
    UINT_8 = 11
    UINT_16 = 12
    UINT_32 = 13
    UINT_64 = 14
    INT_8 = 15
    INT_16 = 16
    INT_32 = 17
    INT_64 = 18
    JSON = 19
    BSON = 20
    INTERVAL = 21

    _VALUES_TO_NAMES = {
        0: "UTF8",
        1: "MAP",
        2: "MAP_KEY_VALUE",
        3: "LIST",
        4: "ENUM",
        5: "DECIMAL",
        6: "DATE",
        7: "TIME_MILLIS",
        8: "TIME_MICROS",
        9: "TIMESTAMP_MILLIS",
        10: "TIMESTAMP_MICROS",
        11: "UINT_8",
        12: "UINT_16",
        13: "UINT_32",
        14: "UINT_64",
        15: "INT_8",
        16: "INT_16",
        17: "INT_32",
        18: "INT_64",
        19: "JSON",
        20: "BSON",
        21: "INTERVAL",
    }

    _NAMES_TO_VALUES = {
        "UTF8": 0,
        "MAP": 1,
        "MAP_KEY_VALUE": 2,
        "LIST": 3,
        "ENUM": 4,
        "DECIMAL": 5,
        "DATE": 6,
        "TIME_MILLIS": 7,
        "TIME_MICROS": 8,
        "TIMESTAMP_MILLIS": 9,
        "TIMESTAMP_MICROS": 10,
        "UINT_8": 11,
        "UINT_16": 12,
        "UINT_32": 13,
        "UINT_64": 14,
        "INT_8": 15,
        "INT_16": 16,
        "INT_32": 17,
        "INT_64": 18,
        "JSON": 19,
        "BSON": 20,
        "INTERVAL": 21,
    }


class FieldRepetitionType(object):
    """
    Representation of Schemas

    """
    REQUIRED = 0
    OPTIONAL = 1
    REPEATED = 2

    _VALUES_TO_NAMES = {
        0: "REQUIRED",
        1: "OPTIONAL",
        2: "REPEATED",
    }

    _NAMES_TO_VALUES = {
        "REQUIRED": 0,
        "OPTIONAL": 1,
        "REPEATED": 2,
    }


class Encoding(object):
    """
    Encodings supported by Parquet.  Not all encodings are valid for all types.  These
    enums are also used to specify the encoding of definition and repetition levels.
    See the accompanying doc for the details of the more complicated encodings.

    """
    PLAIN = 0
    PLAIN_DICTIONARY = 2
    RLE = 3
    BIT_PACKED = 4
    DELTA_BINARY_PACKED = 5
    DELTA_LENGTH_BYTE_ARRAY = 6
    DELTA_BYTE_ARRAY = 7
    RLE_DICTIONARY = 8
    BYTE_STREAM_SPLIT = 9

    _VALUES_TO_NAMES = {
        0: "PLAIN",
        2: "PLAIN_DICTIONARY",
        3: "RLE",
        4: "BIT_PACKED",
        5: "DELTA_BINARY_PACKED",
        6: "DELTA_LENGTH_BYTE_ARRAY",
        7: "DELTA_BYTE_ARRAY",
        8: "RLE_DICTIONARY",
        9: "BYTE_STREAM_SPLIT",
    }

    _NAMES_TO_VALUES = {
        "PLAIN": 0,
        "PLAIN_DICTIONARY": 2,
        "RLE": 3,
        "BIT_PACKED": 4,
        "DELTA_BINARY_PACKED": 5,
        "DELTA_LENGTH_BYTE_ARRAY": 6,
        "DELTA_BYTE_ARRAY": 7,
        "RLE_DICTIONARY": 8,
        "BYTE_STREAM_SPLIT": 9,
    }


class CompressionCodec(object):
    """
    Supported compression algorithms.

    Codecs added in format version X.Y can be read by readers based on X.Y and later.
    Codec support may vary between readers based on the format version and
    libraries available at runtime.

    See Compression.md for a detailed specification of these algorithms.

    """
    UNCOMPRESSED = 0
    SNAPPY = 1
    GZIP = 2
    LZO = 3
    BROTLI = 4
    LZ4 = 5
    ZSTD = 6
    LZ4_RAW = 7

    _VALUES_TO_NAMES = {
        0: "UNCOMPRESSED",
        1: "SNAPPY",
        2: "GZIP",
        3: "LZO",
        4: "BROTLI",
        5: "LZ4",
        6: "ZSTD",
        7: "LZ4_RAW",
    }

    _NAMES_TO_VALUES = {
        "UNCOMPRESSED": 0,
        "SNAPPY": 1,
        "GZIP": 2,
        "LZO": 3,
        "BROTLI": 4,
        "LZ4": 5,
        "ZSTD": 6,
        "LZ4_RAW": 7,
    }


class PageType(object):
    DATA_PAGE = 0
    INDEX_PAGE = 1
    DICTIONARY_PAGE = 2
    DATA_PAGE_V2 = 3

    _VALUES_TO_NAMES = {
        0: "DATA_PAGE",
        1: "INDEX_PAGE",
        2: "DICTIONARY_PAGE",
        3: "DATA_PAGE_V2",
    }

    _NAMES_TO_VALUES = {
        "DATA_PAGE": 0,
        "INDEX_PAGE": 1,
        "DICTIONARY_PAGE": 2,
        "DATA_PAGE_V2": 3,
    }


class BoundaryOrder(object):
    """
    Enum to annotate whether lists of min/max elements inside ColumnIndex
    are ordered and if so, in which direction.

    """
    UNORDERED = 0
    ASCENDING = 1
    DESCENDING = 2

    _VALUES_TO_NAMES = {
        0: "UNORDERED",
        1: "ASCENDING",
        2: "DESCENDING",
    }

    _NAMES_TO_VALUES = {
        "UNORDERED": 0,
        "ASCENDING": 1,
        "DESCENDING": 2,
    }

