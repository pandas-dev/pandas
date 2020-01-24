from typing import Dict, List, Tuple

magic: bytes = (
    b"\x00\x00\x00\x00\x00\x00\x00\x00"
    b"\x00\x00\x00\x00\xc2\xea\x81\x60"
    b"\xb3\x14\x11\xcf\xbd\x92\x08\x00"
    b"\x09\xc7\x31\x8c\x18\x1f\x10\x11"
)

align_1_checker_value: bytes = b"3"
align_1_offset: int = 32
align_1_length: int = 1
align_1_value: int = 4
u64_byte_checker_value: bytes = b"3"
align_2_offset: int = 35
align_2_length: int = 1
align_2_value: int = 4
endianness_offset: int = 37
endianness_length: int = 1
platform_offset: int = 39
platform_length: int = 1
encoding_offset: int = 70
encoding_length: int = 1
dataset_offset: int = 92
dataset_length: int = 64
file_type_offset: int = 156
file_type_length: int = 8
date_created_offset: int = 164
date_created_length: int = 8
date_modified_offset: int = 172
date_modified_length: int = 8
header_size_offset: int = 196
header_size_length: int = 4
page_size_offset: int = 200
page_size_length: int = 4
page_count_offset: int = 204
page_count_length: int = 4
sas_release_offset: int = 216
sas_release_length: int = 8
sas_server_type_offset: int = 224
sas_server_type_length: int = 16
os_version_number_offset: int = 240
os_version_number_length: int = 16
os_maker_offset: int = 256
os_maker_length: int = 16
os_name_offset: int = 272
os_name_length: int = 16
page_bit_offset_x86: int = 16
page_bit_offset_x64: int = 32
subheader_pointer_length_x86: int = 12
subheader_pointer_length_x64: int = 24
page_type_offset: int = 0
page_type_length: int = 2
block_count_offset: int = 2
block_count_length: int = 2
subheader_count_offset: int = 4
subheader_count_length: int = 2
page_meta_type: int = 0
page_data_type: int = 256
page_amd_type: int = 1024
page_metc_type: int = 16384
page_comp_type: int = -28672
page_mix_types: List[int] = [512, 640]
subheader_pointers_offset: int = 8
truncated_subheader_id: int = 1
compressed_subheader_id: int = 4
compressed_subheader_type: int = 1
text_block_size_length: int = 2
row_length_offset_multiplier: int = 5
row_count_offset_multiplier: int = 6
col_count_p1_multiplier: int = 9
col_count_p2_multiplier: int = 10
row_count_on_mix_page_offset_multiplier: int = 15
column_name_pointer_length: int = 8
column_name_text_subheader_offset: int = 0
column_name_text_subheader_length: int = 2
column_name_offset_offset: int = 2
column_name_offset_length: int = 2
column_name_length_offset: int = 4
column_name_length_length: int = 2
column_data_offset_offset: int = 8
column_data_length_offset: int = 8
column_data_length_length: int = 4
column_type_offset: int = 14
column_type_length: int = 1
column_format_text_subheader_index_offset: int = 22
column_format_text_subheader_index_length: int = 2
column_format_offset_offset: int = 24
column_format_offset_length: int = 2
column_format_length_offset: int = 26
column_format_length_length: int = 2
column_label_text_subheader_index_offset: int = 28
column_label_text_subheader_index_length: int = 2
column_label_offset_offset: int = 30
column_label_offset_length: int = 2
column_label_length_offset: int = 32
column_label_length_length: int = 2
rle_compression: bytes = b"SASYZCRL"
rdc_compression: bytes = b"SASYZCR2"

compression_literals: List[bytes] = [rle_compression, rdc_compression]

# Incomplete list of encodings, using SAS nomenclature:
# http://support.sas.com/documentation/cdl/en/nlsref/61893/HTML/default/viewer.htm#a002607278.htm
encoding_names: Dict[int, str] = {
    29: "latin1",
    20: "utf-8",
    33: "cyrillic",
    60: "wlatin2",
    61: "wcyrillic",
    62: "wlatin1",
    90: "ebcdic870",
}


class SASIndex:
    row_size_index: int = 0
    column_size_index: int = 1
    subheader_counts_index: int = 2
    column_text_index: int = 3
    column_name_index: int = 4
    column_attributes_index: int = 5
    format_and_label_index: int = 6
    column_list_index: int = 7
    data_subheader_index: int = 8


subheader_signature_to_index: Dict[bytes, int] = {
    b"\xF7\xF7\xF7\xF7": SASIndex.row_size_index,
    b"\x00\x00\x00\x00\xF7\xF7\xF7\xF7": SASIndex.row_size_index,
    b"\xF7\xF7\xF7\xF7\x00\x00\x00\x00": SASIndex.row_size_index,
    b"\xF7\xF7\xF7\xF7\xFF\xFF\xFB\xFE": SASIndex.row_size_index,
    b"\xF6\xF6\xF6\xF6": SASIndex.column_size_index,
    b"\x00\x00\x00\x00\xF6\xF6\xF6\xF6": SASIndex.column_size_index,
    b"\xF6\xF6\xF6\xF6\x00\x00\x00\x00": SASIndex.column_size_index,
    b"\xF6\xF6\xF6\xF6\xFF\xFF\xFB\xFE": SASIndex.column_size_index,
    b"\x00\xFC\xFF\xFF": SASIndex.subheader_counts_index,
    b"\xFF\xFF\xFC\x00": SASIndex.subheader_counts_index,
    b"\x00\xFC\xFF\xFF\xFF\xFF\xFF\xFF": SASIndex.subheader_counts_index,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFC\x00": SASIndex.subheader_counts_index,
    b"\xFD\xFF\xFF\xFF": SASIndex.column_text_index,
    b"\xFF\xFF\xFF\xFD": SASIndex.column_text_index,
    b"\xFD\xFF\xFF\xFF\xFF\xFF\xFF\xFF": SASIndex.column_text_index,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFD": SASIndex.column_text_index,
    b"\xFF\xFF\xFF\xFF": SASIndex.column_name_index,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF": SASIndex.column_name_index,
    b"\xFC\xFF\xFF\xFF": SASIndex.column_attributes_index,
    b"\xFF\xFF\xFF\xFC": SASIndex.column_attributes_index,
    b"\xFC\xFF\xFF\xFF\xFF\xFF\xFF\xFF": SASIndex.column_attributes_index,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFC": SASIndex.column_attributes_index,
    b"\xFE\xFB\xFF\xFF": SASIndex.format_and_label_index,
    b"\xFF\xFF\xFB\xFE": SASIndex.format_and_label_index,
    b"\xFE\xFB\xFF\xFF\xFF\xFF\xFF\xFF": SASIndex.format_and_label_index,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFB\xFE": SASIndex.format_and_label_index,
    b"\xFE\xFF\xFF\xFF": SASIndex.column_list_index,
    b"\xFF\xFF\xFF\xFE": SASIndex.column_list_index,
    b"\xFE\xFF\xFF\xFF\xFF\xFF\xFF\xFF": SASIndex.column_list_index,
    b"\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFE": SASIndex.column_list_index,
}


# List of frequently used SAS date and datetime formats
# http://support.sas.com/documentation/cdl/en/etsug/60372/HTML/default/viewer.htm#etsug_intervals_sect009.htm
# https://github.com/epam/parso/blob/master/src/main/java/com/epam/parso/impl/SasFileConstants.java
sas_date_formats: Tuple[str, ...] = (
    "DATE",
    "DAY",
    "DDMMYY",
    "DOWNAME",
    "JULDAY",
    "JULIAN",
    "MMDDYY",
    "MMYY",
    "MMYYC",
    "MMYYD",
    "MMYYP",
    "MMYYS",
    "MMYYN",
    "MONNAME",
    "MONTH",
    "MONYY",
    "QTR",
    "QTRR",
    "NENGO",
    "WEEKDATE",
    "WEEKDATX",
    "WEEKDAY",
    "WEEKV",
    "WORDDATE",
    "WORDDATX",
    "YEAR",
    "YYMM",
    "YYMMC",
    "YYMMD",
    "YYMMP",
    "YYMMS",
    "YYMMN",
    "YYMON",
    "YYMMDD",
    "YYQ",
    "YYQC",
    "YYQD",
    "YYQP",
    "YYQS",
    "YYQN",
    "YYQR",
    "YYQRC",
    "YYQRD",
    "YYQRP",
    "YYQRS",
    "YYQRN",
    "YYMMDDP",
    "YYMMDDC",
    "E8601DA",
    "YYMMDDN",
    "MMDDYYC",
    "MMDDYYS",
    "MMDDYYD",
    "YYMMDDS",
    "B8601DA",
    "DDMMYYN",
    "YYMMDDD",
    "DDMMYYB",
    "DDMMYYP",
    "MMDDYYP",
    "YYMMDDB",
    "MMDDYYN",
    "DDMMYYC",
    "DDMMYYD",
    "DDMMYYS",
    "MINGUO",
)

sas_datetime_formats: Tuple[str, ...] = (
    "DATETIME",
    "DTWKDATX",
    "B8601DN",
    "B8601DT",
    "B8601DX",
    "B8601DZ",
    "B8601LX",
    "E8601DN",
    "E8601DT",
    "E8601DX",
    "E8601DZ",
    "E8601LX",
    "DATEAMPM",
    "DTDATE",
    "DTMONYY",
    "DTMONYY",
    "DTWKDATX",
    "DTYEAR",
    "TOD",
    "MDYAMPM",
)
