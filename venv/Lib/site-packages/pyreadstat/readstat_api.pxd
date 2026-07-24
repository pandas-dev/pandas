# cython: c_string_type=unicode, c_string_encoding=utf8, language_level=2

# #############################################################################
# Copyright 2018 Hoffmann-La Roche
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# #############################################################################

"""
These are C definitions coming from the readstat API (readstat.h) plus some
dependecies of that to the C standard library
"""

cdef extern from "<sys/types.h>":
    ctypedef long off_t

cdef extern from "<stdint.h>" nogil:

    # 7.18.1 Integer types
    # 7.18.1.1 Exact-width integer types
    ctypedef   signed char  int8_t
    ctypedef   signed short int16_t
    ctypedef   signed int   int32_t
    ctypedef   signed long  int64_t
    ctypedef unsigned char  uint8_t

#cdef extern from "conditional_includes.h":
#    int close(int fd)
#    ssize_t write(int fd, const void *buf, size_t nbyte)



cdef extern from "readstat.h":

    ctypedef enum:
        READSTAT_HANDLER_OK,
        READSTAT_HANDLER_ABORT,
        READSTAT_HANDLER_SKIP_VARIABLE

    ctypedef struct readstat_variable_t:
        int index

    ctypedef struct readstat_label_set_t:
        pass
        
    ctypedef struct readstat_metadata_t:
        pass

    ctypedef struct mr_set_t:
        char type
        char *name
        char *label
        int is_dichotomy
        int counted_value
        char **subvariables
        int num_subvars

    ctypedef enum readstat_io_flags_t:
        READSTAT_SEEK_SET,
        READSTAT_SEEK_CUR,
        READSTAT_SEEK_END
    
    ctypedef enum readstat_error_t:
        READSTAT_OK,
        READSTAT_ERROR_OPEN = 1,
        READSTAT_ERROR_READ,
        READSTAT_ERROR_MALLOC,
        READSTAT_ERROR_USER_ABORT,
        READSTAT_ERROR_PARSE,
        READSTAT_ERROR_UNSUPPORTED_COMPRESSION,
        READSTAT_ERROR_UNSUPPORTED_CHARSET,
        READSTAT_ERROR_COLUMN_COUNT_MISMATCH,
        READSTAT_ERROR_ROW_COUNT_MISMATCH,
        READSTAT_ERROR_ROW_WIDTH_MISMATCH,
        READSTAT_ERROR_BAD_FORMAT_STRING,
        READSTAT_ERROR_VALUE_TYPE_MISMATCH,
        READSTAT_ERROR_WRITE,
        READSTAT_ERROR_WRITER_NOT_INITIALIZED,
        READSTAT_ERROR_SEEK,
        READSTAT_ERROR_CONVERT,
        READSTAT_ERROR_CONVERT_BAD_STRING,
        READSTAT_ERROR_CONVERT_SHORT_STRING,
        READSTAT_ERROR_CONVERT_LONG_STRING,
        READSTAT_ERROR_NUMERIC_VALUE_IS_OUT_OF_RANGE,
        READSTAT_ERROR_TAGGED_VALUE_IS_OUT_OF_RANGE,
        READSTAT_ERROR_STRING_VALUE_IS_TOO_LONG,
        READSTAT_ERROR_TAGGED_VALUES_NOT_SUPPORTED,
        READSTAT_ERROR_UNSUPPORTED_FILE_FORMAT_VERSION,
        READSTAT_ERROR_NAME_BEGINS_WITH_ILLEGAL_CHARACTER,
        READSTAT_ERROR_NAME_CONTAINS_ILLEGAL_CHARACTER,
        READSTAT_ERROR_NAME_IS_RESERVED_WORD,
        READSTAT_ERROR_NAME_IS_TOO_LONG,
        READSTAT_ERROR_BAD_TIMESTAMP,
        READSTAT_ERROR_BAD_FREQUENCY_WEIGHT,
        READSTAT_ERROR_TOO_MANY_MISSING_VALUE_DEFINITIONS,
        READSTAT_ERROR_NOTE_IS_TOO_LONG,
        READSTAT_ERROR_STRING_REFS_NOT_SUPPORTED,
        READSTAT_ERROR_STRING_REF_IS_REQUIRED,
        READSTAT_ERROR_ROW_IS_TOO_WIDE_FOR_PAGE
        
    ctypedef enum readstat_type_t:
        READSTAT_TYPE_STRING,
        READSTAT_TYPE_INT8,
        READSTAT_TYPE_INT16,
        READSTAT_TYPE_INT32,
        READSTAT_TYPE_FLOAT,
        READSTAT_TYPE_DOUBLE,
        READSTAT_TYPE_STRING_REF
        
    ctypedef struct readstat_value_t:
        pass
    
    ctypedef struct readstat_parser_t:
        pass


    ctypedef enum readstat_measure_t:
        READSTAT_MEASURE_UNKNOWN,
        READSTAT_MEASURE_NOMINAL = 1,
        READSTAT_MEASURE_ORDINAL,
        READSTAT_MEASURE_SCALE

    ctypedef enum readstat_alignment_t:
        READSTAT_ALIGNMENT_UNKNOWN,
        READSTAT_ALIGNMENT_LEFT = 1,
        READSTAT_ALIGNMENT_CENTER,
        READSTAT_ALIGNMENT_RIGHT

    ctypedef enum readstat_compress_t:
        READSTAT_COMPRESS_NONE,
        READSTAT_COMPRESS_ROWS,
        READSTAT_COMPRESS_BINARY
        

    ctypedef off_t readstat_off_t

    cdef readstat_parser_t *readstat_parser_init()
    cdef void readstat_parser_free(readstat_parser_t *parser)

    ctypedef int (*readstat_open_handler)(const char *path, void *io_ctx);
    ctypedef int (*readstat_close_handler)(void *io_ctx);
    ctypedef readstat_off_t (*readstat_seek_handler)(readstat_off_t offset, readstat_io_flags_t whence, void *io_ctx);
    ctypedef ssize_t (*readstat_read_handler)(void *buf, size_t nbyte, void *io_ctx);
    ctypedef int (*readstat_metadata_handler)(readstat_metadata_t *metadata, void *ctx);
    ctypedef int (*readstat_variable_handler)(int index, readstat_variable_t *variable, char *val_labels, void *ctx);
    ctypedef int (*readstat_value_handler)(int obs_index, readstat_variable_t *variable, readstat_value_t value, void *ctx);
    ctypedef int (*readstat_value_label_handler)(const char *val_labels, readstat_value_t value, const char *label, void *ctx);
    ctypedef int (*readstat_note_handler)(int note_index, const char *note, void *ctx);

    cdef readstat_error_t readstat_set_open_handler(readstat_parser_t *parser, readstat_open_handler open_handler);
    cdef readstat_error_t readstat_set_close_handler(readstat_parser_t *parser, readstat_close_handler close_handler);
    cdef readstat_error_t readstat_set_seek_handler(readstat_parser_t *parser, readstat_seek_handler seek_handler);
    cdef readstat_error_t readstat_set_read_handler(readstat_parser_t *parser, readstat_read_handler read_handler);
    cdef readstat_error_t readstat_set_io_ctx(readstat_parser_t *parser, void *io_ctx);
    cdef readstat_error_t readstat_set_metadata_handler(readstat_parser_t *parser, readstat_metadata_handler metadata_handler);
    cdef readstat_error_t readstat_set_note_handler(readstat_parser_t *parser, readstat_note_handler note_handler);
    cdef readstat_error_t readstat_set_variable_handler(readstat_parser_t *parser, readstat_variable_handler variable_handler)
    cdef readstat_error_t readstat_set_value_handler(readstat_parser_t *parser, readstat_value_handler value_handler);
    cdef readstat_error_t readstat_set_value_label_handler(readstat_parser_t *parser, readstat_value_label_handler value_label_handler);

    cdef readstat_error_t readstat_set_file_character_encoding(readstat_parser_t *parser, const char *encoding);
    
    cdef  int readstat_get_row_count(readstat_metadata_t *metadata);
    cdef int readstat_get_var_count(readstat_metadata_t *metadata);
    cdef char *readstat_get_file_label(readstat_metadata_t *metadata);
    cdef size_t readstat_get_multiple_response_sets_length(readstat_metadata_t *metadata);
    cdef mr_set_t *readstat_get_multiple_response_sets(readstat_metadata_t *metadata);
    cdef char *readstat_get_file_encoding(readstat_metadata_t *metadata);
    cdef char *readstat_get_table_name(readstat_metadata_t *metadata);
    cdef int readstat_get_creation_time(readstat_metadata_t *metadata);
    cdef int readstat_get_modified_time(readstat_metadata_t *metadata);

    cdef int readstat_variable_get_index_after_skipping(readstat_variable_t *variable)
    cdef char *readstat_variable_get_name(readstat_variable_t *variable)
    cdef char *readstat_variable_get_label(readstat_variable_t *variable)
    cdef char *readstat_variable_get_format(readstat_variable_t *variable)
    cdef readstat_type_t readstat_variable_get_type(const readstat_variable_t *variable);
    cdef readstat_measure_t readstat_variable_get_measure(const readstat_variable_t *variable);
    cdef readstat_alignment_t readstat_variable_get_alignment(const readstat_variable_t *variable)
    cdef int readstat_variable_get_display_width(const readstat_variable_t *variable);
    cdef size_t readstat_variable_get_storage_width(const readstat_variable_t *variable)

    cdef int readstat_value_is_missing(readstat_value_t value, readstat_variable_t *variable);
    cdef int readstat_value_is_system_missing(readstat_value_t value);
    cdef int readstat_value_is_tagged_missing(readstat_value_t value);
    cdef int readstat_value_is_defined_missing(readstat_value_t value, readstat_variable_t *variable);
    cdef char readstat_value_tag(readstat_value_t value);

    cdef int readstat_variable_get_missing_ranges_count(const readstat_variable_t *variable);
    cdef readstat_value_t readstat_variable_get_missing_range_lo(const readstat_variable_t *variable, int i);
    cdef readstat_value_t readstat_variable_get_missing_range_hi(const readstat_variable_t *variable, int i);

    cdef readstat_type_t readstat_value_type(readstat_value_t value);

    cdef readstat_error_t readstat_set_row_limit(readstat_parser_t *parser, long row_limit);
    cdef readstat_error_t readstat_set_row_offset(readstat_parser_t *parser, long row_offset);
    
    cdef readstat_error_t readstat_parse_dta(readstat_parser_t *parser, const char *path, void *user_ctx);
    cdef readstat_error_t readstat_parse_sav(readstat_parser_t *parser, const char *path, void *user_ctx);
    cdef readstat_error_t readstat_parse_por(readstat_parser_t *parser, const char *path, void *user_ctx);
    cdef readstat_error_t readstat_parse_sas7bdat(readstat_parser_t *parser, const char *path, void *user_ctx);
    cdef readstat_error_t readstat_parse_sas7bcat(readstat_parser_t *parser, const char *path, void *user_ctx);
    cdef readstat_error_t readstat_parse_xport(readstat_parser_t *parser, const char *path, void *user_ctx);
    
    char *readstat_string_value(readstat_value_t value);
    char readstat_int8_value(readstat_value_t value);
    int16_t readstat_int16_value(readstat_value_t value);
    int32_t readstat_int32_value(readstat_value_t value);
    float readstat_float_value(readstat_value_t value);
    double readstat_double_value(readstat_value_t value);

    char *readstat_error_message(readstat_error_t error_code);

    # Write API
    ctypedef struct readstat_writer_t
    ctypedef ssize_t (*readstat_data_writer)(const void *data, size_t len, void *ctx)
    ctypedef struct readstat_string_ref_t

    cdef readstat_writer_t *readstat_writer_init()
    cdef readstat_error_t readstat_set_data_writer(readstat_writer_t *writer, readstat_data_writer data_writer)
    
    cdef readstat_label_set_t *readstat_add_label_set(readstat_writer_t *writer, readstat_type_t type, const char *name);
    cdef void readstat_label_double_value(readstat_label_set_t *label_set, double value, const char *label);
    cdef void readstat_label_int32_value(readstat_label_set_t *label_set, int32_t value, const char *label);
    cdef void readstat_label_string_value(readstat_label_set_t *label_set, const char *value, const char *label);
    cdef void readstat_label_tagged_value(readstat_label_set_t *label_set, char tag, const char *label);
    
    cdef readstat_error_t readstat_writer_set_file_label(readstat_writer_t *writer, const char *file_label);
    cdef readstat_error_t readstat_writer_set_file_format_version(readstat_writer_t *writer, uint8_t file_format_version)
    cdef readstat_error_t readstat_writer_set_compression(readstat_writer_t *writer, readstat_compress_t compression)

    cdef void readstat_add_note(readstat_writer_t *writer, const char *note);

    cdef readstat_variable_t *readstat_add_variable(readstat_writer_t *writer, const char *name, readstat_type_t type,
        size_t storage_width)

    cdef void readstat_variable_set_label(readstat_variable_t *variable, const char *label)
    cdef void readstat_variable_set_format(readstat_variable_t *variable, const char *format);
    cdef void readstat_variable_set_label_set(readstat_variable_t *variable, readstat_label_set_t *label_set);
    cdef void readstat_variable_set_measure(readstat_variable_t *variable, readstat_measure_t measure);
    cdef void readstat_variable_set_alignment(readstat_variable_t *variable, readstat_alignment_t alignment);
    cdef void readstat_variable_set_display_width(readstat_variable_t *variable, int display_width);
    cdef readstat_error_t readstat_variable_add_missing_double_value(readstat_variable_t *variable, double value);
    cdef readstat_error_t readstat_variable_add_missing_double_range(readstat_variable_t *variable, double lo, double hi);
    cdef readstat_error_t readstat_variable_add_missing_string_value(readstat_variable_t *variable, const char *value);
    cdef readstat_error_t readstat_variable_add_missing_string_range(readstat_variable_t *variable, const char *lo, const char *hi);

    cdef readstat_error_t readstat_writer_set_table_name(readstat_writer_t *writer, const char *table_name)

    cdef readstat_variable_t *readstat_get_variable(readstat_writer_t *writer, int index)

    cdef readstat_error_t readstat_begin_writing_dta(readstat_writer_t *writer, void *user_ctx, long row_count);
    cdef readstat_error_t readstat_begin_writing_por(readstat_writer_t *writer, void *user_ctx, long row_count);
    cdef readstat_error_t readstat_begin_writing_sas7bcat(readstat_writer_t *writer, void *user_ctx);
    cdef readstat_error_t readstat_begin_writing_sas7bdat(readstat_writer_t *writer, void *user_ctx, long row_count);
    cdef readstat_error_t readstat_begin_writing_sav(readstat_writer_t *writer, void *user_ctx, long row_count);
    cdef readstat_error_t readstat_begin_writing_xport(readstat_writer_t *writer, void *user_ctx, long row_count);

    cdef readstat_error_t readstat_validate_metadata(readstat_writer_t *writer)
    cdef readstat_error_t readstat_validate_variable(readstat_writer_t *writer, const readstat_variable_t *variable)

    cdef readstat_error_t readstat_begin_row(readstat_writer_t *writer);

    cdef readstat_error_t readstat_insert_int8_value(readstat_writer_t *writer, const readstat_variable_t *variable, int8_t value);
    cdef readstat_error_t readstat_insert_int16_value(readstat_writer_t *writer, const readstat_variable_t *variable, int16_t value);
    cdef readstat_error_t readstat_insert_int32_value(readstat_writer_t *writer, const readstat_variable_t *variable, int32_t value);
    cdef readstat_error_t readstat_insert_float_value(readstat_writer_t *writer, const readstat_variable_t *variable, float value);
    cdef readstat_error_t readstat_insert_double_value(readstat_writer_t *writer, const readstat_variable_t *variable, double value);
    cdef readstat_error_t readstat_insert_string_value(readstat_writer_t *writer, const readstat_variable_t *variable, const char *value);
    cdef readstat_error_t readstat_insert_missing_value(readstat_writer_t *writer, const readstat_variable_t *variable);
    cdef readstat_error_t readstat_insert_tagged_missing_value(readstat_writer_t *writer, const readstat_variable_t *variable, char tag);
    # String refs are used for creating a READSTAT_TYPE_STRING_REF column,
    # which is only supported in Stata. String references can be shared
    # across columns, and inserted with readstat_insert_string_ref().
    cdef readstat_error_t readstat_insert_string_ref(readstat_writer_t *writer, const readstat_variable_t *variable, readstat_string_ref_t *ref);
    cdef readstat_string_ref_t *readstat_add_string_ref(readstat_writer_t *writer, const char *string);
    cdef readstat_string_ref_t *readstat_get_string_ref(readstat_writer_t *writer, int index);

    cdef readstat_error_t readstat_end_row(readstat_writer_t *writer);

    cdef readstat_error_t readstat_end_writing(readstat_writer_t *writer);
    cdef void readstat_writer_free(readstat_writer_t *writer);

cdef extern from "readstat_io_unistd.h":
    cdef struct unistd_io_ctx_t "unistd_io_ctx_s":
        int fd
