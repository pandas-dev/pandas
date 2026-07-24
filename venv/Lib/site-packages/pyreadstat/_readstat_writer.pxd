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
from libc.stddef cimport wchar_t

from readstat_api cimport *
#cdef int write_test()

cdef extern from "readstat_io_unistd.h":
    cdef struct unistd_io_ctx_t "unistd_io_ctx_s":
        pass

cdef extern from "conditional_includes.h":
    int _wsopen(const wchar_t *filename, int oflag, int shflag, int pmode)
    int O_RDONLY
    int O_BINARY
    int O_WRONLY
    int O_CREAT
    int O_TRUNC
    int _O_RDONLY
    int _O_BINARY
    int _O_WRONLY
    int _O_CREAT
    int _O_TRUNC
    int _SH_DENYRW  # Denies read and write access to a file.
    int _SH_DENYWR  # Denies write access to a file.
    int _SH_DENYRD  # Denies read access to a file.
    int _SH_DENYNO
    int _S_IWRITE
    int _S_IREAD
    int open(const char *path, int oflag, int mode)
    int _close(int fd)

cdef extern from "Python.h":
    wchar_t* PyUnicode_AsWideCharString(object, Py_ssize_t *) except NULL

cdef extern from "conditional_includes.h":
    int _close(int fd) 
    ssize_t _write(int fd, const void *buf, size_t nbyte)
    int close(int fd)
    ssize_t write(int fd, const void *buf, size_t nbyte)

ctypedef enum dst_file_format:
    FILE_FORMAT_SAS7BDAT
    FILE_FORMAT_SAS7BCAT
    FILE_FORMAT_XPORT
    FILE_FORMAT_SAV
    FILE_FORMAT_DTA
    FILE_FORMAT_POR

ctypedef enum pywriter_variable_type:
    PYWRITER_DOUBLE
    PYWRITER_INTEGER
    PYWRITER_CHARACTER
    PYWRITER_LOGICAL
    PYWRITER_OBJECT
    PYWRITER_DATE
    PYWRITER_DATETIME
    PYWRITER_TIME
    PYWRITER_DATE64
    PYWRITER_DATETIME64
    PYWRITER_TIME64
    PYWRITER_DTA_STR_REF

cdef object vectorized_convert_datetime_to_number(object df, dst_file_format file_format, list pywriter_types, list pywriter_timeunits, int col_count)
cdef object vectorized_convert_date_to_number(object df, dst_file_format file_format, list pywriter_types,  int col_count)
cdef object vectorized_convert_time_to_number(object df, dst_file_format file_format, list pywriter_types,  int col_count)
cdef double convert_datetimelike_to_number(dst_file_format file_format, pywriter_variable_type curtype, object curval) except *
cdef char * get_datetimelike_format_for_readstat(dst_file_format file_format, pywriter_variable_type curtype)
cdef int get_narwhals_str_series_max_length(object series, dict value_labels, bint isobject)
cdef int check_series_all_same_types(object series, object type_to_check)
cdef list get_narwhals_column_types(object df, dict missing_user_values, dict variable_value_labels, int dta_str_max_len)
cdef ssize_t write_bytes(const void *data, size_t _len, void *ctx)
#cdef void check_exit_status(readstat_error_t retcode) except *
cdef int open_file(bytes filename_path)
cdef int close_file(int fd)
cdef bytes filepath_to_bytes(object filename_path)
cdef void initial_checks(bint is_pandas, bint is_polars, dict variable_value_labels, dict missing_user_values,
                        dst_file_format file_format, list col_names, bytes filename_bytes) except *

cdef int run_write(df, object filename_path, dst_file_format file_format, str file_label, object column_labels,
                   int file_format_version, object note, str table_name, dict variable_value_labels, 
                   dict missing_ranges, dict missing_user_values, dict variable_alignment,
                   dict variable_display_width, dict variable_measure, dict variable_format, bint row_compression) except *
