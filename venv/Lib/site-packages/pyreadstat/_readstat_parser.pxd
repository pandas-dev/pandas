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

# Definitions of enum types
ctypedef enum py_file_extension:
    FILE_EXT_SAV
    FILE_EXT_SAS7BDAT
    FILE_EXT_DTA
    FILE_EXT_XPORT
    FILE_EXT_POR
    FILE_EXT_SAS7BCAT

ctypedef enum py_file_format:
    FILE_FORMAT_SAS
    FILE_FORMAT_SPSS
    FILE_FORMAT_STATA
    
ctypedef enum py_datetime_format:
    DATE_FORMAT_NOTADATE
    DATE_FORMAT_DATE
    DATE_FORMAT_DATETIME
    DATE_FORMAT_TIME
    
ctypedef enum py_variable_format:
    VAR_FORMAT_STRING
    VAR_FORMAT_LONG
    VAR_FORMAT_FLOAT
    VAR_FORMAT_MISSING
    
# Definitions of extension types
    
cdef class data_container:
    """
    This extension type holds all the data we need to get from our file from 
    the readstat parser and will be used later to compose our pandas data frame
    """
    cdef int n_obs
    cdef int n_vars
    cdef int max_n_obs
    cdef list col_data
    cdef list col_data_len
    cdef list col_names
    cdef list col_labels
    cdef list col_dtypes
    cdef list col_numpy_dtypes
    cdef list col_dtypes_isobject
    cdef list col_dytpes_isfloat
    cdef list col_formats
    cdef list col_formats_original
    cdef object origin
    cdef double unix_to_origin_secs
    cdef py_file_format file_format
    cdef bint is_unkown_number_rows
    cdef str file_label
    cdef str file_encoding
    cdef bint metaonly
    cdef bint dates_as_pandas
    cdef object label_to_var_name
    cdef object labels_raw
    cdef object notes
    cdef str user_encoding
    cdef str table_name
    cdef bint filter_cols
    cdef list use_cols
    cdef bint usernan
    cdef dict missing_ranges
    cdef dict missing_user_values
    cdef dict variable_storage_width
    cdef variable_display_width
    cdef variable_alignment
    cdef variable_measure
    cdef bint no_datetime_conversion
    cdef int ctime
    cdef int mtime
    cdef dict mr_sets
    cdef str output_format
    
cdef dict readstat_to_numpy_types

# definitions of functions
cdef py_datetime_format transform_variable_format(str var_format, py_file_format file_format)
cdef object transform_datetime(py_datetime_format var_format, double tstamp, py_file_format file_format, object origin,
                               bint dates_as_pandas, str output_format, double unix_to_origin_secs)

cdef int handle_metadata(readstat_metadata_t *metadata, void *ctx) except READSTAT_HANDLER_ABORT
cdef int handle_variable(int index, readstat_variable_t *variable, 
                         char *val_labels, void *ctx) except READSTAT_HANDLER_ABORT
cdef int handle_value(int obs_index, readstat_variable_t * variable, readstat_value_t value, void *ctx) except READSTAT_HANDLER_ABORT
cdef int handle_value_label(char *val_labels, readstat_value_t value, char *label, void *ctx) except READSTAT_HANDLER_ABORT
cdef int handle_note (int note_index, char *note, void *ctx) except READSTAT_HANDLER_ABORT

cdef void check_exit_status(readstat_error_t retcode) except *

cdef void run_readstat_parser(char * filename, data_container data, py_file_extension file_extension, long row_limit, long row_offset) except *
cdef object data_container_to_dict(data_container data)
cdef object dict_to_dataframe(object dict_data, data_container dc)
cdef object data_container_extract_metadata(data_container data)
cdef object run_conversion(object filename_path, py_file_format file_format, py_file_extension file_extension,
                           str encoding, bint metaonly, bint dates_as_pandas, list usecols, bint usernan,
                           bint no_datetime_conversion, long row_limit, long row_offset, str output_format, list extra_datetime_formats, 
			   list extra_date_formats, list extra_time_formats)

# definitions for stuff about dates
cdef list sas_date_formats 
cdef list sas_datetime_formats 
cdef list sas_time_formats 
cdef list sas_all_formats 
cdef object sas_origin

cdef list spss_datetime_formats 
cdef list spss_date_formats 
cdef list spss_time_formats 
cdef list spss_all_formats 
cdef object spss_origin

cdef list stata_datetime_formats
cdef list stata_date_formats
cdef list stata_time_formats 
cdef list stata_all_formats 
cdef object stata_origin

# Stuff for opening files on windows in order to handle international characters
# Courtesy of Jonathon Love
# works only in python 3
#cdef extern from "readstat_io_unistd.h":
#    cdef struct unistd_io_ctx_t "unistd_io_ctx_s":
#        int fd
    
cdef extern from "Python.h":
    wchar_t* PyUnicode_AsWideCharString(object, Py_ssize_t *) except NULL
    
# these ones would make the c file produced by cython not portable between windows and unix
# therefore the conditional including of the libraries is handled in C
cdef extern from "conditional_includes.h":
    int _wsopen(const wchar_t *filename, int oflag, int shflag, int pmode)
    int _O_RDONLY
    int _O_BINARY
    int _O_WRONLY
    int _O_CREAT
    int _SH_DENYRW  # Denies read and write access to a file.
    int _SH_DENYWR  # Denies write access to a file.
    int _SH_DENYRD  # Denies read access to a file.
    int _SH_DENYNO
    void assign_fd(void *io_ctx, int fd)
    long seek_fd(readstat_off_t offset, readstat_io_flags_t whence, void *io_ctx)

