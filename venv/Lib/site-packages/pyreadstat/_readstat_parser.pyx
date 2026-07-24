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

## if want to profile: # cython: profile=True

from cpython.datetime cimport import_datetime, timedelta_new, datetime_new, total_seconds
from cpython.exc cimport PyErr_Occurred
from cpython.object cimport PyObject
from libc.math cimport floor
from libc.string cimport memcpy

from collections import OrderedDict
import datetime
import os
import warnings
import sys
import re

import narwhals.stable.v2 as nw
import numpy as np

from readstat_api cimport *

from pyclasses import metadata_container

# necessary to work with the datetime C API
import_datetime()

cdef object unix_origin = datetime_new(1970, 1, 1, 0, 0, 0, 0, None)

cdef object format_regex = re.compile(r"^([A-Z][A-Z0-9]+[A-Z])(\d+)?(?(2)(?:\.\d+)?$|$)")
cdef list sas_date_formats = ["WEEKDATE", "MMDDYY", "DDMMYY", "YYMMDD", "DATE", "YYMMDD",
                                "DDMMYYB", "DDMMYYC", "DDMMYYD",
                                "DDMMYYN", "DDMMYYP", "DDMMYYS",
                                "MMDDYYB", "MMDDYYC", "MMDDYYD",
                                "MMDDYYN", "MMDDYYP", "MMDDYYS",
                                #"MONNAME", "MONTH",  "WEEKDAY", "QTR", "QTRR", "YEAR","DAY",  "DOWNAME" # these do not print as full dates in sas
                                "WEEKDATX", "DTDATE",
                                "IS8601DA", "E8601DA", "B8601DA",
                                "YYMMDDB", "YYMMDDD", "YYMMDDN", "YYMMDDP", "YYMMDDS"]
cdef list sas_datetime_formats = ["DATETIME", "E8601DT", "DATEAMPM", "MDYAMPM", "IS8601DT", "B8601DT", "B8601DN"]
cdef list sas_time_formats = ["TIME", "HHMM", "TOD", "TIMEAMPM", "IS8601TM", "E8601TM", "B8601TM"]
# "HOUR" # these do not print as full time formats in sas 
#cdef list sas_all_formats = sas_date_formats + sas_datetime_formats + sas_time_formats
cdef list sas_all_formats
cdef object sas_origin = datetime_new(1960, 1, 1, 0, 0, 0, 0, None)
cdef object sas_secs_from_unix = total_seconds(unix_origin - sas_origin)

cdef list spss_datetime_formats = ["DATETIME", "YMDHMS"]
cdef list spss_date_formats = ["DATE", "ADATE", "EDATE", "JDATE", "SDATE"]
cdef list spss_time_formats = ["TIME", "DTIME"]
#cdef list spss_all_formats = spss_date_formats + spss_datetime_formats + spss_time_formats
cdef list spss_all_formats
cdef object spss_origin = datetime_new(1582, 10, 14, 0, 0, 0, 0, None)
cdef object spss_secs_from_unix = total_seconds(unix_origin - spss_origin)

cdef list stata_datetime_formats = ["%tC", "%tc"]
cdef list stata_date_formats = ["%td", "%d", "%tdD_m_Y", "%tdCCYY-NN-DD"]
cdef list stata_time_formats = ["%tcHH:MM:SS", "%tcHH:MM"]
#cdef list stata_all_formats = stata_datetime_formats + stata_date_formats + stata_time_formats
cdef list stata_all_formats
cdef object stata_origin = datetime_new(1960, 1, 1, 0, 0, 0, 0, None)
cdef object stata_secs_from_unix = total_seconds(unix_origin - stata_origin)

cdef dict readstat_to_numpy_types = {READSTAT_TYPE_STRING: object, READSTAT_TYPE_STRING_REF: object,
                                     READSTAT_TYPE_INT8: np.int64, READSTAT_TYPE_INT16: np.int64, READSTAT_TYPE_INT32:np.int64,
                                     READSTAT_TYPE_FLOAT: np.float64, READSTAT_TYPE_DOUBLE: np.float64}

cdef class data_container:
    """
    This extension type holds all the data we need to get from our file from 
    the readstat parser and will be used later to compose our pandas data frame
    """
    
    def __cinit__(self):
        self.n_obs = 0
        self.n_vars = 0
        self.max_n_obs = 0
        self.col_data = list()
        self.col_data_len = list()
        self.col_names = list()
        self.col_labels = list()
        self.col_dtypes = list()
        self.col_numpy_dtypes = list()
        self.col_dtypes_isobject = list()
        self.col_dytpes_isfloat = list()
        self.col_formats = list()
        self.col_formats_original = list()
        self.origin = None
        self.unix_to_origin_secs = 0
        self.is_unkown_number_rows = 0
        self.file_encoding = None
        self.file_label = None
        self.metaonly = 0
        self.dates_as_pandas = 0
        self.label_to_var_name = dict()
        self.labels_raw = dict()
        self.notes = list()
        self.user_encoding = None
        self.table_name = None
        self.filter_cols = 0
        self.use_cols = list()
        self.usernan = 0
        self.missing_ranges = dict()
        self.missing_user_values = dict()
        self.variable_storage_width = dict()
        self.variable_display_width = dict()
        self.variable_alignment = dict()
        self.variable_measure = dict()
        self.no_datetime_conversion = 0
        self.ctime = 0
        self.mtime = 0
        self.mr_sets = dict()
        self.output_format = ""


class ReadstatError(Exception):
    """
    Just defining a custom exception to raise when readstat gives an error return code.
    """
    pass


class PyreadstatError(Exception):
    """
    Just defining a custom exception to raise when pyreadstat raises an exception.
    """
    pass



cdef py_datetime_format transform_variable_format(str var_format, py_file_format file_format):
    """
    Transforms a readstat var_format to a date, datetime or time format label
    """
    if file_format == FILE_FORMAT_SAS:
        if var_format:
            format_match = format_regex.match(var_format)
            if format_match:
                var_format_name = format_match.group(1)
                if var_format_name in sas_all_formats:
                    if var_format_name in sas_date_formats:
                        return DATE_FORMAT_DATE
                    elif var_format_name in sas_datetime_formats:
                        return DATE_FORMAT_DATETIME
                    elif var_format_name in sas_time_formats:
                        return DATE_FORMAT_TIME
        return DATE_FORMAT_NOTADATE
        
    elif file_format == FILE_FORMAT_SPSS:
        if var_format:
            format_match = format_regex.match(var_format)
            if format_match:
                var_format_name = format_match.group(1)
                if var_format_name in spss_all_formats:
                    if var_format_name in spss_date_formats:
                        return DATE_FORMAT_DATE
                    elif var_format_name in spss_datetime_formats:
                        return DATE_FORMAT_DATETIME
                    elif var_format_name in spss_time_formats:
                        return DATE_FORMAT_TIME
        return DATE_FORMAT_NOTADATE
        
    elif file_format == FILE_FORMAT_STATA:
        if var_format in stata_all_formats:
            if var_format in stata_date_formats:
                return DATE_FORMAT_DATE
            elif var_format in stata_datetime_formats:
                return DATE_FORMAT_DATETIME
            elif var_format in stata_time_formats:
                return DATE_FORMAT_TIME
        else:
            return DATE_FORMAT_NOTADATE

cdef object transform_datetime(py_datetime_format var_format, double tstamp, py_file_format file_format, object origin,
                               bint dates_as_pandas, str output_format, double unix_to_origin_secs):
    """
    Transforms a tstamp integer value to a date, time or datetime pyton object.
    tstamp could represent number of days, seconds or milliseconds
    """
    
    cdef object tdelta
    cdef int days
    cdef int secs
    cdef double msecs
    cdef int usecs
    cdef object mydat

    # For polars we are going to return an epoch from unix origin, 
    # later we will transform that to datetime or date with polars from_epoch
    # for any other we return a python datetime, date or time object. 
    # polars also works with these objects, but it is slower
    if var_format == DATE_FORMAT_DATE:
        if output_format == "polars":
            # we want to return days from unix
            if file_format == FILE_FORMAT_SPSS:
                # tstamp is in seconds
                return (tstamp - unix_to_origin_secs)/86400
            else:
                # tstamp is in days
                return tstamp - (unix_to_origin_secs/86400)

        if file_format == FILE_FORMAT_SPSS:
            # tstamp is in seconds
            days = <int> (floor(tstamp / 86400))
            secs = <int> (tstamp % 86400)
            tdelta = timedelta_new(days, secs, 0)
            #tdelta = timedelta(seconds=tstamp)
        else:
            # tstamp is in days
            days = <int> tstamp
            tdelta = timedelta_new(days, 0, 0)
            #tdelta = timedelta(days=tstamp)
        mydat = origin + tdelta
        if dates_as_pandas:
            return mydat
        else:
            return mydat.date()
    elif var_format == DATE_FORMAT_DATETIME:
        if output_format == "polars":
            # we want to return timestamp in seconds
            if file_format == FILE_FORMAT_STATA:
                # tstamp is in millisecons
                return (tstamp/1000)
            else:
                # tstamp in seconds
                return tstamp

        if file_format == FILE_FORMAT_STATA:
            # tstamp is in millisecons
            days = <int> (floor(tstamp / 86400000))
            msecs = tstamp % 86400000
            secs = <int> (msecs/1000)
            usecs = <int> ((msecs % 1000) * 1000 )
            tdelta = timedelta_new(days, secs, usecs)
            #tdelta = timedelta(milliseconds=tstamp)
        else:
            # tstamp in seconds
            days = <int> (floor(tstamp / 86400))
            secs = <int> (tstamp % 86400)
            usecs = <int> (round(tstamp % 1 * 1e6))
            tdelta = timedelta_new(days, secs, usecs)
            #tdelta = timedelta(seconds=tstamp)
        mydat = origin + tdelta
        return mydat
    elif var_format == DATE_FORMAT_TIME:
        if file_format == FILE_FORMAT_STATA:
            # tstamp is in millisecons
            days = <int> (floor(tstamp / 86400000))
            msecs = tstamp % 86400000
            secs = <int> (msecs/1000)
            usecs = <int> ((msecs % 1000) * 1000 )
            tdelta = timedelta_new(days, secs, usecs)
            #tdelta = timedelta(milliseconds=tstamp)
        else:
            # tstamp in seconds
            days = <int> (floor(tstamp / 86400))
            secs = <int> (tstamp % 86400)
            usecs = <int> (round(tstamp % 1 * 1e6))
            tdelta = timedelta_new(days, secs, usecs)
            #tdelta = timedelta(seconds=tstamp)
        mydat = origin + tdelta
        return mydat.time()


cdef object convert_readstat_to_python_value(readstat_value_t value, int index, data_container dc):
    """
    Converts a readstat value to a python value. 
    """

    #py_file_format file_format, py_datetime_format var_format

    cdef readstat_type_t var_type
    cdef py_datetime_format var_format
    cdef py_variable_format pyformat
    cdef object origin
    cdef bint dates_as_pandas
    cdef py_file_format file_format
    cdef object result

    cdef char * c_str_value
    cdef str py_str_value
    cdef int8_t c_int8_value
    cdef int16_t c_int16_value
    cdef int32_t c_int32_value
    cdef int64_t c_int64_value
    cdef float c_float_value
    cdef double c_double_value
    cdef long py_long_value
    cdef double py_float_value
    cdef double tstamp
    cdef str output_format

    var_type = dc.col_dtypes[index]
    var_format = dc.col_formats[index]
    origin = dc.origin
    unix_to_origin_secs = dc.unix_to_origin_secs
    dates_as_pandas = dc.dates_as_pandas
    file_format = dc.file_format
    output_format = dc.output_format

    # transform to values cython can deal with
    if var_type == READSTAT_TYPE_STRING or var_type == READSTAT_TYPE_STRING_REF:
        c_str_value = readstat_string_value(value)
        if c_str_value != NULL:
            py_str_value = <str> c_str_value
        else:
            py_str_value = ''
        pyformat = VAR_FORMAT_STRING
    elif var_type == READSTAT_TYPE_INT8:
        c_int8_value = readstat_int8_value(value)
        py_long_value = <long> c_int8_value
        pyformat = VAR_FORMAT_LONG
    elif var_type == READSTAT_TYPE_INT16:
        c_int16_value = readstat_int16_value(value)
        py_long_value = <long> c_int16_value
        pyformat = VAR_FORMAT_LONG
    elif var_type == READSTAT_TYPE_INT32:
        c_int32_value = readstat_int32_value(value)
        py_long_value = <long> c_int32_value
        pyformat = VAR_FORMAT_LONG
    elif var_type == READSTAT_TYPE_FLOAT:
        c_float_value = readstat_float_value(value)
        py_float_value = <double> c_float_value
        pyformat = VAR_FORMAT_FLOAT
    elif var_type == READSTAT_TYPE_DOUBLE:
        c_double_value = readstat_double_value(value);
        py_float_value = <double> c_double_value
        pyformat = VAR_FORMAT_FLOAT
    else:
        raise PyreadstatError("Unkown data type")

    # final transformation and storage

    if pyformat == VAR_FORMAT_STRING:
        if var_format == DATE_FORMAT_NOTADATE or dc.no_datetime_conversion:
            result = py_str_value
        else:
            #str_byte_val = py_str_value.encode("UTF-8")
            raise PyreadstatError("STRING type with value '%s' with date type in column '%s'" % (py_str_value, dc.col_names[index]))
    elif pyformat == VAR_FORMAT_LONG:
        if var_format == DATE_FORMAT_NOTADATE or dc.no_datetime_conversion:
            result = py_long_value
        else:
            tstamp = <double> py_long_value
            result = transform_datetime(var_format, tstamp, file_format, origin, dates_as_pandas, output_format, unix_to_origin_secs)
    elif pyformat == VAR_FORMAT_FLOAT:
        if var_format == DATE_FORMAT_NOTADATE or dc.no_datetime_conversion:
            result = py_float_value
        else:
            #tstamp = <int> py_float_value
            tstamp = py_float_value
            result = transform_datetime(var_format, tstamp, file_format, origin, dates_as_pandas, output_format, unix_to_origin_secs)
    #elif pyformat == VAR_FORMAT_MISSING:
    #    pass
    else:
        raise PyreadstatError("Failed convert C to python value")

    return result


cdef int handle_metadata(readstat_metadata_t *metadata, void *ctx) except READSTAT_HANDLER_ABORT:
    """
    This function sets the number of observations(rows), number of variables
    (columns) in data container and initializes the col_data which will store
    all the data we need.
    It also reads other metadata from the file.
    """
        
    cdef int var_count, obs_count, mr_len
    cdef  data_container dc = <data_container> ctx
    #cdef object row
    cdef char * flabel_orig
    cdef char * fencoding_orig
    cdef str flabel, fencoding
    cdef bint metaonly
    cdef char * table
    cdef int ctime
    cdef int mtime
    cdef int i = 0
    cdef mr_set_t * mr_sets_orig
    cdef dict mr_sets = {}
    cdef str name
    cdef list variable_list = []

    metaonly = dc.metaonly
    
    var_count = readstat_get_var_count(metadata)
    if var_count<0:
        raise PyreadstatError("Failed to read number of variables")
    obs_count = readstat_get_row_count(metadata)
    if obs_count <0:
        # if <0 it means the number of rows is not known, allocate 100 000
        obs_count = 100000
        dc.is_unkown_number_rows = 1
    
    dc.n_obs = obs_count
    dc.n_vars = var_count

    mr_len = <int>readstat_get_multiple_response_sets_length(metadata);
    mr_sets_orig = readstat_get_multiple_response_sets(metadata);
    if mr_sets_orig != NULL:
        i = 0
        while i < mr_len:
            name = <str>mr_sets_orig[i].name
            variable_list = []
            for j in range(mr_sets_orig[i].num_subvars):
                variable_list.append(<str>mr_sets_orig[i].subvariables[j])
            mr_sets[name] = {
                'type': chr(mr_sets_orig[i].type),
                'is_dichotomy': bool(mr_sets_orig[i].is_dichotomy),
                'counted_value': mr_sets_orig[i].counted_value if mr_sets_orig[i].counted_value != -1 else None,
                'label': <str>mr_sets_orig[i].label,
                'variable_list': variable_list
            } 
            i += 1
        dc.mr_sets = mr_sets
    
    dc.col_data_len = [obs_count] * var_count
    dc.col_numpy_dtypes = [None] * var_count
    dc.col_dytpes_isfloat = [0] * var_count
    dc.col_dtypes_isobject = [0] * var_count
    
    # read other metadata
    flabel_orig = readstat_get_file_label(metadata);
    fencoding_orig = readstat_get_file_encoding(metadata)
    if flabel_orig != NULL and flabel_orig[0]:
        flabel = <str> flabel_orig
    else:
        flabel = None
    if fencoding_orig != NULL and fencoding_orig[0]:
        fencoding = <str> fencoding_orig
    else:
        fencoding = None
        
    dc.file_encoding = fencoding
    dc.file_label = flabel

    table = readstat_get_table_name(metadata)
    if table != NULL and table[0]:
        dc.table_name = <str> table

    dc.ctime = readstat_get_creation_time(metadata)
    dc.mtime = readstat_get_modified_time(metadata)

    return READSTAT_HANDLER_OK

cdef int handle_variable(int index, readstat_variable_t *variable, 
                         char *val_labels, void *ctx) except READSTAT_HANDLER_ABORT:
    """
    This function extracts the name, label, type and format from a variable
    and stores it in data container for a later access. It also extracts the label set to which the variable is associated,
    if any.
    """

    cdef char * var_name, 
    cdef char * var_label
    cdef char * var_format
    cdef str col_name, col_label, label_name, col_format_original, output_format
    cdef py_datetime_format col_format_final
    cdef readstat_type_t var_type
    cdef py_file_format file_format
    cdef readstat_label_set_t * lset
    cdef int n_ranges
    cdef readstat_value_t loval, hival
    cdef object pyloval, pyhival
    cdef list missing_ranges
    cdef object curnptype
    cdef object row
    cdef bint metaonly
    cdef int obs_count
    cdef bint iscurnptypefloat
    cdef bint iscurnptypeobject
    cdef str newcolname
    cdef int dupcolcnt

    cdef  data_container dc = <data_container> ctx
    output_format = dc.output_format
    
    # get variable name, label, format and type and put into our data container
    var_name = readstat_variable_get_name(variable)
    if var_name == NULL:
        col_name = None
    else:
        col_name = <str>var_name

    # if the user introduced a list of columns to include, continue only if the column is in the list
    if dc.filter_cols and not (col_name in dc.use_cols):
        dc.n_vars -= 1
        return READSTAT_HANDLER_SKIP_VARIABLE

    index = readstat_variable_get_index_after_skipping(variable)

    if col_name in dc.col_names:
        dupcolcnt = 1
        while True:
            newcolname = col_name + "_duplicated" + str(dupcolcnt)
            if newcolname in col_name:
                dupcolcnt += 1
                continue
            else:
                msg = "column '{0}' is duplicated, renamed to '{1}'".format(col_name, newcolname)
                warnings.warn(msg)
                col_name = newcolname
                break

    dc.col_names.append(col_name)

    # the name of the value label for the variable
    if val_labels != NULL:
        label_name = <str> val_labels
        if label_name:
            dc.label_to_var_name[col_name] = label_name
    
    var_label = readstat_variable_get_label(variable)
    if var_label == NULL:
        col_label = None
    else:
        col_label = <str>var_label
    dc.col_labels.append(col_label)

    # format, we have to transform it in something more usable
    var_format = readstat_variable_get_format(variable)
    if var_format == NULL:
        #col_format_original = "NULL"
        col_format_original = None
    else:
        col_format_original = <str>var_format
    file_format = dc.file_format
    dc.col_formats_original.append(col_format_original)
    col_format_final = transform_variable_format(col_format_original, file_format)
    dc.col_formats.append(col_format_final)
    # readstat type
    var_type = readstat_variable_get_type(variable)
    dc.col_dtypes.append(var_type)
    # equivalent numpy type
    # if it's a date then we need object
    if col_format_final != DATE_FORMAT_NOTADATE and dc.no_datetime_conversion == 0: 
        curnptype = object
    else:
        curnptype = readstat_to_numpy_types[var_type]
    iscurnptypefloat = 0
    iscurnptypeobject = 0
    # book keeping numpy types
    dc.col_numpy_dtypes[index] = curnptype
    if curnptype == object:
        iscurnptypeobject = 1
    if curnptype == np.float64:
        iscurnptypefloat = 1
    dc.col_dtypes_isobject[index] = iscurnptypeobject
    dc.col_dytpes_isfloat[index] = iscurnptypefloat
    metaonly = dc.metaonly
    # pre-allocate data
    if metaonly:
        if output_format == "pandas":
            row = np.empty(1, dtype=curnptype)
        else:
            row = list()
    else:
        obs_count = dc.n_obs
        if output_format == "pandas":
            row = np.empty(obs_count, dtype=curnptype)
            if iscurnptypeobject or iscurnptypefloat:
                row.fill(np.nan)
        else:
            row = [None] * obs_count
    dc.col_data.append(row)
    
    # missing values
    if dc.usernan:
        n_ranges = readstat_variable_get_missing_ranges_count(variable)
        if n_ranges>0:
            missing_ranges = list()
            for i in range(0, n_ranges):
                loval = readstat_variable_get_missing_range_lo(variable, i)
                pyloval = convert_readstat_to_python_value(loval, index, dc)
                hival = readstat_variable_get_missing_range_hi(variable, i)
                pyhival = convert_readstat_to_python_value(hival, index, dc)
                missing_ranges.append({'lo':pyloval, 'hi':pyhival})
            dc.missing_ranges[col_name] = missing_ranges

    cdef size_t storage_width
    storage_width = readstat_variable_get_storage_width(variable)
    dc.variable_storage_width[col_name] = <int>storage_width

    dc.variable_display_width[col_name] = readstat_variable_get_display_width(variable)

    cdef readstat_alignment_t align
    cdef str pyalign
    align = readstat_variable_get_alignment(variable)
    if align == READSTAT_ALIGNMENT_UNKNOWN:
        pyalign = "unknown"
    elif align == READSTAT_ALIGNMENT_LEFT:
        pyalign = "left"
    elif align == READSTAT_ALIGNMENT_CENTER:
        pyalign = "center"
    elif align == READSTAT_ALIGNMENT_RIGHT:
        pyalign = "right"
    else:
        pyalign = "undetermined"

    dc.variable_alignment[col_name] = pyalign

    cdef readstat_measure_t measure
    cdef str pymeasure
    measure = readstat_variable_get_measure(variable)
    if measure == READSTAT_MEASURE_UNKNOWN:
        pymeasure = "unknown"
    elif measure == READSTAT_MEASURE_NOMINAL:
        pymeasure = "nominal"
    elif measure == READSTAT_MEASURE_ORDINAL:
        pymeasure = "ordinal"
    elif measure == READSTAT_MEASURE_SCALE:
        pymeasure = "scale"
    else:
        pymeasure = "undetermined"

    dc.variable_measure[col_name] = pymeasure

    return READSTAT_HANDLER_OK


cdef int handle_value(int obs_index, readstat_variable_t * variable, readstat_value_t value, void *ctx) except READSTAT_HANDLER_ABORT:
    """
    This function transforms every value to python types, and to datetime if 
    necessary and stores it into the arrays in data container col_data
    """

    cdef data_container dc
    cdef int index
    cdef int max_n_obs
    cdef bint is_unkown_number_rows
    cdef int var_max_rows
    cdef object buf_list
    cdef bint iscurnptypeobject
    cdef bint iscurnptypefloat
    cdef int missing_tag
    cdef object pyvalue
    cdef set curset
    cdef object curnptype
    cdef str output_format
    
    # extract variables we need from data container
    dc = <data_container> ctx
    output_format = dc.output_format
    index = readstat_variable_get_index_after_skipping(variable)
    max_n_obs = dc.max_n_obs
    is_unkown_number_rows = dc.is_unkown_number_rows
    
    # check that we still have enough room in our pre-allocated lists
    # if not, add more room
    iscurnptypeobject = dc.col_dtypes_isobject[index]
    iscurnptypefloat = dc.col_dytpes_isfloat[index]
    if is_unkown_number_rows:
        if max_n_obs <= obs_index:
            dc.max_n_obs = obs_index + 1
        var_max_rows = dc.col_data_len[index]
        if var_max_rows <= obs_index:
            if output_format == "pandas":
                curnptype = dc.col_numpy_dtypes[index]
                buf_list = np.empty(100000, dtype=curnptype)
                if iscurnptypeobject or iscurnptypefloat:
                    buf_list.fill(np.nan)
                dc.col_data[index] = np.append(dc.col_data[index], buf_list)
            else:
                buf_list = [None] * 100000
                dc.col_data[index].extend(buf_list) 
            var_max_rows += 100000
            dc.col_data_len[index] = var_max_rows

    # transform to python value types
    if readstat_value_is_missing(value, variable):
        # The user does not want to retrieve missing values
        if not dc.usernan or readstat_value_is_system_missing(value):
            if output_format == "pandas":
                if iscurnptypefloat == 1 or iscurnptypeobject == 1: 
                    # already allocated
                    pass
                    #dc.col_data[index][obs_index] = np.nan
                # for any type except float, the numpy type will be object as now we have nans
                else:
                    dc.col_numpy_dtypes[index] = object
                    dc.col_dtypes_isobject[index] = 1
                    iscurnptypeobject = 1
                    dc.col_data[index] = dc.col_data[index].astype(object, copy=False)
                    dc.col_data[index][obs_index:] = np.nan
                    #dc.col_data[index][obs_index] = NAN
        elif readstat_value_is_defined_missing(value, variable):
            # SPSS missing values
            pyvalue = convert_readstat_to_python_value(value, index, dc)
            dc.col_data[index][obs_index] = pyvalue
        elif readstat_value_is_tagged_missing(value):
            iscurnptypeobject = dc.col_dtypes_isobject[index]
            # SAS and Stata missing values
            missing_tag = <int> readstat_value_tag(value)
            # In SAS missing values are A to Z or _ in stata a to z
            # if (missing_tag >=65 and missing_tag <= 90) or missing_tag == 95 or (missing_tag >=61 and missing_tag <= 122):
            if output_format == "pandas":
                if iscurnptypeobject == 1:
                    dc.col_data[index][obs_index] =  chr(missing_tag) 
                else:
                    dc.col_numpy_dtypes[index] = object
                    dc.col_dtypes_isobject[index] = 1
                    dc.col_dytpes_isfloat[index] = 0
                    iscurnptypeobject = 1
                    dc.col_data[index] = dc.col_data[index].astype(object, copy=False)
                    dc.col_data[index][obs_index] =  chr(missing_tag)
            else:
                dc.col_data[index][obs_index] =  chr(missing_tag)
            curset = dc.missing_user_values.get(index)
            if curset is None:
                curset = set()
            curset.add(chr(missing_tag))
            dc.missing_user_values[index] = curset
    else:
        pyvalue = convert_readstat_to_python_value(value, index, dc)
        dc.col_data[index][obs_index] = pyvalue
        
    return READSTAT_HANDLER_OK


cdef int handle_value_label(char *val_labels, readstat_value_t value, char *label, void *ctx) except READSTAT_HANDLER_ABORT:
    """
    Reads the label for the value that belongs to the label set val_labels. In Handle variable we need to do a map
    from variable name to val_label so that later we can match both things.    
    """

    cdef  data_container dc = <data_container> ctx

    cdef char * c_str_value
    cdef str py_str_value
    cdef int8_t c_int8_value
    cdef int16_t c_int16_value
    cdef int32_t c_int32_value
    cdef int64_t c_int64_value
    cdef float c_float_value
    cdef double c_double_value
    cdef long py_long_value
    cdef double py_float_value
    cdef py_variable_format pyformat
    cdef object labels_raw
    cdef str var_label
    cdef object cur_dir
    cdef str value_label_name

    var_label = <str> val_labels
    value_label_name = <str> label

    cdef readstat_type_t value_type
    value_type = readstat_value_type(value)
    
    labels_raw = dc.labels_raw
    cur_dict = labels_raw.get(var_label)
    if not cur_dict:
        cur_dict = dict()
    
    if readstat_value_is_tagged_missing(value):
        # SAS and Stata missing values
        missing_tag = <int> readstat_value_tag(value)
        # In SAS missing values are A to Z or _ in stata a to z
        cur_dict[chr(missing_tag)] = value_label_name
    else:

        if value_type == READSTAT_TYPE_STRING or value_type == READSTAT_TYPE_STRING_REF:
            c_str_value = readstat_string_value(value)
            py_str_value = <str> c_str_value
            pyformat = VAR_FORMAT_STRING
        elif value_type == READSTAT_TYPE_INT8:
            c_int8_value = readstat_int8_value(value)
            py_long_value = <long> c_int8_value
            pyformat = VAR_FORMAT_LONG
        elif value_type == READSTAT_TYPE_INT16:
            c_int16_value = readstat_int16_value(value)
            py_long_value = <long> c_int16_value
            pyformat = VAR_FORMAT_LONG
        elif value_type == READSTAT_TYPE_INT32:
            c_int32_value = readstat_int32_value(value)
            py_long_value = <long> c_int32_value
            pyformat = VAR_FORMAT_LONG
        elif value_type == READSTAT_TYPE_FLOAT:
            c_float_value = readstat_float_value(value)
            py_float_value = <double> c_float_value
            pyformat = VAR_FORMAT_FLOAT
        elif value_type == READSTAT_TYPE_DOUBLE:
            c_double_value = readstat_double_value(value);
            py_float_value = <double> c_double_value
            pyformat = VAR_FORMAT_FLOAT
        else:
            raise PyreadstatError("Unkown data type")


        if pyformat == VAR_FORMAT_STRING:
            cur_dict[py_str_value] = value_label_name
        elif pyformat == VAR_FORMAT_LONG:
            cur_dict[py_long_value] = value_label_name
        elif pyformat == VAR_FORMAT_FLOAT:
            cur_dict[py_float_value] = value_label_name
        elif pyformat == VAR_FORMAT_MISSING:
            pass
        else:
            raise PyreadstatError("Failed convert C to python value")

    dc.labels_raw[var_label] = cur_dict

    return READSTAT_HANDLER_OK

cdef int handle_note (int note_index, char *note, void *ctx) except READSTAT_HANDLER_ABORT:
    """
    Collects notes (text annotations) attached to the documents. It happens for spss and stata
    """

    cdef str pynote
    cdef  data_container dc = <data_container> ctx

    pynote = <str> note
    dc.notes.append(pynote)

    return READSTAT_HANDLER_OK

cdef int handle_open(const char *u8_path, void *io_ctx) except READSTAT_HANDLER_ABORT:
    """
    Special open handler for windows in order to be able to handle paths with international characters
    Courtesy of Jonathon Love.
    """
    cdef int fd
    cdef Py_ssize_t length

    if not os.path.isfile(u8_path):
        return -1

    if os.name == "nt":
        u16_path = PyUnicode_AsWideCharString(u8_path, &length)
        fd = _wsopen(u16_path, _O_RDONLY | _O_BINARY, _SH_DENYWR, 0)
        assign_fd(io_ctx, fd)
        return fd
    else:
        return -1


cdef object _file_object_ctx = None

cdef int pyobject_open_handler(const char *path, void *io_ctx) noexcept:
    """File is already open - this is a no-op"""
    return 0

cdef int pyobject_close_handler(void *io_ctx) noexcept:
    """User manages file lifetime - this is a no-op"""
    return 0

cdef ssize_t pyobject_read_handler(void *buf, size_t nbyte, void *io_ctx) noexcept:
    """Bridge Python file.read() to C read operation"""
    global _file_object_ctx
    cdef bytes data
    cdef ssize_t bytes_read
    cdef char *data_ptr
    
    try:
        file_obj = _file_object_ctx
        data = file_obj.read(nbyte)
        bytes_read = len(data)
        if bytes_read > 0:
            data_ptr = <char*>data
            memcpy(buf, data_ptr, bytes_read)
        return bytes_read
    except:
        return -1

cdef readstat_off_t pyobject_seek_handler(readstat_off_t offset, readstat_io_flags_t whence, void *io_ctx) noexcept:
    """Bridge Python file.seek() to C seek operation"""
    global _file_object_ctx
    cdef int py_whence
    
    try:
        file_obj = _file_object_ctx
        if whence == READSTAT_SEEK_SET:
            py_whence = 0
        elif whence == READSTAT_SEEK_CUR:
            py_whence = 1
        else:  # READSTAT_SEEK_END
            py_whence = 2
        
        file_obj.seek(offset, py_whence)
        return file_obj.tell()
    except:
        return -1


cdef void check_exit_status(readstat_error_t retcode) except *:
    """
    transforms a readstat exit status to a python error if status is not READSTAT OK
    """

    cdef char * err_readstat
    cdef str err_message
    if retcode != READSTAT_OK:
        err_readstat = readstat_error_message(retcode)
        err_message = <str> err_readstat
        raise ReadstatError(err_message)


cdef void run_readstat_parser(char * filename, data_container data, py_file_extension file_extension, long row_limit, long row_offset, object file_obj=None) except *:
    """
    Runs the parsing of the file by readstat library.
    
    If file_obj is provided, it will be used instead of filename for I/O operations.
    """
    global _file_object_ctx
    
    cdef readstat_parser_t *parser
    cdef readstat_error_t error
    cdef readstat_metadata_handler metadata_handler
    cdef readstat_variable_handler variable_handler
    cdef readstat_value_handler value_handler
    cdef readstat_value_label_handler value_label_handler
    cdef readstat_note_handler note_handler
    cdef readstat_open_handler open_handler
    cdef readstat_close_handler close_handler
    cdef readstat_read_handler read_handler
    cdef readstat_seek_handler seek_handler

    cdef void *ctx
    cdef str err_message
    cdef PyObject *pyerr
    cdef bint metaonly
    cdef char *err_readstat
    cdef bytes encoding_byte

    metaonly = data.metaonly
    ctx = <void *>data
    
    #readstat_error_t error = READSTAT_OK;
    parser = readstat_parser_init()
    metadata_handler = <readstat_metadata_handler> handle_metadata
    variable_handler = <readstat_variable_handler> handle_variable
    value_handler = <readstat_value_handler> handle_value
    value_label_handler = <readstat_value_label_handler> handle_value_label
    note_handler = <readstat_note_handler> handle_note
    
    
    check_exit_status(readstat_set_metadata_handler(parser, metadata_handler))
    check_exit_status(readstat_set_variable_handler(parser, variable_handler))
    check_exit_status(readstat_set_value_label_handler(parser, value_label_handler))
    check_exit_status(readstat_set_note_handler(parser, note_handler))

    # Set up custom I/O handlers for file objects
    if file_obj is not None:
        _file_object_ctx = file_obj
        open_handler = <readstat_open_handler> pyobject_open_handler
        close_handler = <readstat_close_handler> pyobject_close_handler
        read_handler = <readstat_read_handler> pyobject_read_handler
        seek_handler = <readstat_seek_handler> pyobject_seek_handler
        readstat_set_open_handler(parser, open_handler)
        readstat_set_close_handler(parser, close_handler)
        readstat_set_read_handler(parser, read_handler)
        readstat_set_seek_handler(parser, seek_handler)
    elif os.name == "nt":
        # on windows we need a custom open handler in order to deal with internation characters in the path.
        open_handler = <readstat_open_handler> handle_open
        readstat_set_open_handler(parser, open_handler)

    if not metaonly:
        check_exit_status(readstat_set_value_handler(parser, value_handler))

    # if the user set the encoding manually
    if data.user_encoding:
        encoding_bytes = data.user_encoding.encode("utf-8")
        readstat_set_file_character_encoding(parser, <char *> encoding_bytes)

    if row_limit:
        check_exit_status(readstat_set_row_limit(parser, row_limit))
    
    if row_offset:
        check_exit_status(readstat_set_row_offset(parser, row_offset))

    # parse!
    if file_extension == FILE_EXT_SAV:
        error = readstat_parse_sav(parser, filename, ctx);
    elif file_extension == FILE_EXT_SAS7BDAT:
        error = readstat_parse_sas7bdat(parser, filename, ctx);
    elif file_extension == FILE_EXT_DTA:
        error = readstat_parse_dta(parser, filename, ctx);
    elif file_extension == FILE_EXT_XPORT:
        error = readstat_parse_xport(parser, filename, ctx);
    elif file_extension == FILE_EXT_POR:
        error = readstat_parse_por(parser, filename, ctx);
    elif file_extension == FILE_EXT_SAS7BCAT:
        error = readstat_parse_sas7bcat(parser, filename, ctx);
    #error = parse_func(parser, filename, ctx);
    readstat_parser_free(parser)
    # check if a python error ocurred, if yes, it will be printed by the interpreter, 
    # if not, make sure that the return from parse_func is OK, if not print
    pyerr = PyErr_Occurred()
    if <void *>pyerr == NULL:
        check_exit_status(error)
        

cdef object data_container_to_dict(data_container data):
    """
    Transforms a data container object to a pandas data frame
    """
    
    cdef object final_container
    cdef list col_data
    cdef list col_names
    cdef str cur_name_str
    cdef int fc_cnt
    cdef bint is_unkown_number_rows
    cdef int max_n_obs
    cdef bint metaonly

    final_container = OrderedDict()
    col_data = data.col_data
    col_names = data.col_names
    is_unkown_number_rows = data.is_unkown_number_rows
    max_n_obs = data.max_n_obs
    metaonly = data.metaonly
    
    for fc_cnt in range(0, len(col_names)):
        cur_name_str = col_names[fc_cnt]
        cur_data = col_data[fc_cnt]
        if is_unkown_number_rows and not metaonly:
            cur_data = cur_data[0:max_n_obs]
        if not metaonly:
            final_container[cur_name_str] = cur_data
        else:
            final_container[cur_name_str] = list() 

    return final_container

cdef object dict_to_dataframe(object dict_data, data_container dc):
    """
    Transforms a dict of numpy arrays to a pandas data frame
    """

    cdef bint dates_as_pandas, allsame
    cdef int index
    cdef str column, output_format, col_name
    cdef py_datetime_format var_format
    cdef list dtypes, sertypes, datetime_cols, date_cols
    cdef dict schema = None

    dates_as_pandas = dc.dates_as_pandas
    output_format = dc.output_format

    if dict_data:
        #schema = None
        # in polars if missing user values we need to explicitly set the type 
        # of that column to Object if not all the elements are of the same type
        if output_format != "pandas" and dc.missing_user_values:
            schema = dict()
            for indx in range(0, len(dc.col_names)):
                col_name = dc.col_names[indx]
                if indx in dc.missing_user_values.keys():
                    sertypes = [type(x) for x in dict_data[col_name] if x is not None]
                    # all missing, let polars decide
                    if not sertypes:
                        schema[col_name] = None
                        continue
                    allsame = all([x==sertypes[0] for x in sertypes])
                    # all the same: let polars decide
                    if allsame:
                        schema[col_name] = None
                    # not all the same type, has to be Object
                    else:
                        schema[col_name] = nw.Object
                else:
                    schema[col_name] = None

        data_frame = nw.from_dict(dict_data, backend=output_format, schema=schema)
        natnamespace = nw.get_native_namespace(data_frame)
        data_frame = data_frame.to_native()

        if dates_as_pandas and output_format=="pandas":
            pd = natnamespace
            dtypes = data_frame.dtypes.tolist()
            # check that datetime columns are datetime type
            # this is needed in case all date values are nan
            for index, column in enumerate(data_frame.columns):
                var_format = dc.col_formats[index]
                if dtypes[index] != '<M8[ns]' and (var_format == DATE_FORMAT_DATE or var_format == DATE_FORMAT_DATETIME):
                    data_frame.loc[:, column] = pd.to_datetime(data_frame[column])

        if output_format == "polars" and not dc.no_datetime_conversion:
            # datetime and date vectorized conversion
            pl = natnamespace
            datetime_cols = list()
            date_cols = list()
            for index, column in enumerate(data_frame.columns):
                var_format = dc.col_formats[index]
                if var_format == DATE_FORMAT_DATETIME:
                    datetime_cols.append(column)
                if var_format == DATE_FORMAT_DATE:
                    date_cols.append(column)
            if datetime_cols:
                data_frame = data_frame.with_columns(
                    [
                        pl.from_epoch(
                            (pl.col(c) % 1 * 1e6).round().cast(pl.Int64) + (
                                    pl.col(c).floor() * 1e6).cast(pl.Int64) - (
                                    pl.lit(dc.unix_to_origin_secs) * 1e6).cast(pl.Int64),
                            time_unit='us')
                        for c in datetime_cols if data_frame[c].len() > 0
                    ]
                )
            if date_cols:
                data_frame = data_frame.with_columns(pl.from_epoch(pl.col(*date_cols), time_unit='d'))

    else:
        data_frame = nw.from_dict(dict_data, backend=output_format).to_native()

    return data_frame

cdef object data_container_extract_metadata(data_container data):
    """
    Extracts metadata from a data container and puts it into a metadata 
    object
    """
    cdef list col_names_byte, col_labels_byte
    cdef str colstr
    cdef bint metaonly
    cdef bint is_unkown_number_rows
    cdef object label_to_var_name
    cdef object labels_raw
    cdef str var_name, var_label
    cdef object current_labels
    cdef object labels_str
    cdef object original_types
    cdef readstat_type_t var_type

    metaonly = data.metaonly
    is_unkown_number_rows = data.is_unkown_number_rows
    
    cdef object metadata = metadata_container()

    # mr sets
    metadata.mr_sets = data.mr_sets

    # number of rows
    metadata.number_columns = data.n_vars
    if is_unkown_number_rows:
        if not metaonly:
            metadata.number_rows = data.max_n_obs
    else:
        metadata.number_rows = data.n_obs

    # value labels
    labels_raw = data.labels_raw
    label_to_var_name = data.label_to_var_name
    variable_value_labels = dict()

    if labels_raw:
        for var_name, var_label in label_to_var_name.items():
            current_labels = labels_raw.get(var_label)
            if current_labels:
                variable_value_labels[var_name] = current_labels

    original_types = dict()
    readstat_types = dict()
    for indx in range(metadata.number_columns):
        cur_col = data.col_names[indx]
        cur_type = data.col_formats_original[indx]
        original_types[cur_col] = cur_type
        var_type = data.col_dtypes[indx]
        if var_type == READSTAT_TYPE_STRING or var_type == READSTAT_TYPE_STRING_REF:
            readstat_types[cur_col] = "string"
        elif var_type == READSTAT_TYPE_INT8:
            readstat_types[cur_col] = "int8"
        elif var_type == READSTAT_TYPE_INT16:
            readstat_types[cur_col] = "int16"
        elif var_type == READSTAT_TYPE_INT32:
            readstat_types[cur_col] = "int32"
        elif var_type == READSTAT_TYPE_FLOAT:
            readstat_types[cur_col] = "float"
        elif var_type == READSTAT_TYPE_DOUBLE:
            readstat_types[cur_col] = "double"
        else:
            raise PyreadstatError("Unkown data type")

    for indx, curset in data.missing_user_values.items():
        cur_col = data.col_names[indx]
        metadata.missing_user_values[cur_col] = sorted(list(curset))    

    metadata.notes = data.notes
    metadata.column_names = data.col_names
    metadata.column_labels = data.col_labels
    metadata.column_names_to_labels = {k:v for k,v in zip(data.col_names, data.col_labels)}
    metadata.file_encoding = data.file_encoding
    metadata.file_label = data.file_label
    metadata.variable_value_labels = variable_value_labels
    metadata.value_labels = labels_raw
    metadata.variable_to_label = label_to_var_name
    metadata.original_variable_types = original_types
    metadata.readstat_variable_types = readstat_types
    metadata.table_name = data.table_name
    metadata.missing_ranges = data.missing_ranges
    metadata.variable_storage_width = data.variable_storage_width
    metadata.variable_display_width = data.variable_display_width
    metadata.variable_alignment = data.variable_alignment
    metadata.variable_measure = data.variable_measure
    metadata.creation_time = datetime.datetime.fromtimestamp(data.ctime)
    metadata.modification_time = datetime.datetime.fromtimestamp(data.mtime)
    
    return metadata


cdef object run_conversion(object filename_path, py_file_format file_format, py_file_extension file_extension,
                           str encoding, bint metaonly, bint dates_as_pandas, list usecols, bint usernan,
                           bint no_datetime_conversion, long row_limit, long row_offset, str output_format, 
                           list extra_datetime_formats, list extra_date_formats, list extra_time_formats):
    """
    Coordinates the activities to parse a file. This is the entry point 
    for the public methods.
    
    filename_path can be a file path (str, bytes, Path-like) or a file-like object
    with read() and seek() methods.
    """
    
    cdef bytes filename_bytes
    cdef char * filename    
    cdef data_container data
    cdef object origin
    cdef set allowed_formats
    cdef object data_dict
    cdef object data_frame
    cdef object file_obj = None

    # Check if filename_path is a file-like object
    if hasattr(filename_path, 'read') and hasattr(filename_path, 'seek'):
        file_obj = filename_path
        filename_bytes = b""  # Empty path - will be ignored when file_obj is used
    elif hasattr(os, 'fsencode'):
        try:
            filename_bytes = os.fsencode(filename_path)
        except UnicodeError:
            warnings.warn("file path could not be encoded with %s which is set as your system encoding, trying to encode it as utf-8. Please set your system encoding correctly." % sys.getfilesystemencoding())
            filename_bytes = os.fsdecode(filename_path).encode("utf-8", "surrogateescape")
    else:
        if type(filename_path) == str:
            filename_bytes = filename_path.encode('utf-8')
        elif type(filename_path) == bytes:
            filename_bytes = filename_path
        else:
            raise PyreadstatError("path must be either str or bytes")
        if type(filename_path) not in (str, bytes, unicode):
            raise PyreadstatError("path must be str, bytes or unicode")
        filename_bytes = filename_path.encode('utf-8')

    # Only check file existence for path-based reads
    if file_obj is None:
        filename_bytes = os.path.expanduser(filename_bytes)
        if not os.path.isfile(filename_bytes):
            raise PyreadstatError("File {0} does not exist!".format(filename_path))

    if output_format is None:
        output_format = 'pandas'
    allowed_formats = {'pandas', 'dict', 'polars'}
    if output_format not in allowed_formats:
        raise PyreadstatError("output format must be one of {allowed_formats}, '{output_format}' was given".format(allowed_formats=allowed_formats, output_format=output_format))
    if output_format == "pandas":
        try:
            import pandas
        except:
            raise PyreadstatError("You requested pandas as output_format but cannot import pandas")
    if output_format == "polars":
        try:
            import polars
        except:
            raise PyreadstatError("You requested polars as output_format but cannot import polars")


    if extra_date_formats is not None:
        if file_format == FILE_FORMAT_SAS:
            sas_date_formats.extend([format_regex.match(edf).group(1) for edf in extra_date_formats if format_regex.match(edf)])
        elif file_format == FILE_FORMAT_SPSS:
            spss_date_formats.extend([format_regex.match(edf).group(1) for edf in extra_date_formats if format_regex.match(edf)])
        elif file_format == FILE_FORMAT_STATA:
            stata_date_formats.extend(extra_date_formats)
        else:
            raise PyreadstatError("Unknown file format")
    if extra_datetime_formats is not None:
        if file_format == FILE_FORMAT_SAS:
            sas_datetime_formats.extend([format_regex.match(edtf).group(1) for edtf in extra_datetime_formats if format_regex.match(edtf)])
        elif file_format == FILE_FORMAT_SPSS:
            spss_datetime_formats.extend([format_regex.match(edtf).group(1) for edtf in extra_datetime_formats if format_regex.match(edtf)])
        elif file_format == FILE_FORMAT_STATA:
            stata_datetime_formats.extend(extra_datetime_formats)
        else:
            raise PyreadstatError("Unknown file format")
    if extra_time_formats is not None:
        if file_format == FILE_FORMAT_SAS:
            sas_time_formats.extend([format_regex.match(etf).group(1) for etf in extra_time_formats if format_regex.match(etf)])
        elif file_format == FILE_FORMAT_SPSS:
            spss_time_formats.extend([format_regex.match(etf).group(1) for etf in extra_time_formats if format_regex.match(etf)])
        elif file_format == FILE_FORMAT_STATA:
            stata_time_formats.extend(extra_time_formats)
        else:
            raise PyreadstatError("Unknown file format")
    global sas_all_formats, spss_all_formats, stata_all_formats
    sas_all_formats = sas_date_formats + sas_datetime_formats + sas_time_formats
    stata_all_formats = stata_date_formats + stata_datetime_formats + stata_time_formats
    spss_all_formats = spss_date_formats + spss_datetime_formats + spss_time_formats

    filename = <char *> filename_bytes
    
    data = data_container()
    ctx = <void *>data        
    
    data.file_format = file_format
    data.metaonly = metaonly
    data.dates_as_pandas = dates_as_pandas
    data.output_format = output_format

    if encoding:
        data.user_encoding = encoding
    
    if file_format == FILE_FORMAT_SAS:
        origin = sas_origin
        unix_to_origin_secs = sas_secs_from_unix 
    elif file_format == FILE_FORMAT_SPSS:
        origin = spss_origin
        unix_to_origin_secs = spss_secs_from_unix 
    elif file_format == FILE_FORMAT_STATA:
        origin = stata_origin
        unix_to_origin_secs = stata_secs_from_unix 
    else:
        raise PyreadstatError("Unknown file format")
    
    data.origin = origin
    data.unix_to_origin_secs = unix_to_origin_secs

    if usecols is not None:
        data.filter_cols = 1
        data.use_cols = usecols

    data.usernan = usernan
    data.no_datetime_conversion = no_datetime_conversion
    
    # go!
    run_readstat_parser(filename, data, file_extension, row_limit, row_offset, file_obj)    
    data_dict = data_container_to_dict(data)
    if output_format == 'dict':
        data_frame = data_dict
    else:
        data_frame = dict_to_dataframe(data_dict, data)
    metadata = data_container_extract_metadata(data)

    return data_frame, metadata
    
def parser_entry_point(filename_path, str parser_format=None,
                       metadataonly=False, dates_as_pandas_datetime=False, 
             formats_as_category=True, formats_as_ordered_category=False, str encoding=None, list usecols=None, user_missing=False,
             disable_datetime_conversion=False, int row_limit=0, int row_offset=0, str output_format=None, list extra_datetime_formats=None, 
             list extra_date_formats=None, list extra_time_formats=None):


    cdef py_file_format file_format
    cdef py_file_extension file_extension

    if parser_format == "sav/zsav":
        file_format = FILE_FORMAT_SPSS
        file_extension = FILE_EXT_SAV
    elif parser_format == "sas7bdat":
        file_format = FILE_FORMAT_SAS
        file_extension = FILE_EXT_SAS7BDAT
    elif parser_format == "xport":
        file_format = FILE_FORMAT_SAS
        file_extension = FILE_EXT_XPORT
    elif parser_format == "dta":
        file_format = FILE_FORMAT_STATA
        file_extension = FILE_EXT_DTA
    elif parser_format == "por":
        file_format = FILE_FORMAT_SPSS
        file_extension = FILE_EXT_POR
    elif parser_format == "sas7bcat":
        file_format = FILE_FORMAT_SAS
        file_extension = FILE_EXT_SAS7BCAT
    else:
        raise PyreadstatError("wrong parser format")

    cdef bint metaonly = 0
    if metadataonly:
        metaonly = 1

    cdef bint dates_as_pandas = 0
    if dates_as_pandas_datetime:
        dates_as_pandas = 1

    cdef bint usernan = 0
    if user_missing:
        usernan = 1

    cdef bint no_datetime_conversion = 0
    if disable_datetime_conversion:
        no_datetime_conversion = 1
    
    data_frame, metadata = run_conversion(filename_path, file_format, file_extension, encoding, metaonly,
                                          dates_as_pandas, usecols, usernan, no_datetime_conversion, <long>row_limit, <long>row_offset,
                                          output_format, extra_datetime_formats, extra_date_formats, extra_time_formats)

    return data_frame, metadata

