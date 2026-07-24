# cython: c_string_type=unicode, c_string_encoding=utf8, language_level=2
# ###cython: profile=True
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
import os
import warnings
import sys
import datetime
#import calendar
from datetime import timezone
#from datetime import timezone as _timezone
#from libc.math cimport round, NAN

import numpy as np
import narwhals.stable.v2 as nw

from readstat_api cimport *
from _readstat_parser import ReadstatError, PyreadstatError
from _readstat_parser cimport check_exit_status

cdef set int_types = {nw.Int32, nw.Int16, nw.Int8, nw.UInt16, nw.UInt8, }
cdef set float_types = {nw.Float64, nw.Float32, nw.Decimal, nw.Int128, nw.Int64, nw.UInt128, nw.UInt64, nw.UInt32}
cdef set nat_types = {datetime.datetime, np.datetime64, datetime.time, datetime.date} #pd._libs.tslibs.timestamps.Timestamp,
cdef set pyrwriter_datetimelike_types = {PYWRITER_DATE, PYWRITER_DATETIME, PYWRITER_TIME, PYWRITER_DATETIME64, PYWRITER_DATE64, PYWRITER_TIME64}
cdef set pywriter_numeric_types = {PYWRITER_DOUBLE, PYWRITER_INTEGER, PYWRITER_LOGICAL, PYWRITER_DATE, PYWRITER_DATETIME, PYWRITER_TIME, 
                                   PYWRITER_DATETIME64, PYWRITER_DATE64, PYWRITER_TIME64 }
cdef dict narwhals_to_readstat_types = {PYWRITER_DOUBLE: READSTAT_TYPE_DOUBLE, PYWRITER_INTEGER: READSTAT_TYPE_INT32,
                                      PYWRITER_CHARACTER: READSTAT_TYPE_STRING, PYWRITER_LOGICAL: READSTAT_TYPE_INT32,
                                      PYWRITER_DTA_STR_REF: READSTAT_TYPE_STRING_REF,
                                      PYWRITER_OBJECT: READSTAT_TYPE_STRING, PYWRITER_DATE: READSTAT_TYPE_DOUBLE,
                                      PYWRITER_DATETIME: READSTAT_TYPE_DOUBLE, PYWRITER_TIME: READSTAT_TYPE_DOUBLE,
                                      PYWRITER_DATETIME64: READSTAT_TYPE_DOUBLE, PYWRITER_DATE64: READSTAT_TYPE_DOUBLE, 
                                      PYWRITER_TIME64: READSTAT_TYPE_DOUBLE, }

cdef double spss_offset_secs = 12219379200
cdef double sas_offset_secs = 315619200
cdef double spss_offset_days = 141428
cdef double sas_offset_days = 3653
cdef object date_0 = datetime.datetime(1970,1,1).date()

cdef valid_user_missing_sas = [chr(x) for x in range(ord("A"), ord("Z")+1)] + ["_"]
cdef valid_user_missing_stata = [chr(x) for x in range(ord("a"), ord("z")+1)]

# max lenght of a string in dta before it has to use string_ref
cdef int dta_old_max_width = 128
cdef int dta_111_max_width = 244
cdef int dta_117_max_width = 2045

cdef object vectorized_convert_datetime_to_number(object df, dst_file_format file_format, list pywriter_types, list pywriter_timeunits, int col_count):
    """
    transforms datetime64 columns in the dataframe to floats
    """
    cdef dict convfacs
    cdef double offset_secs
    cdef double mulfac = 1.0
    cdef int col_indx
    cdef list col_indxs
    cdef double convfac

    if file_format == FILE_FORMAT_SAV or file_format == FILE_FORMAT_POR:
        offset_secs = spss_offset_secs
    else:
        offset_secs = sas_offset_secs
    if file_format == FILE_FORMAT_DTA:
        # stata stores in milliseconds
        mulfac = 1000.0
    convfacs = {'ns': 1e9, 'us': 1e6, 'ms': 1e3}

    col_indxs = list()
    for col_indx in range(col_count):
        if pywriter_types[col_indx] == PYWRITER_DATETIME64: 
            col_indxs.append(col_indx)

    df = df.with_columns(nw.nth(col_indxs).cast(nw.Int64))
    for col_indx in col_indxs:
        convfac = convfacs[pywriter_timeunits[col_indx]]
        df = df.with_columns(nw.when(nw.nth(col_indx)!=-9223372036854775808).then(nw.nth(col_indx)))
        df = df.with_columns((((nw.nth(col_indx).cast(nw.Float64))/convfac) + offset_secs).round() * mulfac)
    return df


cdef object vectorized_convert_date_to_number(object df, dst_file_format file_format, list pywriter_types,  int col_count):
    """
    transforms date64 columns in the dataframe to floats
    """
    cdef double offset_days
    cdef double mulfac = 1.0
    cdef int col_indx
    cdef list col_indxs

    if file_format == FILE_FORMAT_SAV or file_format == FILE_FORMAT_POR:
        offset_days = spss_offset_days
        # spss stores in seconds
        mulfac = 86400
    else:
        offset_days = sas_offset_days
    col_indxs = list()
    for col_indx in range(col_count):
        if pywriter_types[col_indx] == PYWRITER_DATE64: 
            col_indxs.append(col_indx)

    df = df.with_columns(nw.nth(col_indxs).cast(nw.Int64))
    for col_indx in col_indxs:
        df = df.with_columns(nw.when(nw.nth(col_indx)!=-9223372036854775808).then(nw.nth(col_indx)))
        df = df.with_columns((nw.nth(col_indx).cast(nw.Float64) + offset_days).round() * mulfac)
    return df

cdef object vectorized_convert_time_to_number(object df, dst_file_format file_format, list pywriter_types,  int col_count):
    """
    transforms time64 columns in the dataframe to floats
    """
    cdef double mulfac = 1.0
    cdef int col_indx
    cdef list col_indxs

    if file_format == FILE_FORMAT_DTA:
        # stata stores in milliseconds
        mulfac = 1000.0

    col_indxs = list()
    for col_indx in range(col_count):
        if pywriter_types[col_indx] == PYWRITER_TIME64: 
            col_indxs.append(col_indx)

    df = df.with_columns(nw.nth(col_indxs).cast(nw.Int64))
    for col_indx in col_indxs:
        df = df.with_columns(nw.when(nw.nth(col_indx)!=-9223372036854775808).then(nw.nth(col_indx)))
        df = df.with_columns((nw.nth(col_indx).cast(nw.Float64)/1e9).round() * mulfac)
    return df

cdef double convert_datetimelike_to_number(dst_file_format file_format, pywriter_variable_type curtype, object curval) except *:
    """
    converts a datime like python/pandas object to a float
    """

    cdef double offset_days, tstamp

    if file_format == FILE_FORMAT_SAV or file_format == FILE_FORMAT_POR:
        offset_days = spss_offset_days
        offset_secs = spss_offset_secs
    else:
        offset_days = sas_offset_days
        offset_secs = sas_offset_secs

    if curtype == PYWRITER_DATETIME or curtype == PYWRITER_DATETIME64:
        # get timestamp in seconds
        if type(curval) == datetime.datetime:
            tstamp = curval.replace(tzinfo=timezone.utc).timestamp() # works only in python 3
            #tstamp = calendar.timegm(curval.replace(tzinfo=_timezone.utc).timetuple())

        tstamp += offset_secs
        if file_format == FILE_FORMAT_DTA:
            # stata stores in milliseconds
            tstamp *= 1000

    elif curtype == PYWRITER_DATE:
        if type(curval) == datetime.date:
            days = curval - date_0
            tstamp = days.days
            tstamp += offset_days
            if file_format == FILE_FORMAT_SAV or file_format == FILE_FORMAT_POR:
                # spss stores in seconds
                tstamp *= 86400

    elif curtype == PYWRITER_TIME:
        if type(curval) == datetime.time:
            tdelta = datetime.datetime.combine(datetime.date.min, curval) - datetime.datetime.min
            tstamp = tdelta.total_seconds()
            #tstamp += offset * 86400
            if file_format == FILE_FORMAT_DTA:
                # stata stores in milliseconds
                tstamp *= 1000

    return tstamp

cdef char * get_datetimelike_format_for_readstat(dst_file_format file_format, pywriter_variable_type curtype):
    """
    gives back a string with the format of the variable (according to the final application) to be used by readstat 
    """

    if curtype == PYWRITER_DATE or curtype == PYWRITER_DATE64:
        if file_format == FILE_FORMAT_DTA:
            return "%td"
        #elif file_format == FILE_FORMAT_SAV:
        #    return "DATE11"
        else:
            return "DATE"
    elif curtype == PYWRITER_DATETIME or curtype == PYWRITER_DATETIME64:
        if file_format == FILE_FORMAT_DTA:
            return "%tc"
        #elif file_format == FILE_FORMAT_SAV:
        #    return "DATETIME20"
        else:
            return "DATETIME"
    elif curtype == PYWRITER_TIME or curtype == PYWRITER_TIME64:
        if file_format == FILE_FORMAT_DTA:
            return "%tcHH:MM:SS"
        else:
            return "TIME"
    else:
        raise PyreadstatError("Unknown pywriter variable format")

cdef int get_narwhals_str_series_max_length(object series, dict value_labels, bint isobject):
    """ For a string series get the max length of the strings. Assumes there is no NaN among the elements. 
    """
    cdef object val
    cdef bytes temp
    cdef int max_length = 1
    cdef int curlen
    cdef list labels
    for val in series:
        if isobject:
            val = str(val)
        temp = val.encode("utf-8")
        curlen = len(temp)
        if curlen > max_length:
            max_length = curlen
    if value_labels:
        labels = list(value_labels.keys())
        for lab in labels:
            curlen = len(str(lab))
            if curlen > max_length:
                max_length = curlen

    return max_length


cdef int check_series_all_same_types(object series, object type_to_check):
    """
    1 if all elements in a series are of type type_to_check, 0 otherwise
    """
    for val in series:
        if type(val) != type_to_check:
            return 0
    return 1

cdef list get_narwhals_column_types(object df, dict missing_user_values, dict variable_value_labels, int dta_str_max_len):
    """
    From a narwhals data frame, get a list with tuples column types as first element, max_length as second, is_missing
    as third and time units as fourth.
    max_lenght is the max length of a string or string representation of an object, 0 for numeric types. is_missing flags
    wether the series has missing values (1) or not (0)
    dta_str_max_len is the max length for a dta string, 0 if the file format is not dta
    """

    cdef int max_length, isobject
    cdef bint has_missing
    cdef list result = list()
    cdef int equal, is_missing

    if variable_value_labels is None:
        variable_value_labels = dict()

    for curseries in df.iter_columns():
        col_name = curseries.name
        col_type = curseries.dtype

        # if categorical, let's use the type of the categories
        if col_type == nw.Categorical:
            col_type = curseries.cat.get_categories().dtype

        # we need to remove nulls for series inspection
        has_missing = <bint>curseries.null_count()>0
        if has_missing:
            curseries = curseries.drop_nulls()
        
        # missing_user_values are strings that represent a missing value in a numeric column
        # they are object because it is a mix of string and number.
        # it could also be that the type of the series is string because it contains only
        # missing values and missing_user_values.
        # we need to take out those before inspecting the object series futher
        curuser_missing = None
        if missing_user_values:
            curuser_missing = missing_user_values.get(col_name)
        if curuser_missing:
            if not df.implementation.is_pandas() and col_type == nw.Object:
                curseries = [x for x in curseries if x not in curuser_missing]
            else:
                curseries = curseries.filter(~curseries.is_in(curuser_missing))
            if not len(curseries):
                result.append((PYWRITER_DOUBLE, 0, 1, None))
                continue
        
        max_length = 0
        curtype = None
        equal = True
        # let's deal first with object type, Enum could also be anything
        if col_type == nw.Object or col_type==nw.Enum:
            curtype = type(curseries[0])
            equal = check_series_all_same_types(curseries, curtype)
            # if all elements are equal, they could be a few we expect to be an object class
            # or it could be that they are some other common types (numeric) after removing the missing_user_values
            # therefore if one of these conditions are not met they continuing flowing into the next if
            # if not all the elements of the series are from the same type, then we deal with it at the end
            # as other types we don't know what to do with
            if equal:
                # types expected to be object
                if curtype == datetime.date:
                    result.append((PYWRITER_DATE, 0, has_missing, None))
                    continue
                elif curtype == datetime.datetime:
                    result.append((PYWRITER_DATETIME, 0, has_missing, None))
                    continue
                elif curtype == datetime.time:
                    result.append((PYWRITER_TIME, 0, has_missing, None))
                    continue

        # numeric types: they could contain missing_user_values
        if equal and (col_type in float_types or curtype == float):
            result.append((PYWRITER_DOUBLE, 0, has_missing, None))
            continue
        elif equal and (col_type in int_types or curtype == int):
            result.append((PYWRITER_INTEGER, 0,has_missing, None))
            continue
        elif equal and (col_type == nw.Boolean or curtype == bool):
            result.append((PYWRITER_LOGICAL, 0,has_missing, None))
            continue
        # these types here should not contain missing_user_values,
        # for string we still check, as later we will raise an error
        elif col_type == nw.String or curtype == str:
            isobject = 0
            if not equal:
                isobject = 1
            max_length = get_narwhals_str_series_max_length(curseries, variable_value_labels.get(col_name), isobject)
            if dta_str_max_len and max_length >= dta_str_max_len:
                result.append((PYWRITER_DTA_STR_REF, max_length, has_missing, None))
                continue
            else:
                if isobject:
                    result.append((PYWRITER_OBJECT, max_length, has_missing, None))
                else:
                    result.append((PYWRITER_CHARACTER, max_length, has_missing, None))
                continue
        elif col_type == nw.Datetime:
            if col_type.time_unit == 'us':
                result.append((PYWRITER_DATETIME64, 0,has_missing, 'us'))
                continue
            elif col_type.time_unit == 'ns':
                result.append((PYWRITER_DATETIME64, 0,has_missing, 'ns'))
                continue
            elif col_type.time_unit == 'ms':
                result.append((PYWRITER_DATETIME64, 0,has_missing, 'ms'))
                continue
            else:
                result.append((PYWRITER_DATETIME, 0, has_missing, None))
                continue
        elif col_type == nw.Date:
            result.append((PYWRITER_DATE64, 0, has_missing, None))
            continue
        elif col_type == nw.Time:
            result.append((PYWRITER_TIME64, 0, has_missing, None))
            continue

        # if the object was not captured by any of the previous cases, we transform it to string
        max_length = get_narwhals_str_series_max_length(curseries, variable_value_labels.get(col_name), 1)
        if dta_str_max_len and max_length >= dta_str_max_len:
            result.append((PYWRITER_DTA_STR_REF, max_length, has_missing, None))
        else:
            result.append((PYWRITER_OBJECT, max_length, has_missing, None))
        continue

    return result

cdef readstat_label_set_t *set_value_label(readstat_writer_t *writer, dict value_labels, str labelset_name,
                        pywriter_variable_type curpytype, dst_file_format file_format, str variable_name, 
                        list user_missing_tags) except *:
    """
    Sets value labels for normal values and also tagged missing values (user defined missing for stata and sas)
    """

    cdef readstat_label_set_t *label_set
    cdef readstat_type_t curtype
    cdef double double_val
    
    curtype = narwhals_to_readstat_types[curpytype]
    label_set = readstat_add_label_set(writer, curtype, labelset_name.encode("utf-8"))

    for value, label in value_labels.items():

        if type(label) != str:
            msg = "variable_value_labels: type of Label %s in variable %s must be string" % (str(label), variable_name)
            raise PyreadstatError(msg)

        if user_missing_tags and value in user_missing_tags:
            if curpytype == PYWRITER_CHARACTER or curpytype == PYWRITER_OBJECT or curpytype==PYWRITER_DTA_STR_REF:
                msg = "missing_user_values not allowed for character variable %s" % variable_name
                raise PyreadstatError(msg)
            
            readstat_label_tagged_value(label_set, ord(value), label.encode("utf-8"))
            continue


        if curpytype == PYWRITER_DOUBLE:
            if type(value) != float and type(value) != int:
                msg = "variable_value_labels: type of Value %s in variable %s must be numeric" % (str(value), variable_name)
                raise PyreadstatError(msg)
            readstat_label_double_value(label_set, value, label.encode("utf-8"))

        elif curpytype == PYWRITER_INTEGER:
            if type(value) != int:
                #if type(value) not in int_types:
                msg = "variable_value_labels: type of Value %s in variable %s must be int" % (str(value), variable_name)
                raise PyreadstatError(msg)
            readstat_label_int32_value(label_set, value, label.encode("utf-8"))

        elif curpytype == PYWRITER_LOGICAL:
            if type(value) != bool and (value != 0 and value != 1):
                msg = "variable_value_labels: type of Value %s in variable %s must be boolean or be 1 or 0" % (str(value), variable_name)
                raise PyreadstatError(msg)
            readstat_label_int32_value(label_set, int(value), label.encode("utf-8"))

        elif curpytype == PYWRITER_CHARACTER or curpytype == PYWRITER_OBJECT or curpytype==PYWRITER_DTA_STR_REF:
            value = str(value)
            readstat_label_string_value(label_set, value.encode("utf-8"), label.encode("utf-8"))

        elif curpytype in (PYWRITER_DATE, PYWRITER_DATETIME, PYWRITER_TIME, PYWRITER_DATETIME64, PYWRITER_DATE64, PYWRITER_TIME64):
            if type(value) not in nat_types:
                msg = "variable_value_labels: type of Value %s in variable %s must match the type of the column in dataframe and be of type date, datetime or time" % (str(value), variable_name)
                raise PyreadstatError(msg)
            double_val = convert_datetimelike_to_number(file_format, curpytype, value) 
            readstat_label_double_value(label_set, double_val, label.encode("utf-8"))

    return label_set

cdef void add_missing_ranges(list cur_ranges, readstat_variable_t *variable, pywriter_variable_type vartype, str variablename) except *:
    """
    Adding missing ranges, user defined missing discrete values both numeric and character,
     this happens for SPSS
    """

    cdef int range_values = 0
    cdef int discrete_values = 0
    cdef int discrete_strings = 0

    for cur_range in cur_ranges:
        if isinstance(cur_range, dict):
            hi = cur_range.get("hi")
            lo = cur_range.get("lo")
            if hi is None or lo is None:
                msg = "dictionaries in missing_ranges must have the keys hi and lo"
                raise PyreadstatError(msg)
            if type(hi) in (int, float)  and type(lo) in (int, float):
                if vartype not in pywriter_numeric_types:
                    msg = "numeric missing_ranges value given for non numeric variable %s" %variablename
                    raise PyreadstatError(msg)
                if hi == lo:
                    check_exit_status(readstat_variable_add_missing_double_value(variable, hi))
                    discrete_values += 1
                else:
                    check_exit_status(readstat_variable_add_missing_double_range(variable, lo, hi))
                    range_values += 1
            elif type(hi) == str and type(lo) == str:
                if vartype != PYWRITER_CHARACTER and vartype != PYWRITER_OBJECT and vartype !=PYWRITER_DTA_STR_REF:
                    msg = "character missing_ranges value given for non character variable %s" %variablename
                    raise PyreadstatError(msg)
                if hi == lo:
                    if len(hi) > 8:
                        msg = "missing_ranges: string values length must not be larger than 8"
                        raise PyreadstatError(msg)
                    check_exit_status(readstat_variable_add_missing_string_value(variable, hi))#.encode("utf-8")))
                    discrete_strings += 1
                else:
                    #check_exit_status(readstat_variable_add_missing_string_range(variable, lo, hi))
                    msg = "missing_ranges: hi and lo values must be both the same for string type"
                    raise PyreadstatError(msg)
            else:
                msg = "missing_ranges: hi and lo values must be both either of numeric or string type"
                raise PyreadstatError(msg)
        else:
            if type(cur_range) in (int, float):
                if vartype not in pywriter_numeric_types:
                    msg = "numeric missing_ranges value given for non numeric variable %s" %variablename
                    raise PyreadstatError(msg)
                check_exit_status(readstat_variable_add_missing_double_value(variable, cur_range))
                discrete_values += 1
            elif type(cur_range) == str:
                if vartype != PYWRITER_CHARACTER and vartype != PYWRITER_OBJECT or vartype==PYWRITER_DTA_STR_REF:
                    msg = "character missing_ranges value given for non character variable %s" %variablename
                    raise PyreadstatError(msg)
                if len(cur_range) > 8:
                        msg = "missing_ranges: string values length must not be larger than 8"
                        raise PyreadstatError(msg)
                check_exit_status(readstat_variable_add_missing_string_value(variable, cur_range))#.encode("utf-8")))
                discrete_strings += 1
            else:
                msg = "missing_ranges: values must be both either of numeric or string type"
                raise PyreadstatError(msg)
               
        if discrete_strings > 3:
            msg = "missing_ranges: max 3 string values per variable allowed"
            raise PyreadstatError(msg)
        if range_values:
            if range_values > 1:
                msg = "missing_ranges: max 1 range value per variable allowed"
                raise PyreadstatError(msg)
            if discrete_values > 1:
                msg = "missing_ranges: max 1 discrete numeric value if combined with 1 range value per variable allowed"
                raise PyreadstatError(msg)
        if discrete_values >3:
            msg = "missing_ranges: max 3 discrete numeric values per variable allowed"
            raise PyreadstatError(msg)

cdef void set_variable_alignment(readstat_variable_t *variable, str alignment_str, str var_name) except *:
    """
    Sets the variable alignment, ineffective on SPSS, STATA and XPORT (what about SAS7bdat?)
    """

    cdef readstat_alignment_t alignment

    if alignment_str == "right":
        alignment = READSTAT_ALIGNMENT_RIGHT
    elif alignment_str == "left":
        alignment = READSTAT_ALIGNMENT_LEFT
    elif alignment_str == "center":
        alignment = READSTAT_ALIGNMENT_CENTER
    elif alignment_str == "unknown":
        alignment = READSTAT_ALIGNMENT_UNKNOWN
    else:
        msg = "alignment for variable %s must be either right, center, left or unknown got %s instead" % (var_name, alignment_str)
        raise PyreadstatError(msg)

    readstat_variable_set_alignment(variable, alignment)

cdef void set_variable_display_width(readstat_variable_t *variable, int display_width, str var_name) except *:
    """
    Sets the variable display width (SPSS). Not effective on STATA. (what about SAS7BDAT?)
    """

    readstat_variable_set_display_width(variable, display_width)

cdef void set_variable_measure(readstat_variable_t *variable, str measure_str, str var_name) except *:
    """
    sets the variable measure type (SPSS). Not effective on STATA.
    """

    cdef readstat_measure_t measure

    if measure_str == "nominal":
        measure = READSTAT_MEASURE_NOMINAL
    elif measure_str == "ordinal":
        measure = READSTAT_MEASURE_ORDINAL
    elif measure_str == "scale":
        measure = READSTAT_MEASURE_SCALE
    elif measure_str == "unknown":
        measure = READSTAT_MEASURE_UNKNOWN
    else:
        msg = "measure for variable %s must be either nominal, ordinal, scale or unknown got %s instead" % (var_name, measure_str)
        raise PyreadstatError(msg)

    readstat_variable_set_measure(variable, measure);


cdef ssize_t write_bytes(const void *data, size_t _len, void *ctx) noexcept:
    """
    for the writer an explicit function to write must be defined 
    """
    cdef int fd
    fd = (<int *>ctx)[0]
    if os.name=='nt':
        return _write(fd, data, _len)
    else:
        return write(fd, data, _len)

cdef void _check_exit_status(readstat_error_t retcode) except *:
    """
    transforms a readstat exit status to a python error if status is not READSTAT OK
    """

    cdef char * err_readstat
    cdef str err_message
    if retcode != READSTAT_OK:
        err_readstat = readstat_error_message(retcode)
        err_message = <str> err_readstat
        raise ReadstatError(err_message)

cdef int open_file(bytes filename_bytes):

    cdef int fd
    cdef int flags
    cdef Py_ssize_t length

    cdef char *path

    if os.name == "nt":
        filename_str = os.fsdecode(filename_bytes)
        u16_path = PyUnicode_AsWideCharString(filename_str, &length)
        flags = _O_WRONLY | _O_CREAT | _O_BINARY | _O_TRUNC
        fd = _wsopen(u16_path, flags, _SH_DENYRW, _S_IREAD | _S_IWRITE)
    else:
        path = <char *> filename_bytes
        flags = O_WRONLY | O_CREAT | O_TRUNC
        fd = open(path, flags, 0644)

    return fd

cdef int close_file(int fd):
    if os.name == "nt":
        return _close(fd)
    else:
        return close(fd)

cdef void initial_checks(bint is_pandas, bint is_polars, dict variable_value_labels, dict missing_user_values,
                        dst_file_format file_format, list col_names, bytes filename_bytes) except *:
    """
    Running some checks before starting writing
    """

    if not is_pandas and not is_polars:
        msg = "dataframe must be pandas or polars dataframe"
        raise PyreadstatError(msg)

    if variable_value_labels:
        for k,v in variable_value_labels.items():
            if type(v) != dict:
                msg = "variable_value_labels: value for key %s must be dict, got %s" % (k, str(type(v)))
                raise PyreadstatError(msg)

    cdef list valid_user_missing
    if missing_user_values:
        if file_format == FILE_FORMAT_DTA:
            valid_user_missing = valid_user_missing_stata
        elif file_format == FILE_FORMAT_SAS7BDAT or file_format == FILE_FORMAT_SAS7BCAT:
            valid_user_missing = valid_user_missing_sas

        for key, missing_values in missing_user_values.items():
            if not isinstance(missing_values, list):
                msg = "missing_user_values: values in dictionary must be list"
                raise PyreadstatError(msg)
            for val in missing_values:
                if val not in valid_user_missing:
                    msg = "missing_user_values supports values a to z for Stata and A to Z and _ for SAS, got %s instead" % str(val)
                    raise PyreadstatError(msg)

    if len(col_names) != len(set(col_names)):
        msg = "Non unique column names detected in the dataframe!"
        raise PyreadstatError(msg)

    for variable_name in col_names:
        if type(variable_name) != str:
                raise PyreadstatError("variable name '%s' is of type %s and it must be str (not starting with numbers!)" % (variable_name, str(type(variable_name))))
        if len(variable_name) == 0:
            raise PyreadstatError("variable names must be non-empty strings, not starting with numbers")
        if not variable_name[0].isalpha():
            raise PyreadstatError("variable name '%s' starts with an illegal (non-alphabetic) character: '%s' (ordinal %s)" % (variable_name, variable_name[0], ord(variable_name[0])))
        if " " in variable_name:
            raise PyreadstatError("variable name '%s' contains a space, which is not allowed" % variable_name)

    dirname = os.path.dirname(filename_bytes)
    if dirname and  not os.path.isdir(dirname):
        raise PyreadstatError(f"the destination folder {dirname} does not exist!")

cdef bytes filepath_to_bytes(object filename_path):
    """
    transforms an object with the path to the filename to bytes
    """
    if hasattr(os, 'fsencode'):
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
    return filename_bytes


cdef int run_write(df, object filename_path, dst_file_format file_format, str file_label, object column_labels,
                   int file_format_version, object note, str table_name, dict variable_value_labels, 
                   dict missing_ranges, dict missing_user_values, dict variable_alignment,
                   dict variable_display_width, dict variable_measure, dict variable_format, bint row_compression) except *:
    """
    main entry point for writing all formats. Some parameters are specific for certain file type
    and are even incompatible between them. This function relies on the caller to select the right
    combination of parameters, not checking them otherwise.
    """

    cdef object natnamespace
    cdef object pd = None
    cdef bint is_pandas, is_polars
    cdef bytes filename_bytes
    cdef list col_names

    filename_bytes = filepath_to_bytes(filename_path)
    filename_bytes = os.path.expanduser(filename_bytes)

    df = nw.from_native(df, eager_only=True)
    natnamespace = nw.get_native_namespace(df)
    is_pandas = df.implementation.is_pandas()
    is_polars = df.implementation.is_polars()
    col_names = df.columns
    if is_pandas:
        pd = natnamespace

    initial_checks(is_pandas, is_polars, variable_value_labels, missing_user_values, file_format,
                         col_names, filename_bytes) 

    cdef readstat_error_t retcode
    cdef char *err_readstat
    cdef str err_message

    cdef readstat_writer_t *writer

    cdef bytes file_label_bytes
    cdef char *file_labl
    cdef int dta_str_max_len = 0


    if file_format == FILE_FORMAT_POR:
        col_names = [x.upper() for x in col_names]

    if file_format == FILE_FORMAT_DTA:
        if file_format_version >= 117:
            dta_str_max_len = dta_117_max_width
        elif file_format_version >= 111:
            dta_str_max_len = dta_111_max_width
        else:
            dta_str_max_len = dta_old_max_width

    cdef list col_types = get_narwhals_column_types(df, missing_user_values, variable_value_labels, dta_str_max_len)
    cdef int row_count = len(df)
    cdef int col_count = len(col_names)
    cdef dict col_names_to_types = {k:v[0] for k,v in zip(col_names, col_types)}

    cdef readstat_variable_t *variable
    cdef pywriter_variable_type curtype
    cdef int max_length
    cdef char *curformat
    cdef str tempformat
    cdef int col_indx
    cdef bytes cur_col_label
    cdef int col_label_count = 0

    cdef readstat_variable_t *tempvar
    cdef int row_indx
    cdef str curvalstr
    cdef double dtimelikeval
    #cdef np.ndarray values
    cdef object values
    cdef dict value_labels
    cdef int lblset_cnt = 0
    cdef readstat_label_set_t *label_set
    cdef list col_label_temp 
    cdef bint hasdatetime64
    cdef list pywriter_types, pywriter_timeunits
    cdef object df2
    cdef float mulfac, conv2secs
    cdef readstat_string_ref_t* strref
    cdef dict strref_map = dict()
    cdef int strref_cnt 
    cdef object strref_indx


    cdef int fd = open_file(filename_bytes)
    writer = readstat_writer_init()

    try:

        check_exit_status(readstat_set_data_writer(writer, write_bytes))

        if file_label:
            file_label_bytes = file_label.encode("utf-8")
            file_labl = <char *> file_label_bytes
            check_exit_status(readstat_writer_set_file_label(writer, file_labl))

        if note:
            if type(note) == str:
                note = [note]
            if type(note) == list:
                for line in note:
                    readstat_add_note(writer, line.encode("utf-8"))
            else:
                raise PyreadstatError(f"note should be either str or list, got {type(note)}")

        if file_format_version > -1:
            check_exit_status(readstat_writer_set_file_format_version(writer, file_format_version))

        if row_compression:
            check_exit_status(readstat_writer_set_compression(writer, READSTAT_COMPRESS_ROWS))

        # table name is used only for xpt files
        if table_name:
            table_name_bytes = table_name.encode("utf-8")
            tab_name = <char *> table_name_bytes
            check_exit_status(readstat_writer_set_table_name(writer, tab_name))

        # add variables
        if column_labels:
            if type(column_labels) != list and type(column_labels) != dict:
                raise PyreadstatError("column_labels must be either list or dict!")
            if type(column_labels) == dict:
                col_label_temp = list()
                for col_indx in range(col_count):
                    variable_name = col_names[col_indx]
                    if variable_name in column_labels.keys():
                        col_label_temp.append(column_labels[variable_name])
                    else:
                        col_label_temp.append(None)
                column_labels = col_label_temp

            col_label_count = len(column_labels)
            if col_label_count != col_count:
                raise PyreadstatError("length of column labels must be the same as number of columns")
     
        strref_cnt = 0
        for col_indx in range(col_count):
            curtype, max_length, _,_ = col_types[col_indx]
            variable_name = col_names[col_indx]
            # add variable
            variable = readstat_add_variable(writer, variable_name.encode("utf-8"), narwhals_to_readstat_types[curtype], max_length)
            # add format
            if variable_format:
                tempformat = variable_format.get(variable_name)
                if tempformat:
                   readstat_variable_set_format(variable, tempformat.encode("utf-8")) 
            if curtype in pyrwriter_datetimelike_types and (variable_format is None or variable_name not in variable_format.keys()):
                curformat = get_datetimelike_format_for_readstat(file_format, curtype)
                readstat_variable_set_format(variable, curformat)
            # prepare string_ref
            # for STRING_REF we have to add to a dict here before start writing
            if curtype == PYWRITER_DTA_STR_REF:
                for curval in df[variable_name]:
                    if curval not in strref_map:
                        curvalstr = str(curval)
                        strref = readstat_add_string_ref(writer, curvalstr.encode("utf-8"))
                        strref_map[curvalstr] = strref_cnt
                        strref_cnt += 1
            # labels
            if col_label_count:
                if column_labels[col_indx] is not None:
                    if type(column_labels[col_indx]) != str:
                        raise PyreadstatError("Column labels must be strings")
                    cur_col_label = column_labels[col_indx].encode("utf-8")
                    readstat_variable_set_label(variable, cur_col_label)
            if variable_value_labels:
                value_labels = variable_value_labels.get(variable_name)
                if value_labels:
                    labelset_name = variable_name + str(lblset_cnt)
                    lblset_cnt += 1
                    curuser_missing = None
                    if missing_user_values:
                        curuser_missing = missing_user_values.get(variable_name)
                    label_set = set_value_label(writer, value_labels, labelset_name,
                        col_names_to_types[variable_name], file_format, variable_name, curuser_missing)
                    readstat_variable_set_label_set(variable, label_set)
            # missing ranges
            if missing_ranges:
                cur_ranges = missing_ranges.get(variable_name)
                if cur_ranges:
                    if not isinstance(cur_ranges, list):
                        msg = "missing_ranges: values in dictionary must be list"
                        raise PyreadstatError(msg)
                    add_missing_ranges(cur_ranges, variable, curtype, variable_name)
            if variable_alignment:
                # At the moment this is ineffective for sav and dta (the function runs but in
                # the resulting file all alignments are still unknown)
                cur_alignment = variable_alignment.get(variable_name)
                if cur_alignment:
                    set_variable_alignment(variable, cur_alignment, variable_name)
            if variable_display_width:
                cur_display_width = variable_display_width.get(variable_name)
                if cur_display_width:
                    set_variable_display_width(variable, cur_display_width, variable_name)
            if variable_measure:
                cur_measure = variable_measure.get(variable_name)
                if cur_measure:
                    set_variable_measure(variable, cur_measure, variable_name)

        # start writing
        if file_format == FILE_FORMAT_SAS7BCAT:
            check_exit_status(readstat_begin_writing_sas7bcat(writer, &fd))
        elif file_format == FILE_FORMAT_DTA:
            check_exit_status(readstat_begin_writing_dta(writer, &fd, row_count))
        elif file_format == FILE_FORMAT_SAV:
            check_exit_status(readstat_begin_writing_sav(writer, &fd, row_count))
        elif file_format == FILE_FORMAT_POR:
            check_exit_status(readstat_begin_writing_por(writer, &fd, row_count))
        elif file_format == FILE_FORMAT_SAS7BDAT:
            check_exit_status(readstat_begin_writing_sas7bdat(writer, &fd, row_count))
        elif file_format == FILE_FORMAT_XPORT:
            check_exit_status(readstat_begin_writing_xport(writer, &fd, row_count))
        else:
            raise PyreadstatError("unknown file format")

        # validation
        check_exit_status(readstat_validate_metadata(writer))
        for col_indx in range(col_count):
            tempvar = readstat_get_variable(writer, col_indx)
            check_exit_status(readstat_validate_variable(writer, tempvar))

        # vectorized transform of datetime64ns columns
        pywriter_types = [x[0] for x in col_types]
        pywriter_timeunits = [x[3] for x in col_types]
        hasdatetime64 = PYWRITER_DATETIME64 in pywriter_types 
        hasdate64 = PYWRITER_DATE64 in pywriter_types 
        hastime64 = PYWRITER_TIME64 in pywriter_types 
        if hasdatetime64 or hasdate64 or hastime64:
            df2 = df.clone()
            if hasdatetime64:
                df2 = vectorized_convert_datetime_to_number(df2, file_format, pywriter_types, pywriter_timeunits, col_count)
            if hasdate64:
                df2 = vectorized_convert_date_to_number(df2, file_format, pywriter_types, col_count)
            if hastime64:
                df2 = vectorized_convert_time_to_number(df2, file_format, pywriter_types, col_count)
        else:
            df2 = df


        # inserting
        rowcnt = 0

        for row in df2.iter_rows():
            check_exit_status(readstat_begin_row(writer))

            for col_indx in range(col_count):

                tempvar = readstat_get_variable(writer, col_indx)
                curval = row[col_indx]
                curtype = pywriter_types[col_indx]
                is_missing = col_types[col_indx][2]
                curuser_missing = None
                if missing_user_values:
                    curuser_missing = missing_user_values.get(col_names[col_indx])
                
                if is_missing:
                    # For pandas we need isna because values can be Nan, NA, NAT, None
                    if is_pandas:
                        check_if_missing = pd.isna(curval)
                    # for other libraries the value would be None
                    else:
                        check_if_missing = curval is None
                    if check_if_missing:
                        check_exit_status(readstat_insert_missing_value(writer, tempvar))
                        continue

                if curuser_missing and curtype in pywriter_numeric_types:
                    if curval in curuser_missing:
                        check_exit_status(readstat_insert_tagged_missing_value(writer, tempvar, ord(curval)))
                        continue

                if curtype == PYWRITER_DOUBLE:
                    check_exit_status(readstat_insert_double_value(writer, tempvar, <double>curval))
                elif curtype == PYWRITER_INTEGER:
                    check_exit_status(readstat_insert_int32_value(writer, tempvar, curval))
                elif curtype == PYWRITER_LOGICAL:
                    check_exit_status(readstat_insert_int32_value(writer, tempvar, <int>curval))
                elif curtype == PYWRITER_CHARACTER:
                    check_exit_status(readstat_insert_string_value(writer, tempvar, curval.encode("utf-8")))
                elif curtype == PYWRITER_OBJECT:
                    curvalstr = str(curval)
                    check_exit_status(readstat_insert_string_value(writer, tempvar, curvalstr.encode("utf-8")))
                elif curtype == PYWRITER_DTA_STR_REF:
                    curvalstr = str(curval)
                    strref_indx = strref_map[curvalstr]
                    strref = readstat_get_string_ref(writer, strref_indx)
                    check_exit_status(readstat_insert_string_ref(writer, tempvar, strref))
                elif curtype == PYWRITER_DATETIME64 or curtype == PYWRITER_DATE64 or curtype == PYWRITER_TIME64:
                    check_exit_status(readstat_insert_double_value(writer, tempvar, <double>curval))
                elif curtype in pyrwriter_datetimelike_types:
                    dtimelikeval = convert_datetimelike_to_number(file_format, curtype, curval)
                    check_exit_status(readstat_insert_double_value(writer, tempvar, dtimelikeval))
                else:
                    raise PyreadstatError("Unknown data format to insert")

            check_exit_status(readstat_end_row(writer))
            rowcnt += 1

        check_exit_status(readstat_end_writing(writer))

    except:
        raise
    finally:
        readstat_writer_free(writer)
        close_file(fd)

    return 0

def writer_entry_point(df, dst_path, str writer_format=None, str file_label="",
                       int version=0,
                       str table_name=None,
                       object column_labels=None, compress=False, row_compress=False, object note=None,
                dict variable_value_labels=None, dict missing_ranges=None, dict variable_display_width=None,
                dict variable_measure=None,
                dict missing_user_values=None,
                dict variable_format=None,
                dict variable_alignment = None,
                       ):


    cdef int file_format_version = 0
    cdef bint row_compression = 0

    cdef dst_file_format writer_file_format

    if writer_format == "sav":
        writer_file_format = FILE_FORMAT_SAV
        file_format_version = 2
        if compress and row_compress:
            raise PyreadstatError("compress and row_compress cannot be both True")
        if compress:
            file_format_version = 3
        if row_compress:
            row_compression = 1
    elif writer_format == "dta":
        writer_file_format = FILE_FORMAT_DTA
        if version == 15:
            file_format_version = 119
        elif version == 14:
            file_format_version = 118
        elif version == 13:
            file_format_version = 117
        elif version == 12:
            file_format_version = 115
        elif version in {10, 11}:
            file_format_version = 114
        elif version in {8, 9}:
            file_format_version = 113
        else:
            raise PyreadstatError("Version not supported")
        note = ""
    elif writer_format == "xport":
        writer_file_format = FILE_FORMAT_XPORT
        file_format_version = version
    elif writer_format == "por":
        writer_file_format = FILE_FORMAT_POR
    else:
        raise PyreadstatError("wrong writer format")

    run_write(df, dst_path, writer_file_format, file_label, column_labels, 
        file_format_version, note, table_name, variable_value_labels, missing_ranges, missing_user_values,
        variable_alignment, variable_display_width, variable_measure, variable_format, row_compression)
