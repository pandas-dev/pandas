from cpython.datetime cimport (
    datetime,
    tzinfo,
)
from numpy cimport (
    int32_t,
    int64_t,
    ndarray,
)

from pandas._libs.tslibs.np_datetime cimport (
    NPY_DATETIMEUNIT,
    npy_datetimestruct,
)
from pandas._libs.tslibs.timestamps cimport _Timestamp
from pandas._libs.tslibs.timezones cimport tz_compare


cdef class _TSObject:
    cdef readonly:
        npy_datetimestruct dts      # npy_datetimestruct
        int64_t value               # numpy dt64
        tzinfo tzinfo
        bint fold
        NPY_DATETIMEUNIT creso

    cdef int64_t ensure_reso(
        self, NPY_DATETIMEUNIT creso, val=*, bint round_ok=*
    ) except? -1


cdef _TSObject convert_to_tsobject(object ts, tzinfo tz, str unit,
                                   bint dayfirst, bint yearfirst,
                                   int32_t nanos=*)

cdef _TSObject convert_datetime_to_tsobject(datetime ts, tzinfo tz,
                                            int32_t nanos=*,
                                            NPY_DATETIMEUNIT reso=*)

cdef _TSObject convert_str_to_tsobject(str ts, tzinfo tz,
                                       bint dayfirst=*,
                                       bint yearfirst=*)

cdef int64_t get_datetime64_nanos(object val, NPY_DATETIMEUNIT reso) except? -1

cpdef datetime localize_pydatetime(datetime dt, tzinfo tz)
cdef int64_t cast_from_unit(object ts, str unit, NPY_DATETIMEUNIT out_reso=*) except? -1

cdef maybe_localize_tso(_TSObject obj, tzinfo tz, NPY_DATETIMEUNIT reso)


cdef int64_t parse_pydatetime(
    datetime val,
    npy_datetimestruct *dts,
    NPY_DATETIMEUNIT creso,
) except? -1
