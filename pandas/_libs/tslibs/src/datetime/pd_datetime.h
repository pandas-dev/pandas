/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

*/

#ifndef PANDAS__LIBS_TSLIBS_SRC_DATETIME_PD_DATETIME_H_
#define PANDAS__LIBS_TSLIBS_SRC_DATETIME_PD_DATETIME_H_

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif  // NPY_NO_DEPRECATED_API

#include <numpy/ndarraytypes.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  npy_int64 days;
  npy_int32 hrs, min, sec, ms, us, ns, seconds, microseconds, nanoseconds;
} pandas_timedeltastruct;

const npy_datetimestruct _AS_MIN_DTS = {1969, 12,     31,     23,    59,
                                        50,   776627, 963145, 224193};
const npy_datetimestruct _FS_MIN_DTS = {1969, 12,     31,     21,    26,
                                        16,   627963, 145224, 193000};
const npy_datetimestruct _PS_MIN_DTS = {1969, 9,      16,     5, 57,
                                        7,    963145, 224193, 0};
const npy_datetimestruct _NS_MIN_DTS = {1677, 9,      21,     0, 12,
                                        43,   145224, 193000, 0};
const npy_datetimestruct _US_MIN_DTS = {-290308, 12,     21, 19, 59,
                                        05,      224193, 0,  0};
const npy_datetimestruct _MS_MIN_DTS = {-292275055, 5,      16, 16, 47,
                                        4,          193000, 0,  0};
const npy_datetimestruct _S_MIN_DTS = {
    -292277022657, 1, 27, 8, 29, 53, 0, 0, 0};
const npy_datetimestruct _M_MIN_DTS = {
    -17536621475646, 5, 4, 5, 53, 0, 0, 0, 0};

const npy_datetimestruct _AS_MAX_DTS = {1970, 1,      1,     0,     0,
                                        9,    223372, 36854, 775807};
const npy_datetimestruct _FS_MAX_DTS = {1970, 1,      1,      2,     33,
                                        43,   372036, 854775, 807000};
const npy_datetimestruct _PS_MAX_DTS = {1970, 4,     17,     18, 2,
                                        52,   36854, 775807, 0};
const npy_datetimestruct _NS_MAX_DTS = {2262, 4,      11,     23, 47,
                                        16,   854775, 807000, 0};
const npy_datetimestruct _US_MAX_DTS = {294247, 1, 10, 4, 0, 54, 775807, 0, 0};
const npy_datetimestruct _MS_MAX_DTS = {292278994, 8,      17, 7, 12,
                                        55,        807000, 0,  0};
const npy_datetimestruct _S_MAX_DTS = {292277026596, 12, 4, 15, 30, 7, 0, 0, 0};
const npy_datetimestruct _M_MAX_DTS = {
    17536621479585, 8, 30, 18, 7, 0, 0, 0, 0};

/* 'format_requirement' can be one of three values:
 *      * PARTIAL_MATCH : Only require a partial match with 'format'.
 *           For example, if the string is '2020-01-01 05:00:00' and
 *           'format' is '%Y-%m-%d', then parse '2020-01-01';
 *      * EXACT_MATCH : require an exact match with 'format'. If the
 *           string is '2020-01-01', then the only format which will
 *           be able to parse it without error is '%Y-%m-%d';
 *      * INFER_FORMAT: parse without comparing 'format' (i.e. infer it).
 */
typedef enum { PARTIAL_MATCH, EXACT_MATCH, INFER_FORMAT } FormatRequirement;

typedef struct {
  npy_datetime (*npy_datetimestruct_to_datetime)(NPY_DATETIMEUNIT,
                                                 const npy_datetimestruct *);
  int (*scaleNanosecToUnit)(npy_int64 *, NPY_DATETIMEUNIT);
  char *(*int64ToIso)(int64_t, NPY_DATETIMEUNIT, size_t *);
  npy_datetime (*NpyDateTimeToEpoch)(npy_datetime, NPY_DATETIMEUNIT);
  char *(*PyDateTimeToIso)(PyObject *, NPY_DATETIMEUNIT, size_t *);
  npy_datetime (*PyDateTimeToEpoch)(PyObject *, NPY_DATETIMEUNIT);
  char *(*int64ToIsoDuration)(int64_t, size_t *);
  void (*pandas_datetime_to_datetimestruct)(npy_datetime, NPY_DATETIMEUNIT,
                                            npy_datetimestruct *);
  void (*pandas_timedelta_to_timedeltastruct)(npy_datetime, NPY_DATETIMEUNIT,
                                              pandas_timedeltastruct *);
  int (*convert_pydatetime_to_datetimestruct)(PyObject *, npy_datetimestruct *);
  int (*cmp_npy_datetimestruct)(const npy_datetimestruct *,
                                const npy_datetimestruct *);
  PyArray_DatetimeMetaData (*get_datetime_metadata_from_dtype)(PyArray_Descr *);
  int (*parse_iso_8601_datetime)(const char *, int, int, npy_datetimestruct *,
                                 NPY_DATETIMEUNIT *, int *, int *, const char *,
                                 int, FormatRequirement);
  int (*get_datetime_iso_8601_strlen)(int, NPY_DATETIMEUNIT);
  int (*make_iso_8601_datetime)(npy_datetimestruct *, char *, int, int,
                                NPY_DATETIMEUNIT);
  int (*make_iso_8601_timedelta)(pandas_timedeltastruct *, char *, size_t *);
} PandasDateTime_CAPI;

// The capsule name appears limited to module.attributename; see bpo-32414
// cpython has an open PR gh-6898 to fix, but hasn't had traction for years
#define PandasDateTime_CAPSULE_NAME "pandas.pandas_datetime_CAPI"

/* block used as part of public API */
#ifndef _PANDAS_DATETIME_IMPL
static PandasDateTime_CAPI *PandasDateTimeAPI = NULL;

#define PandasDateTime_IMPORT                                                  \
  PandasDateTimeAPI =                                                          \
      (PandasDateTime_CAPI *)PyCapsule_Import(PandasDateTime_CAPSULE_NAME, 0)

#define npy_datetimestruct_to_datetime(NPY_DATETIMEUNIT, npy_datetimestruct)   \
  PandasDateTimeAPI->npy_datetimestruct_to_datetime((NPY_DATETIMEUNIT),        \
                                                    (npy_datetimestruct))
#define scaleNanosecToUnit(value, unit)                                        \
  PandasDateTimeAPI->scaleNanosecToUnit((value), (unit))
#define int64ToIso(value, base, len)                                           \
  PandasDateTimeAPI->int64ToIso((value), (base), (len))
#define NpyDateTimeToEpoch(dt, base)                                           \
  PandasDateTimeAPI->NpyDateTimeToEpoch((dt), (base))
#define PyDateTimeToIso(obj, base, len)                                        \
  PandasDateTimeAPI->PyDateTimeToIso((obj), (base), (len))
#define PyDateTimeToEpoch(dt, base)                                            \
  PandasDateTimeAPI->PyDateTimeToEpoch((dt), (base))
#define int64ToIsoDuration(value, len)                                         \
  PandasDateTimeAPI->int64ToIsoDuration((value), (len))
#define pandas_datetime_to_datetimestruct(dt, base, out)                       \
  PandasDateTimeAPI->pandas_datetime_to_datetimestruct((dt), (base), (out))
#define pandas_timedelta_to_timedeltastruct(td, base, out)                     \
  PandasDateTimeAPI->pandas_timedelta_to_timedeltastruct((td), (base), (out))
#define convert_pydatetime_to_datetimestruct(dtobj, out)                       \
  PandasDateTimeAPI->convert_pydatetime_to_datetimestruct((dtobj), (out))
#define cmp_npy_datetimestruct(a, b)                                           \
  PandasDateTimeAPI->cmp_npy_datetimestruct((a), (b))
#define get_datetime_metadata_from_dtype(dtype)                                \
  PandasDateTimeAPI->get_datetime_metadata_from_dtype((dtype))
#define parse_iso_8601_datetime(str, len, want_exc, out, out_bestunit,         \
                                out_local, out_tzoffset, format, format_len,   \
                                format_requirement)                            \
  PandasDateTimeAPI->parse_iso_8601_datetime(                                  \
      (str), (len), (want_exc), (out), (out_bestunit), (out_local),            \
      (out_tzoffset), (format), (format_len), (format_requirement))
#define get_datetime_iso_8601_strlen(local, base)                              \
  PandasDateTimeAPI->get_datetime_iso_8601_strlen((local), (base))
#define make_iso_8601_datetime(dts, outstr, outlen, utc, base)                 \
  PandasDateTimeAPI->make_iso_8601_datetime((dts), (outstr), (outlen), (utc),  \
                                            (base))
#define make_iso_8601_timedelta(tds, outstr, outlen)                           \
  PandasDateTimeAPI->make_iso_8601_timedelta((tds), (outstr), (outlen))
#endif /* !defined(_PANDAS_DATETIME_IMPL) */

#ifdef __cplusplus
}
#endif
#endif  // PANDAS__LIBS_TSLIBS_SRC_DATETIME_PD_DATETIME_H_
