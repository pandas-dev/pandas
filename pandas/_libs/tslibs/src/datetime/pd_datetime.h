/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

*/

#ifndef PANDAS__LIBS_TSLIBS_SRC_DATETIME_PD_DATETIME_H
#define PANDAS__LIBS_TSLIBS_SRC_DATETIME_PD_DATETIME_H


#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif  // NPY_NO_DEPRECATED_API

#include <numpy/ndarraytypes.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    npy_datetime (*npy_datetimestruct_to_datetime)(NPY_DATETIMEUNIT, const npy_datetimestruct *);
  } PandasDateTime_CAPI;

  #define PandasDateTime_CAPSULE_NAME "pandas._libs.tslibs.src.datetime.pd_datetime_CAPI"

  /* block used as part of public API */
#ifndef _PANDAS_DATETIME_IMPL
  static PandasDateTime_CAPI *PandasDateTimeAPI = NULL;

  #define PandasDateTime_IMPORT \
    PandasDateTimeAPI = (PandasDateTime_CAPI *)PyCapsule_Import(PandasDateTime_CAPSULE_NAME, 0)
#endif /* !defined(_PANDAS_DATETIME_IMPL) */

#ifdef __cplusplus
}
#endif
#endif /* !defined(PANDAS__LIBS_TSLIBS_SRC_DATETIME_PD_DATETIME_H) */
