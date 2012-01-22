/*
 * NB: This is derived from numpy 1.7 datetime.c, just enough code to
 * do some conversions. Copyrights from that file apply.
 */

#ifndef _PANDAS_DATETIME_H_
#define _PANDAS_DATETIME_H_

int convert_pydatetime_to_datetimestruct(PyObject *obj, npy_datetimestruct *out,
                                         NPY_DATETIMEUNIT *out_bestunit,
                                         int apply_tzinfo);

#endif
