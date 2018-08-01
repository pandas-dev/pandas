/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Borrowed and derived code from scikits.timeseries that we will expose via
Cython to pandas. This primarily concerns interval representation and
frequency conversion routines.
*/

#ifndef PANDAS__LIBS_SRC_PERIOD_HELPER_H_
#define PANDAS__LIBS_SRC_PERIOD_HELPER_H_

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

#include <Python.h>
//#include "../../src/headers/stdint.h"
#include "limits.h"
#include "numpy/ndarraytypes.h"

/*** FREQUENCY CONSTANTS ***/

#define FR_ANN 1000      /* Annual */
#define FR_ANNDEC FR_ANN /* Annual - December year end*/
#define FR_ANNJAN 1001   /* Annual - January year end*/
#define FR_ANNFEB 1002   /* Annual - February year end*/
#define FR_ANNMAR 1003   /* Annual - March year end*/
#define FR_ANNAPR 1004   /* Annual - April year end*/
#define FR_ANNMAY 1005   /* Annual - May year end*/
#define FR_ANNJUN 1006   /* Annual - June year end*/
#define FR_ANNJUL 1007   /* Annual - July year end*/
#define FR_ANNAUG 1008   /* Annual - August year end*/
#define FR_ANNSEP 1009   /* Annual - September year end*/
#define FR_ANNOCT 1010   /* Annual - October year end*/
#define FR_ANNNOV 1011   /* Annual - November year end*/

/* The standard quarterly frequencies with various fiscal year ends
   eg, Q42005 for Q@OCT runs Aug 1, 2005 to Oct 31, 2005 */
#define FR_QTR 2000      /* Quarterly - December year end (default quarterly) */
#define FR_QTRDEC FR_QTR /* Quarterly - December year end */
#define FR_QTRJAN 2001   /* Quarterly - January year end */
#define FR_QTRFEB 2002   /* Quarterly - February year end */
#define FR_QTRMAR 2003   /* Quarterly - March year end */
#define FR_QTRAPR 2004   /* Quarterly - April year end */
#define FR_QTRMAY 2005   /* Quarterly - May year end */
#define FR_QTRJUN 2006   /* Quarterly - June year end */
#define FR_QTRJUL 2007   /* Quarterly - July year end */
#define FR_QTRAUG 2008   /* Quarterly - August year end */
#define FR_QTRSEP 2009   /* Quarterly - September year end */
#define FR_QTROCT 2010   /* Quarterly - October year end */
#define FR_QTRNOV 2011   /* Quarterly - November year end */

#define FR_MTH 3000 /* Monthly */

#define FR_WK 4000     /* Weekly */
#define FR_WKSUN FR_WK /* Weekly - Sunday end of week */
#define FR_WKMON 4001  /* Weekly - Monday end of week */
#define FR_WKTUE 4002  /* Weekly - Tuesday end of week */
#define FR_WKWED 4003  /* Weekly - Wednesday end of week */
#define FR_WKTHU 4004  /* Weekly - Thursday end of week */
#define FR_WKFRI 4005  /* Weekly - Friday end of week */
#define FR_WKSAT 4006  /* Weekly - Saturday end of week */

#define FR_BUS 5000 /* Business days */
#define FR_DAY 6000 /* Daily */
#define FR_HR 7000  /* Hourly */
#define FR_MIN 8000 /* Minutely */
#define FR_SEC 9000 /* Secondly */
#define FR_MS 10000 /* Millisecondly */
#define FR_US 11000 /* Microsecondly */
#define FR_NS 12000 /* Nanosecondly */

#define FR_UND -10000 /* Undefined */

#define INT_ERR_CODE NPY_MIN_INT32

typedef struct asfreq_info {
    int is_end;
    // char relation == 'S' (for START) --> is_end = 0
    // char relation == 'E' (for END) --> is_end = 1

    int from_end;
    int to_end;
    // weekly:
    // from_end --> day the week ends on in the "from" frequency
    // to_end   --> day the week ends on in the "to" frequency
    //
    // annual:
    // from_end --> month the year ends on in the "from" frequency
    // to_end   --> month the year ends on in the "to" frequency
    //
    // quarterly:
    // from_end --> month the year ends on in the "from" frequency
    // to_end   --> month the year ends on in the "to" frequency

    npy_int64 intraday_conversion_factor;
} asfreq_info;

typedef npy_int64 (*freq_conv_func)(npy_int64, asfreq_info *af_info);

/*
 * new pandas API helper functions here
 */

freq_conv_func get_asfreq_func(int fromFreq, int toFreq);

npy_int64 get_daytime_conversion_factor(int from_index, int to_index);
int max_value(int a, int b);

#endif  // PANDAS__LIBS_SRC_PERIOD_HELPER_H_
