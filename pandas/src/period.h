/*
 * Borrowed and derived code from scikits.timeseries that we will expose via
 * Cython to pandas. This primarily concerns interval representation and
 * frequency conversion routines.
 */

#ifndef C_PERIOD_H
#define C_PERIOD_H

#include <Python.h>
#include "numpy/ndarraytypes.h"
#include "headers/stdint.h"
#include "limits.h"

/*
 * declarations from period here
 */

#define GREGORIAN_CALENDAR 0
#define JULIAN_CALENDAR 1

#define SECONDS_PER_DAY ((double) 86400.0)

#define Py_AssertWithArg(x,errortype,errorstr,a1) {if (!(x)) {PyErr_Format(errortype,errorstr,a1);goto onError;}}
#define Py_Error(errortype,errorstr) {PyErr_SetString(errortype,errorstr);goto onError;}

/*** FREQUENCY CONSTANTS ***/

// HIGHFREQ_ORIG is the datetime ordinal from which to begin the second
// frequency ordinal sequence

// typedef int64_t npy_int64;
// begins second ordinal at 1/1/1970 unix epoch

// #define HIGHFREQ_ORIG 62135683200LL
#define BASE_YEAR 1970
#define ORD_OFFSET 719163LL // days until 1970-01-01
#define BDAY_OFFSET 513689LL // days until 1970-01-01
#define WEEK_OFFSET 102737LL
#define HIGHFREQ_ORIG 0 // ORD_OFFSET * 86400LL // days until 1970-01-01

#define FR_ANN  1000  /* Annual */
#define FR_ANNDEC  FR_ANN  /* Annual - December year end*/
#define FR_ANNJAN  1001  /* Annual - January year end*/
#define FR_ANNFEB  1002  /* Annual - February year end*/
#define FR_ANNMAR  1003  /* Annual - March year end*/
#define FR_ANNAPR  1004  /* Annual - April year end*/
#define FR_ANNMAY  1005  /* Annual - May year end*/
#define FR_ANNJUN  1006  /* Annual - June year end*/
#define FR_ANNJUL  1007  /* Annual - July year end*/
#define FR_ANNAUG  1008  /* Annual - August year end*/
#define FR_ANNSEP  1009  /* Annual - September year end*/
#define FR_ANNOCT  1010  /* Annual - October year end*/
#define FR_ANNNOV  1011  /* Annual - November year end*/

/* The standard quarterly frequencies with various fiscal year ends
   eg, Q42005 for Q@OCT runs Aug 1, 2005 to Oct 31, 2005 */
#define FR_QTR  2000       /* Quarterly - December year end (default quarterly) */
#define FR_QTRDEC  FR_QTR  /* Quarterly - December year end */
#define FR_QTRJAN  2001    /* Quarterly - January year end */
#define FR_QTRFEB  2002    /* Quarterly - February year end */
#define FR_QTRMAR  2003    /* Quarterly - March year end */
#define FR_QTRAPR  2004    /* Quarterly - April year end */
#define FR_QTRMAY  2005    /* Quarterly - May year end */
#define FR_QTRJUN  2006    /* Quarterly - June year end */
#define FR_QTRJUL  2007    /* Quarterly - July year end */
#define FR_QTRAUG  2008    /* Quarterly - August year end */
#define FR_QTRSEP  2009    /* Quarterly - September year end */
#define FR_QTROCT  2010    /* Quarterly - October year end */
#define FR_QTRNOV  2011    /* Quarterly - November year end */

#define FR_MTH  3000  /* Monthly */

#define FR_WK   4000  /* Weekly */
#define FR_WKSUN FR_WK /* Weekly - Sunday end of week */
#define FR_WKMON 4001 /* Weekly - Monday end of week */
#define FR_WKTUE 4002 /* Weekly - Tuesday end of week */
#define FR_WKWED 4003 /* Weekly - Wednesday end of week */
#define FR_WKTHU 4004 /* Weekly - Thursday end of week */
#define FR_WKFRI 4005 /* Weekly - Friday end of week */
#define FR_WKSAT 4006 /* Weekly - Saturday end of week */

#define FR_BUS  5000  /* Business days */
#define FR_DAY  6000  /* Daily */
#define FR_HR   7000  /* Hourly */
#define FR_MIN  8000  /* Minutely */
#define FR_SEC  9000  /* Secondly */

#define FR_UND  -10000 /* Undefined */

#define INT_ERR_CODE INT32_MIN

#define MEM_CHECK(item) if (item == NULL) { return PyErr_NoMemory(); }
#define ERR_CHECK(item) if (item == NULL) { return NULL; }

typedef struct asfreq_info {
    int from_week_end;   // day the week ends on in the "from" frequency
    int to_week_end;     // day the week ends on in the "to" frequency

    int from_a_year_end; // month the year ends on in the "from" frequency
    int to_a_year_end;   // month the year ends on in the "to" frequency

    int from_q_year_end; // month the year ends on in the "from" frequency
    int to_q_year_end;   // month the year ends on in the "to" frequency
} asfreq_info;


typedef struct date_info {
    npy_int64 absdate;
    double abstime;

    double second;
    int minute;
    int hour;
    int day;
    int month;
    int quarter;
    int year;
    int day_of_week;
    int day_of_year;
    int calendar;
} date_info;

typedef npy_int64 (*freq_conv_func)(npy_int64, char, asfreq_info*);

/*
 * new pandas API helper functions here
 */

npy_int64 asfreq(npy_int64 period_ordinal, int freq1, int freq2, char relation);

npy_int64 get_period_ordinal(int year, int month, int day,
                      int hour, int minute, int second,
                      int freq);

npy_int64 get_python_ordinal(npy_int64 period_ordinal, int freq);

int get_date_info(npy_int64 ordinal, int freq, struct date_info *dinfo);
freq_conv_func get_asfreq_func(int fromFreq, int toFreq);
void get_asfreq_info(int fromFreq, int toFreq, asfreq_info *af_info);

int pyear(npy_int64 ordinal, int freq);
int pqyear(npy_int64 ordinal, int freq);
int pquarter(npy_int64 ordinal, int freq);
int pmonth(npy_int64 ordinal, int freq);
int pday(npy_int64 ordinal, int freq);
int pweekday(npy_int64 ordinal, int freq);
int pday_of_week(npy_int64 ordinal, int freq);
int pday_of_year(npy_int64 ordinal, int freq);
int pweek(npy_int64 ordinal, int freq);
int phour(npy_int64 ordinal, int freq);
int pminute(npy_int64 ordinal, int freq);
int psecond(npy_int64 ordinal, int freq);

double getAbsTime(int freq, npy_int64 dailyDate, npy_int64 originalDate);
char *c_strftime(struct date_info *dinfo, char *fmt);
int get_yq(npy_int64 ordinal, int freq, int *quarter, int *year);

#endif
