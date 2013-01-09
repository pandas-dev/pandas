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

typedef int64_t i8;

enum CalendarType
{
    GREGORIAN,
    JULIAN
};

static const i8 SECONDS_PER_DAY = 86400;

#define Py_AssertWithArg(x, errortype, errorstr, a1)        \
    {                                                       \
        if (!((x))) {                                       \
            PyErr_Format((errortype), (errorstr), (a1));    \
            goto onError;                                   \
        }                                                   \
    }

#define Py_Error(errortype, errorstr)               \
    {                                               \
        PyErr_SetString((errortype), (errorstr));   \
        goto onError;                               \
    }


/*** FREQUENCY CONSTANTS ***/

// HIGHFREQ_ORIG is the datetime ordinal from which to begin the second
// frequency ordinal sequence

// typedef int64_t i8;
// begins second ordinal at 1/1/1970 unix epoch

// #define HIGHFREQ_ORIG 62135683200LL
static const i8 BASE_YEAR = 1970;
static const i8 ORD_OFFSET = 719163; // days until 1970-01-01
static const i8 BDAY_OFFSET = 513689; // days until 1970-01-01
static const i8 WEEK_OFFSET = 102737;
static const i8 HIGHFREQ_ORIG = 0; // ORD_OFFSET * 86400LL // days until 1970-01-01

enum Annual
{
    FR_ANN = 1000,  /* Annual */
    FR_ANNDEC = FR_ANN,  /* Annual - December year end*/
    FR_ANNJAN,  /* Annual - January year end*/
    FR_ANNFEB,  /* Annual - February year end*/
    FR_ANNMAR,  /* Annual - March year end*/
    FR_ANNAPR,  /* Annual - April year end*/
    FR_ANNMAY,  /* Annual - May year end*/
    FR_ANNJUN,  /* Annual - June year end*/
    FR_ANNJUL,  /* Annual - July year end*/
    FR_ANNAUG,  /* Annual - August year end*/
    FR_ANNSEP,  /* Annual - September year end*/
    FR_ANNOCT,  /* Annual - October year end*/
    FR_ANNNOV  /* Annual - November year end*/
};


/* The standard quarterly frequencies with various fiscal year ends
   eg, Q42005 for Q@OCT runs Aug 1, 2005 to Oct 31, 2005 */
enum Quarterly
{
    FR_QTR = 2000, /* Quarterly - December year end (default quarterly) */
    FR_QTRDEC = FR_QTR, /* Quarterly - December year end */
    FR_QTRJAN, /* Quarterly - January year end */
    FR_QTRFEB, /* Quarterly - February year end */
    FR_QTRMAR, /* Quarterly - March year end */
    FR_QTRAPR, /* Quarterly - April year end */
    FR_QTRMAY, /* Quarterly - May year end */
    FR_QTRJUN, /* Quarterly - June year end */
    FR_QTRJUL, /* Quarterly - July year end */
    FR_QTRAUG, /* Quarterly - August year end */
    FR_QTRSEP, /* Quarterly - September year end */
    FR_QTROCT, /* Quarterly - October year end */
    FR_QTRNOV /* Quarterly - November year end */
};

/* #define FR_MTH  3000  /\* Monthly *\/ */

enum Monthly
{
    FR_MTH = 3000
};


enum Weekly
{
    FR_WK = 4000,  /* Weekly */
    FR_WKSUN = FR_WK, /* Weekly - Sunday end of week */
    FR_WKMON, /* Weekly - Monday end of week */
    FR_WKTUE, /* Weekly - Tuesday end of week */
    FR_WKWED, /* Weekly - Wednesday end of week */
    FR_WKTHU, /* Weekly - Thursday end of week */
    FR_WKFRI, /* Weekly - Friday end of week */
    FR_WKSAT /* Weekly - Saturday end of week */
};

enum BusinessDaily { FR_BUS = 5000 };
enum Daily { FR_DAY = 6000 };
enum Hourly { FR_HR = 7000 };
enum Minutely { FR_MIN = 8000 };
enum Secondly { FR_SEC = 9000 };
enum Microsecondly { FR_USEC = 10000 };

enum Undefined { FR_UND = -10000 };

static const i8 US_PER_SECOND = 1000000L;
static const i8 US_PER_MINUTE = 60 * 1000000L;
static const i8 US_PER_HOUR = 60 * 60 * 1000000L;
static const i8 US_PER_DAY = 24 * 60 * 60 * 1000000L;
static const i8 US_PER_WEEK = 7 * 24 * 60 * 60 * 1000000L;

static const i8 NS_PER_SECOND = 1000000000L;
static const i8 NS_PER_MINUTE = 60 * 1000000000L;
static const i8 NS_PER_HOUR = 60 * 60 * 1000000000L;
static const i8 NS_PER_DAY = 24 * 60 * 60 * 1000000000L;
static const i8 NS_PER_WEEK = 7 * 24 * 60 * 60 * 1000000000L;

// make sure INT64_MIN is a macro!
static const i8 INT_ERR_CODE = INT64_MIN;

#define MEM_CHECK(item)                         \
    {                                           \
        if (item == NULL)                       \
            return PyErr_NoMemory();            \
    }

#define ERR_CHECK(item)                         \
    {                                           \
        if (item == NULL)                       \
            return NULL;                        \
    }

typedef struct asfreq_info
{
    i8 from_week_end;   // day the week ends on in the "from" frequency
    i8 to_week_end;     // day the week ends on in the "to" frequency

    i8 from_a_year_end; // month the year ends on in the "from" frequency
    i8 to_a_year_end;   // month the year ends on in the "to" frequency

    i8 from_q_year_end; // month the year ends on in the "from" frequency
    i8 to_q_year_end;   // month the year ends on in the "to" frequency
} asfreq_info;


typedef struct date_info
{
    i8 absdate;
    i8 abstime;

    i8 attosecond;
    i8 femtosecond;
    i8 picosecond;
    i8 nanosecond;
    i8 microsecond;
    i8 second;
    i8 minute;
    i8 hour;
    i8 day;
    i8 month;
    i8 quarter;
    i8 year;
    i8 day_of_week;
    i8 day_of_year;
    i8 calendar;
} date_info;

typedef i8 (*freq_conv_func)(i8, const char*, asfreq_info*);

/*
 * new pandas API helper functions here
 */

i8 asfreq(i8 period_ordinal, i8 freq1, i8 freq2, const char* relation);

i8 get_period_ordinal(i8 year, i8 month, i8 day, i8 hour, i8 minute, i8 second,
                      i8 microsecond, i8 freq);

i8 get_python_ordinal(i8 period_ordinal, i8 freq);

i8 get_date_info(i8 ordinal, i8 freq, struct date_info *dinfo);
freq_conv_func get_asfreq_func(i8 fromFreq, i8 toFreq);
void get_asfreq_info(i8 fromFreq, i8 toFreq, asfreq_info *af_info);

i8 pyear(i8 ordinal, i8 freq);
i8 pqyear(i8 ordinal, i8 freq);
i8 pquarter(i8 ordinal, i8 freq);
i8 pmonth(i8 ordinal, i8 freq);
i8 pday(i8 ordinal, i8 freq);
i8 pweekday(i8 ordinal, i8 freq);
i8 pday_of_week(i8 ordinal, i8 freq);
i8 pday_of_year(i8 ordinal, i8 freq);
i8 pweek(i8 ordinal, i8 freq);
i8 phour(i8 ordinal, i8 freq);
i8 pminute(i8 ordinal, i8 freq);
i8 psecond(i8 ordinal, i8 freq);
i8 pmicrosecond(i8 ordinal, i8 freq);
i8 pnanosecond(i8 ordinal, i8 freq);
i8 ppicosecond(i8 ordinal, i8 freq);
i8 pfemtosecond(i8 ordinal, i8 freq);
i8 pattosecond(i8 ordinal, i8 freq);

char *c_strftime(struct date_info *dinfo, char *fmt);
i8 get_yq(i8 ordinal, i8 freq, i8 *quarter, i8 *year);

#endif
