/* 
 * This file will contain borrowed and derived code from scikits.timeseries
 * that we will expose via Cython.
 */

#ifndef C_SKTS_H
#define C_SKTS_H

#include <Python.h>

/*
 * declarations from skts here
 */

#define GREGORIAN_CALENDAR 0
#define JULIAN_CALENDAR 1

#define SECONDS_PER_DAY ((double) 86400.0)

#define Py_AssertWithArg(x,errortype,errorstr,a1) {if (!(x)) {PyErr_Format(errortype,errorstr,a1);goto onError;}}
#define Py_Error(errortype,errorstr) {PyErr_SetString(errortype,errorstr);goto onError;}

/*** FREQUENCY CONSTANTS ***/

// datetime ordinal of unix epoch
#define HIGHFREQ_ORIG 719163

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

/* The standard quarterly frequencies. Year is determined by what year the end
   month lies in. */
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

/* End period based quarterly frequencies. Year is determined by what year the
   end month lies in. */
#define FR_QTREDEC  FR_QTRDEC  /* Quarterly - December year end*/
#define FR_QTREJAN  FR_QTRJAN  /* Quarterly - January year end*/
#define FR_QTREFEB  FR_QTRFEB  /* Quarterly - February year end*/
#define FR_QTREMAR  FR_QTRMAR  /* Quarterly - March year end*/
#define FR_QTREAPR  FR_QTRAPR  /* Quarterly - April year end*/
#define FR_QTREMAY  FR_QTRMAY  /* Quarterly - May year end*/
#define FR_QTREJUN  FR_QTRJUN  /* Quarterly - June year end*/
#define FR_QTREJUL  FR_QTRJUL  /* Quarterly - July year end*/
#define FR_QTREAUG  FR_QTRAUG  /* Quarterly - August year end*/
#define FR_QTRESEP  FR_QTRSEP  /* Quarterly - September year end*/
#define FR_QTREOCT  FR_QTROCT  /* Quarterly - October year end*/
#define FR_QTRENOV  FR_QTRNOV  /* Quarterly - November year end*/

/* Starting period based quarterly frequencies. Year is determined by what year
   the starting month lies in. */
#define FR_QTRSDEC  FR_QTRDEC+12  /* Quarterly - December year end*/
#define FR_QTRSJAN  FR_QTRJAN+12  /* Quarterly - January year end*/
#define FR_QTRSFEB  FR_QTRFEB+12  /* Quarterly - February year end*/
#define FR_QTRSMAR  FR_QTRMAR+12  /* Quarterly - March year end*/
#define FR_QTRSAPR  FR_QTRAPR+12  /* Quarterly - April year end*/
#define FR_QTRSMAY  FR_QTRMAY+12  /* Quarterly - May year end*/
#define FR_QTRSJUN  FR_QTRJUN+12  /* Quarterly - June year end*/
#define FR_QTRSJUL  FR_QTRJUL+12  /* Quarterly - July year end*/
#define FR_QTRSAUG  FR_QTRAUG+12  /* Quarterly - August year end*/
#define FR_QTRSSEP  FR_QTRSEP+12  /* Quarterly - September year end*/
#define FR_QTRSOCT  FR_QTROCT+12  /* Quarterly - October year end*/
#define FR_QTRSNOV  FR_QTRNOV+12  /* Quarterly - November year end*/

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

#define INT_ERR_CODE -1

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
    long absdate;
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

typedef long (*freq_conv_func)(long, char, asfreq_info*);

/*
 * new pandas API helper functions here
 */

int frequency_conversion(long dtordinal, int freq1, int freq2, char relation);

long get_skts_ordinal(int year, int month, int day,
                      int hour, int minute, int second,
                      int freq);

long get_python_ordinal(long skts_ordinal, int freq);

#endif
