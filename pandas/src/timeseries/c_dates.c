#include "c_dates.h"
#include <datetime.h>
#include <time.h>


int get_freq_group(int freq) { return (freq/1000)*1000; }

static asfreq_info NULL_AF_INFO;

/*********************************************************
** Python callbacks. These functions must be called by  **
** the module __init__ script                           **
*********************************************************/

static PyObject *DateFromString = NULL;
PyObject *
set_callback_DateFromString(PyObject *dummy, PyObject *args) {
    return set_callback(args, &DateFromString);
}

static PyObject *DateTimeFromString = NULL;
PyObject *
set_callback_DateTimeFromString(PyObject *dummy, PyObject *args) {
    return set_callback(args, &DateTimeFromString);
}

//DERIVED FROM mx.DateTime
/*
    Functions in the following section are borrowed from mx.DateTime version
    2.0.6, and hence this code is subject to the terms of the egenix public
    license version 1.0.0
*/

#define Py_AssertWithArg(x,errortype,errorstr,a1) {if (!(x)) {PyErr_Format(errortype,errorstr,a1);goto onError;}}
#define Py_Error(errortype,errorstr) {PyErr_SetString(errortype,errorstr);goto onError;}

 /* Error Exception objects */
static PyObject *DateCalc_Error;
static PyObject *DateCalc_RangeError;

#define GREGORIAN_CALENDAR 0
#define JULIAN_CALENDAR 1

#define SECONDS_PER_DAY ((double) 86400.0)

/* Table with day offsets for each month (0-based, without and with leap) */
static int month_offset[2][13] = {
    { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 },
    { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 }
};

/* Table of number of days in a month (0-based, without and with leap) */
static int days_in_month[2][12] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

struct date_info {
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
};


/* Return 1/0 iff year points to a leap year in calendar. */
static
int dInfoCalc_Leapyear(register long year,
            int calendar)
{
    if (calendar == GREGORIAN_CALENDAR) {
        return (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));
    } else {
        return (year % 4 == 0);
    }
}

static
int dInfoCalc_ISOWeek(struct date_info *dinfo)
{
    int week;

    /* Estimate */
    week = (dinfo->day_of_year-1) - dinfo->day_of_week + 3;
    if (week >= 0) week = week / 7 + 1;

    /* Verify */
    if (week < 0) {
        /* The day lies in last week of the previous year */
        if ((week > -2) ||
            (week == -2 && dInfoCalc_Leapyear(dinfo->year-1, dinfo->calendar)))
            week = 53;
        else
            week = 52;
    } else if (week == 53) {
    /* Check if the week belongs to year or year+1 */
        if (31-dinfo->day + dinfo->day_of_week < 3) {
            week = 1;
        }
    }

    return week;
}


/* Return the day of the week for the given absolute date. */
static
int dInfoCalc_DayOfWeek(register long absdate)
{
    int day_of_week;

    if (absdate >= 1) {
        day_of_week = (absdate - 1) % 7;
    } else {
        day_of_week = 6 - ((-absdate) % 7);
    }
    return day_of_week;
}

/* Return the year offset, that is the absolute date of the day
   31.12.(year-1) in the given calendar.

   Note:
   For the Julian calendar we shift the absdate (which is measured
   using the Gregorian Epoch) value by two days because the Epoch
   (0001-01-01) in the Julian calendar lies 2 days before the Epoch in
   the Gregorian calendar. */
static
int dInfoCalc_YearOffset(register long year,
              int calendar)
{
    year--;
    if (calendar == GREGORIAN_CALENDAR) {
    if (year >= 0 || -1/4 == -1)
        return year*365 + year/4 - year/100 + year/400;
    else
        return year*365 + (year-3)/4 - (year-99)/100 + (year-399)/400;
    }
    else if (calendar == JULIAN_CALENDAR) {
    if (year >= 0 || -1/4 == -1)
        return year*365 + year/4 - 2;
    else
        return year*365 + (year-3)/4 - 2;
    }
    Py_Error(DateCalc_Error, "unknown calendar");
 onError:
    return -1;
}


/* Set the instance's value using the given date and time. calendar
   may be set to the flags: GREGORIAN_CALENDAR,
   JULIAN_CALENDAR to indicate the calendar to be used. */

static
int dInfoCalc_SetFromDateAndTime(struct date_info *dinfo,
                  int year,
                  int month,
                  int day,
                  int hour,
                  int minute,
                  double second,
                  int calendar)
{

    /* Calculate the absolute date */
    {
        int leap;
        long yearoffset,absdate;

        /* Range check */
        Py_AssertWithArg(year > -(INT_MAX / 366) && year < (INT_MAX / 366),
                 DateCalc_RangeError,
                 "year out of range: %i",
                 year);

        /* Is it a leap year ? */
        leap = dInfoCalc_Leapyear(year,calendar);

        /* Negative month values indicate months relative to the years end */
        if (month < 0) month += 13;
        Py_AssertWithArg(month >= 1 && month <= 12,
                 DateCalc_RangeError,
                 "month out of range (1-12): %i",
                 month);

        /* Negative values indicate days relative to the months end */
        if (day < 0) day += days_in_month[leap][month - 1] + 1;
        Py_AssertWithArg(day >= 1 && day <= days_in_month[leap][month - 1],
                 DateCalc_RangeError,
                 "day out of range: %i",
                 day);

        yearoffset = dInfoCalc_YearOffset(year,calendar);
        if (PyErr_Occurred()) goto onError;

        absdate = day + month_offset[leap][month - 1] + yearoffset;

        dinfo->absdate = absdate;

        dinfo->year = year;
        dinfo->month = month;
        dinfo->quarter = ((month-1)/3)+1;
        dinfo->day = day;

        dinfo->day_of_week = dInfoCalc_DayOfWeek(absdate);
        dinfo->day_of_year = (short)(absdate - yearoffset);

        dinfo->calendar = calendar;
    }

    /* Calculate the absolute time */
    {
    Py_AssertWithArg(hour >= 0 && hour <= 23,
             DateCalc_RangeError,
             "hour out of range (0-23): %i",
             hour);
    Py_AssertWithArg(minute >= 0 && minute <= 59,
             DateCalc_RangeError,
             "minute out of range (0-59): %i",
             minute);
    Py_AssertWithArg(second >= (double)0.0 &&
             (second < (double)60.0 ||
              (hour == 23 && minute == 59 &&
               second < (double)61.0)),
             DateCalc_RangeError,
             "second out of range (0.0 - <60.0; <61.0 for 23:59): %f",
             second);

    dinfo->abstime = (double)(hour*3600 + minute*60) + second;

    dinfo->hour = hour;
    dinfo->minute = minute;
    dinfo->second = second;
    }
    return 0;
 onError:
    return -1;
}

static int monthToQuarter(int month) { return ((month-1)/3)+1; }

/* Sets the date part of the date_info struct using the indicated
   calendar.

   XXX This could also be done using some integer arithmetics rather
       than with this iterative approach... */
static
int dInfoCalc_SetFromAbsDate(register struct date_info *dinfo,
                  long absdate,
                  int calendar)
{
    register long year;
    long yearoffset;
    int leap,dayoffset;
    int *monthoffset;

    /* Approximate year */
    if (calendar == GREGORIAN_CALENDAR) {
        year = (long)(((double)absdate) / 365.2425);
    } else if (calendar == JULIAN_CALENDAR) {
        year = (long)(((double)absdate) / 365.25);
    } else {
        Py_Error(DateCalc_Error, "unknown calendar");
    }
    if (absdate > 0) year++;

    /* Apply corrections to reach the correct year */
    while (1) {
        /* Calculate the year offset */
        yearoffset = dInfoCalc_YearOffset(year,calendar);
        if (PyErr_Occurred())
            goto onError;

        /* Backward correction: absdate must be greater than the
           yearoffset */
        if (yearoffset >= absdate) {
            year--;
            continue;
        }

        dayoffset = absdate - yearoffset;
        leap = dInfoCalc_Leapyear(year,calendar);

        /* Forward correction: non leap years only have 365 days */
        if (dayoffset > 365 && !leap) {
            year++;
            continue;
        }
        break;
    }

    dinfo->year = year;
    dinfo->calendar = calendar;

    /* Now iterate to find the month */
    monthoffset = month_offset[leap];
    {
        register int month;

        for (month = 1; month < 13; month++) {
            if (monthoffset[month] >= dayoffset)
            break;
        }

        dinfo->month = month;
        dinfo->quarter = monthToQuarter(month);
        dinfo->day = dayoffset - month_offset[leap][month-1];
    }


    dinfo->day_of_week = dInfoCalc_DayOfWeek(absdate);
    dinfo->day_of_year = dayoffset;
    dinfo->absdate = absdate;

    return 0;

 onError:
    return -1;
}

/* Sets the time part of the DateTime object. */
static
int dInfoCalc_SetFromAbsTime(struct date_info *dinfo,
                  double abstime)
{
    int inttime;
    int hour,minute;
    double second;

    inttime = (int)abstime;
    hour = inttime / 3600;
    minute = (inttime % 3600) / 60;
    second = abstime - (double)(hour*3600 + minute*60);

    dinfo->hour = hour;
    dinfo->minute = minute;
    dinfo->second = second;

    dinfo->abstime = abstime;

    return 0;
}

/* Set the instance's value using the given date and time. calendar
   may be set to the flags: GREGORIAN_CALENDAR, JULIAN_CALENDAR to
   indicate the calendar to be used. */
static
int dInfoCalc_SetFromAbsDateTime(struct date_info *dinfo,
                  long absdate,
                  double abstime,
                  int calendar)
{

    /* Bounds check */
    Py_AssertWithArg(abstime >= 0.0 && abstime <= SECONDS_PER_DAY,
             DateCalc_Error,
             "abstime out of range (0.0 - 86400.0): %f",
             abstime);

    /* Calculate the date */
    if (dInfoCalc_SetFromAbsDate(dinfo,
                  absdate,
                  calendar))
    goto onError;

    /* Calculate the time */
    if (dInfoCalc_SetFromAbsTime(dinfo,
                  abstime))
    goto onError;

    return 0;
 onError:
    return -1;
}

/*
====================================================
== End of section borrowed from mx.DateTime       ==
====================================================
*/





///////////////////////////////////////////////////////////////////////

// helpers for frequency conversion routines //

static long DtoB_weekday(long fromDate) { return (((fromDate) / 7) * 5) + (fromDate)%7; }

static long DtoB_WeekendToMonday(long absdate, int day_of_week) {

    if (day_of_week > 4) {
        //change to Monday after weekend
        absdate += (7 - day_of_week);
    }
    return DtoB_weekday(absdate);
}

static long DtoB_WeekendToFriday(long absdate, int day_of_week) {

    if (day_of_week > 4) {
        //change to friday before weekend
        absdate -= (day_of_week - 4);
    }
    return DtoB_weekday(absdate);
}

static long absdate_from_ymd(int y, int m, int d) {
    struct date_info tempDate;
    if (dInfoCalc_SetFromDateAndTime(&tempDate, y, m, d, 0, 0, 0, GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    return tempDate.absdate;
}


///////////////////////////////////////////////

// frequency specifc conversion routines
// each function must take an integer fromDate and a char relation ('S' or 'E' for 'START' or 'END')

//************ FROM DAILY ***************

static long asfreq_DtoA(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    if (dinfo.month > af_info->to_a_year_end) { return (long)(dinfo.year + 1); }
    else { return (long)(dinfo.year); }
}

static long DtoQ_yq(long fromDate, asfreq_info *af_info,
                              int *year, int *quarter) {
    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    if (af_info->to_q_year_end != 12) {
        dinfo.month -= af_info->to_q_year_end;
        if (dinfo.month <= 0) { dinfo.month += 12; }
        else { dinfo.year += 1; }
        dinfo.quarter = monthToQuarter(dinfo.month);
    }

    *year = dinfo.year;
    *quarter = dinfo.quarter;

    return 0;
}


static long asfreq_DtoQ(long fromDate, char relation, asfreq_info *af_info) {

    int year, quarter;

    if (DtoQ_yq(fromDate, af_info, &year, &quarter) == INT_ERR_CODE)
    { return INT_ERR_CODE; }

    return (long)((year - 1) * 4 + quarter);
}

static long asfreq_DtoM(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    return (long)((dinfo.year - 1) * 12 + dinfo.month);
}

static long asfreq_DtoW(long fromDate, char relation, asfreq_info *af_info) {
    return (fromDate - (1 + af_info->to_week_end))/7 + 1;
}

static long asfreq_DtoB(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') {
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
    } else {
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
    }
}

static long asfreq_DtoB_forConvert(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (dinfo.day_of_week > 4) {
        return -1;
    } else {
        return DtoB_weekday(fromDate);
    }
}

// needed for getDateInfo function
static long asfreq_DtoD(long fromDate, char relation, asfreq_info *af_info) { return fromDate; }

static long asfreq_DtoHIGHFREQ(long fromDate, char relation, long periodsPerDay) {
    if (fromDate >= HIGHFREQ_ORIG) {
        if (relation == 'S') { return (fromDate - HIGHFREQ_ORIG)*(periodsPerDay) + 1; }
        else                 { return (fromDate - HIGHFREQ_ORIG + 1)*(periodsPerDay); }
    } else { return -1; }
}

static long asfreq_DtoH(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24); }
static long asfreq_DtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24*60); }
static long asfreq_DtoS(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(fromDate, relation, 24*60*60); }

//************ FROM SECONDLY ***************

static long asfreq_StoD(long fromDate, char relation, asfreq_info *af_info)
    { return (fromDate - 1)/(60*60*24) + HIGHFREQ_ORIG; }

static long asfreq_StoA(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_StoQ(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_StoM(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_StoB(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoB_forConvert(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoB_forConvert(asfreq_StoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_StoT(long fromDate, char relation, asfreq_info *af_info)
    { return (fromDate - 1)/60 + 1; }
static long asfreq_StoH(long fromDate, char relation, asfreq_info *af_info)
    { return (fromDate - 1)/(60*60) + 1; }

//************ FROM MINUTELY ***************

static long asfreq_TtoD(long fromDate, char relation, asfreq_info *af_info)
    { return (fromDate - 1)/(60*24) + HIGHFREQ_ORIG; }

static long asfreq_TtoA(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_TtoQ(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_TtoM(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_TtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_TtoB(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_TtoB_forConvert(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoB_forConvert(asfreq_TtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_TtoH(long fromDate, char relation, asfreq_info *af_info)
    { return (fromDate - 1)/60 + 1; }
static long asfreq_TtoS(long fromDate, char relation, asfreq_info *af_info) {
    if (relation == 'S') {  return fromDate*60 - 59; }
    else                 {  return fromDate*60;      }}

//************ FROM HOURLY ***************

static long asfreq_HtoD(long fromDate, char relation, asfreq_info *af_info)
    { return (fromDate - 1)/24 + HIGHFREQ_ORIG; }
static long asfreq_HtoA(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_HtoQ(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_HtoM(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_HtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }
static long asfreq_HtoB(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_HtoB_forConvert(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoB_forConvert(asfreq_HtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

// calculation works out the same as TtoS, so we just call that function for HtoT
static long asfreq_HtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_TtoS(fromDate, relation, &NULL_AF_INFO); }
static long asfreq_HtoS(long fromDate, char relation, asfreq_info *af_info) {
    if (relation == 'S') {  return fromDate*60*60 - 60*60 + 1; }
    else                 {  return fromDate*60*60;             }}

//************ FROM BUSINESS ***************

static long asfreq_BtoD(long fromDate, char relation, asfreq_info *af_info)
    { return ((fromDate-1)/5)*7 + (fromDate-1)%5 + 1; }

static long asfreq_BtoA(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_BtoQ(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_BtoM(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_BtoH(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static long asfreq_BtoS(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_BtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM WEEKLY ***************

static long asfreq_WtoD(long fromDate, char relation, asfreq_info *af_info) {
    if (relation == 'S') { return fromDate * 7 - 6 + af_info->from_week_end;}
    else                 { return fromDate * 7 + af_info->from_week_end; }
}

static long asfreq_WtoA(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_WtoD(fromDate, 'E', af_info), relation, af_info); }
static long asfreq_WtoQ(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_WtoD(fromDate, 'E', af_info), relation, af_info); }
static long asfreq_WtoM(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_WtoD(fromDate, 'E', af_info), relation, &NULL_AF_INFO); }

static long asfreq_WtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_WtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_WtoB(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_WtoD(fromDate, relation, af_info),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static long asfreq_WtoH(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_WtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_WtoS(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_WtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }

//************ FROM MONTHLY ***************

static void MtoD_ym(long fromDate, long *y, long *m) {
    *y = (fromDate - 1) / 12 + 1;
    *m = fromDate - 12 * (*y) - 1;
}

static long asfreq_MtoD(long fromDate, char relation, asfreq_info *af_info) {

    long y, m, absdate;

    if (relation == 'S') {
        MtoD_ym(fromDate, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate;
    } else {
        MtoD_ym(fromDate+1, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate-1;
    }
}

static long asfreq_MtoA(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_MtoD(fromDate, 'E', &NULL_AF_INFO), relation, af_info); }

static long asfreq_MtoQ(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_MtoD(fromDate, 'E', &NULL_AF_INFO), relation, af_info); }

static long asfreq_MtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, af_info); }

static long asfreq_MtoB(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_MtoD(fromDate, relation, &NULL_AF_INFO),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static long asfreq_MtoH(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_MtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static long asfreq_MtoS(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_MtoD(fromDate, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM QUARTERLY ***************

static void QtoD_ym(long fromDate, long *y, long *m, asfreq_info *af_info) {

    *y = (fromDate - 1) / 4 + 1;
    *m = (fromDate + 4) * 3 - 12 * (*y) - 2;

    if (af_info->from_q_year_end != 12) {
        *m += af_info->from_q_year_end;
        if (*m > 12) { *m -= 12; }
        else { *y -= 1; }
    }
}

static long asfreq_QtoD(long fromDate, char relation, asfreq_info *af_info) {

    long y, m, absdate;

    if (relation == 'S') {
        QtoD_ym(fromDate, &y, &m, af_info);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate;
    } else {
        QtoD_ym(fromDate+1, &y, &m, af_info);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate - 1;
    }
}

static long asfreq_QtoQ(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_QtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_QtoA(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_QtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_QtoM(long fromDate, char relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_QtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }

static long asfreq_QtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_QtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_QtoB(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_QtoD(fromDate, relation, af_info),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}


static long asfreq_QtoH(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_QtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_QtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_QtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_QtoS(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_QtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }


//************ FROM ANNUAL ***************

static long asfreq_AtoD(long fromDate, char relation, asfreq_info *af_info) {
    long absdate, year, final_adj;
    int month = (af_info->from_a_year_end) % 12;

    if (month == 0) { month = 1; }
    else { month += 1; }

    if (relation == 'S') {
        if (af_info->from_a_year_end == 12) {year = fromDate;}
        else {year = fromDate - 1;}
        final_adj = 0;
    } else {
        if (af_info->from_a_year_end == 12) {year = fromDate+1;}
        else {year = fromDate;}
        final_adj = -1;
    }
    absdate = absdate_from_ymd(year, month, 1);
    if (absdate  == INT_ERR_CODE) return INT_ERR_CODE;
    return absdate + final_adj;
}

static long asfreq_AtoA(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_AtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_AtoQ(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_AtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_AtoM(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_AtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_AtoW(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_AtoD(fromDate, relation, af_info), relation, af_info); }

static long asfreq_AtoB(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, asfreq_AtoD(fromDate, relation, af_info),
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static long asfreq_AtoH(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_AtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_AtoT(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_AtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }
static long asfreq_AtoS(long fromDate, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_AtoD(fromDate, relation, af_info), relation, &NULL_AF_INFO); }

static long nofunc(long fromDate, char relation, asfreq_info *af_info) { return -1; }

// end of frequency specific conversion routines

// return a pointer to appropriate conversion function
long (*get_asfreq_func(int fromFreq, int toFreq, int forConvert))(long, char, asfreq_info*) {

    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    if (fromGroup == FR_UND) { fromGroup = FR_DAY; }

    switch(fromGroup)
    {
        case FR_ANN:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_AtoA;
                case FR_QTR: return &asfreq_AtoQ;
                case FR_MTH: return &asfreq_AtoM;
                case FR_WK: return &asfreq_AtoW;
                case FR_BUS: return &asfreq_AtoB;
                case FR_DAY: return &asfreq_AtoD;
                case FR_HR: return &asfreq_AtoH;
                case FR_MIN: return &asfreq_AtoT;
                case FR_SEC: return &asfreq_AtoS;
                default: return &nofunc;
            }

        case FR_QTR:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_QtoA;
                case FR_QTR: return &asfreq_QtoQ;
                case FR_MTH: return &asfreq_QtoM;
                case FR_WK: return &asfreq_QtoW;
                case FR_BUS: return &asfreq_QtoB;
                case FR_DAY: return &asfreq_QtoD;
                case FR_HR: return &asfreq_QtoH;
                case FR_MIN: return &asfreq_QtoT;
                case FR_SEC: return &asfreq_QtoS;
                default: return &nofunc;
            }

        case FR_MTH:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_MtoA;
                case FR_QTR: return &asfreq_MtoQ;
                case FR_WK: return &asfreq_MtoW;
                case FR_BUS: return &asfreq_MtoB;
                case FR_DAY: return &asfreq_MtoD;
                case FR_HR: return &asfreq_MtoH;
                case FR_MIN: return &asfreq_MtoT;
                case FR_SEC: return &asfreq_MtoS;
                default: return &nofunc;
            }

        case FR_WK:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_WtoA;
                case FR_QTR: return &asfreq_WtoQ;
                case FR_MTH: return &asfreq_WtoM;
                case FR_WK: return &asfreq_WtoW;
                case FR_BUS: return &asfreq_WtoB;
                case FR_DAY: return &asfreq_WtoD;
                case FR_HR: return &asfreq_WtoH;
                case FR_MIN: return &asfreq_WtoT;
                case FR_SEC: return &asfreq_WtoS;
                default: return &nofunc;
            }

        case FR_BUS:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_BtoA;
                case FR_QTR: return &asfreq_BtoQ;
                case FR_MTH: return &asfreq_BtoM;
                case FR_WK: return &asfreq_BtoW;
                case FR_DAY: return &asfreq_BtoD;
                case FR_HR: return &asfreq_BtoH;
                case FR_MIN: return &asfreq_BtoT;
                case FR_SEC: return &asfreq_BtoS;
                default: return &nofunc;
            }

        case FR_DAY:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_DtoA;
                case FR_QTR: return &asfreq_DtoQ;
                case FR_MTH: return &asfreq_DtoM;
                case FR_WK: return &asfreq_DtoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_DtoB_forConvert; }
                    else            { return &asfreq_DtoB; }
                case FR_DAY: return &asfreq_DtoD;
                case FR_HR: return &asfreq_DtoH;
                case FR_MIN: return &asfreq_DtoT;
                case FR_SEC: return &asfreq_DtoS;
                default: return &nofunc;
            }

        case FR_HR:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_HtoA;
                case FR_QTR: return &asfreq_HtoQ;
                case FR_MTH: return &asfreq_HtoM;
                case FR_WK: return &asfreq_HtoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_HtoB_forConvert; }
                    else            { return &asfreq_HtoB; }
                case FR_DAY: return &asfreq_HtoD;
                case FR_MIN: return &asfreq_HtoT;
                case FR_SEC: return &asfreq_HtoS;
                default: return &nofunc;
            }

        case FR_MIN:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_TtoA;
                case FR_QTR: return &asfreq_TtoQ;
                case FR_MTH: return &asfreq_TtoM;
                case FR_WK: return &asfreq_TtoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_TtoB_forConvert; }
                    else            { return &asfreq_TtoB; }
                case FR_DAY: return &asfreq_TtoD;
                case FR_HR: return &asfreq_TtoH;
                case FR_SEC: return &asfreq_TtoS;
                default: return &nofunc;
            }

        case FR_SEC:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_StoA;
                case FR_QTR: return &asfreq_StoQ;
                case FR_MTH: return &asfreq_StoM;
                case FR_WK: return &asfreq_StoW;
                case FR_BUS:
                    if (forConvert) { return &asfreq_StoB_forConvert; }
                    else            { return &asfreq_StoB; }
                case FR_DAY: return &asfreq_StoD;
                case FR_HR: return &asfreq_StoH;
                case FR_MIN: return &asfreq_StoT;
                default: return &nofunc;
            }
        default: return &nofunc;
    }
}

static int calc_a_year_end(int freq, int group) {
    int result = (freq - group) % 12;
    if (result == 0) {return 12;}
    else {return result;}
}

static int calc_week_end(int freq, int group) {
    return freq - group;
}

void get_asfreq_info(int fromFreq, int toFreq, asfreq_info *af_info) {

    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    switch(fromGroup)
    {
        case FR_WK: {
            af_info->from_week_end = calc_week_end(fromFreq, fromGroup);
        } break;
        case FR_ANN: {
            af_info->from_a_year_end = calc_a_year_end(fromFreq, fromGroup);
        } break;
        case FR_QTR: {
            af_info->from_q_year_end = calc_a_year_end(fromFreq, fromGroup);
        } break;

    }

    switch(toGroup)
    {
        case FR_WK: {
            af_info->to_week_end = calc_week_end(toFreq, toGroup);
        } break;
        case FR_ANN: {
            af_info->to_a_year_end = calc_a_year_end(toFreq, toGroup);
        } break;
        case FR_QTR: {
            af_info->to_q_year_end = calc_a_year_end(toFreq, toGroup);
        } break;
    }

}

static double getAbsTime(int freq, long dailyDate, long originalDate) {

    long startOfDay, periodsPerDay;

    switch(freq)
    {
        case FR_HR:
            periodsPerDay = 24;
            break;
        case FR_MIN:
            periodsPerDay = 24*60;
            break;
        case FR_SEC:
            periodsPerDay = 24*60*60;
            break;
        default:
            return 24*60*60 - 1;
    }

    startOfDay = asfreq_DtoHIGHFREQ(dailyDate, 'S', periodsPerDay);
    return (24*60*60)*((double)(originalDate - startOfDay))/((double)periodsPerDay);
}

/************************************************************
** Date type definition
************************************************************/

typedef struct {
    PyObject_HEAD
    int freq; /* frequency of date */
    long value; /* integer representation of date */
    PyObject* cached_vals;
} DateObject;

/* Forward declarations */
static PyTypeObject DateType;
#define DateObject_Check(op) PyObject_TypeCheck(op, &DateType)

static void
DateObject_dealloc(DateObject* self) {
    Py_XDECREF(self->cached_vals);
    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *freq_dict, *freq_dict_rev, *freq_constants;

#define DICT_SETINT_STRKEY(dict, key, val) \
    {PyObject *pyval = PyInt_FromLong(val); \
     PyDict_SetItemString(dict, key, pyval); \
     Py_DECREF(pyval); }

#define ADD_FREQ_CONSTANT(const_name, val) \
    DICT_SETINT_STRKEY(freq_constants, const_name, val)

#define INIT_FREQ(const_name, key, aliases) \
    {PyObject *pykey = PyInt_FromLong(key); \
     PyDict_SetItem(freq_dict, pykey, aliases); \
     PyDict_SetItemString(freq_constants, const_name, pykey); \
     Py_DECREF(pykey); \
     Py_DECREF(aliases); }


static int init_freq_group(int num_items, int num_roots, int base_const,
                            char item_abbrevs[][2][10], char group_prefixes[][15],
                            char item_const_names[][15]) {

    int i;

    for (i = 0; i < num_items; i++) {

        PyObject *aliases;
        int j, size, k;

        if (i == 0) { k = 3; } else { k = 2; }

        size = num_roots * k;

        aliases = PyTuple_New(size);

        for (j = 0; j < num_roots; j++) {
            PyObject *alias_v1, *alias_v2;
            char *root, *alt;

            if ((root = PyArray_malloc((30) * sizeof(char))) == NULL) return INT_ERR_CODE;
            if ((alt = PyArray_malloc((30) * sizeof(char))) == NULL) return INT_ERR_CODE;

            strcpy(root, group_prefixes[j]);
            strcpy(alt, group_prefixes[j]);

            if (i == 0) {
                PyObject *alias = PyString_FromString(root);
                PyTuple_SET_ITEM(aliases, j*k + 2, alias);
            }

            strcat(root, "-");
            strcat(root, item_abbrevs[i][0]);
            strcat(alt, "-");
            strcat(alt, item_abbrevs[i][1]);

            alias_v1 = PyString_FromString(root);
            alias_v2 = PyString_FromString(alt);

            free(root);
            free(alt);

            PyTuple_SET_ITEM(aliases, j*k, alias_v1);
            PyTuple_SET_ITEM(aliases, j*k + 1, alias_v2);
        }

        INIT_FREQ(item_const_names[i], base_const+i, aliases);
    }

    return 0;
}

/* take a dictionary with integer keys and tuples of strings for values,
   and populate a dictionary with all the strings as keys and integers
   for values */
static int reverse_dict(PyObject *source, PyObject *dest) {

    PyObject *key, *value;

    Py_ssize_t pos = 0;

    while (PyDict_Next(source, &pos, &key, &value)) {
        PyObject *tuple_iter;
        PyObject *item;

        if((tuple_iter = PyObject_GetIter(value)) == NULL) return INT_ERR_CODE;

        while ((item = PyIter_Next(tuple_iter)) != NULL) {
            PyDict_SetItem(dest, item, key);
            Py_DECREF(item);
        }
        Py_DECREF(tuple_iter);
    }
    return 0;
}

static int build_freq_dict(void) {

    char ANN_prefixes[8][15] = { "A", "Y", "ANN", "ANNUAL", "ANNUALLY",
                                 "YR", "YEAR", "YEARLY" };

    char QTRE_prefixes[8][15] = { "Q", "QTR", "QUARTER", "QUARTERLY", "Q-E",
                                  "QTR-E", "QUARTER-E", "QUARTERLY-E"};
    char QTRS_prefixes[4][15] = { "Q-S", "QTR-S", "QUARTER-S", "QUARTERLY-S" };

    char WK_prefixes[4][15] =  { "W", "WK", "WEEK", "WEEKLY" };

    /* Note: order of this array must match up with how the Annual
       frequency constants are lined up */
    char month_names[12][2][10] = {
        { "DEC", "DECEMBER" },
        { "JAN", "JANUARY" },
        { "FEB", "FEBRUARY" },
        { "MAR", "MARCH" },
        { "APR", "APRIL" },
        { "MAY", "MAY" },
        { "JUN", "JUNE" },
        { "JUL", "JULY" },
        { "AUG", "AUGUST" },
        { "SEP", "SEPTEMBER" },
        { "OCT", "OCTOBER" },
        { "NOV", "NOVEMBER" }};

    char day_names[7][2][10] = {
        { "SUN", "SUNDAY" },
        { "MON", "MONDAY" },
        { "TUE", "TUESDAY" },
        { "WED", "WEDNESDAY" },
        { "THU", "THURSDAY" },
        { "FRI", "FRIDAY" },
        { "SAT", "SATURDAY" }};

    char ANN_const_names[12][15] = {
        "FR_ANNDEC",
        "FR_ANNJAN",
        "FR_ANNFEB",
        "FR_ANNMAR",
        "FR_ANNAPR",
        "FR_ANNMAY",
        "FR_ANNJUN",
        "FR_ANNJUL",
        "FR_ANNAUG",
        "FR_ANNSEP",
        "FR_ANNOCT",
        "FR_ANNNOV"};

    char QTRE_const_names[12][15] = {
        "FR_QTREDEC",
        "FR_QTREJAN",
        "FR_QTREFEB",
        "FR_QTREMAR",
        "FR_QTREAPR",
        "FR_QTREMAY",
        "FR_QTREJUN",
        "FR_QTREJUL",
        "FR_QTREAUG",
        "FR_QTRESEP",
        "FR_QTREOCT",
        "FR_QTRENOV"};

    char QTRS_const_names[12][15] = {
        "FR_QTRSDEC",
        "FR_QTRSJAN",
        "FR_QTRSFEB",
        "FR_QTRSMAR",
        "FR_QTRSAPR",
        "FR_QTRSMAY",
        "FR_QTRSJUN",
        "FR_QTRSJUL",
        "FR_QTRSAUG",
        "FR_QTRSSEP",
        "FR_QTRSOCT",
        "FR_QTRSNOV"};

    char WK_const_names[7][15] = {
        "FR_WKSUN",
        "FR_WKMON",
        "FR_WKTUE",
        "FR_WKWED",
        "FR_WKTHU",
        "FR_WKFRI",
        "FR_WKSAT"};

    PyObject *aliases;

    freq_dict = PyDict_New();
    freq_dict_rev = PyDict_New();
    freq_constants = PyDict_New();

    aliases = Py_BuildValue("(ssss)", "M", "MTH", "MONTH", "MONTHLY");
    INIT_FREQ("FR_MTH", FR_MTH, aliases);

    aliases = Py_BuildValue("(ssss)", "B", "BUS", "BUSINESS", "BUSINESSLY");
    INIT_FREQ("FR_BUS", FR_BUS, aliases);

    aliases = Py_BuildValue("(ssss)", "D", "DAY", "DLY", "DAILY");
    INIT_FREQ("FR_DAY", FR_DAY, aliases);

    aliases = Py_BuildValue("(sssss)", "H", "HR", "HOUR", "HRLY", "HOURLY");
    INIT_FREQ("FR_HR", FR_HR, aliases);

    aliases = Py_BuildValue("(ssss)", "T", "MIN", "MINUTE", "MINUTELY");
    INIT_FREQ("FR_MIN", FR_MIN, aliases);

    aliases = Py_BuildValue("(ssss)", "S", "SEC", "SECOND", "SECONDLY");
    INIT_FREQ("FR_SEC", FR_SEC, aliases);

    aliases = Py_BuildValue("(ssss)", "U", "UND", "UNDEF", "UNDEFINED");
    INIT_FREQ("FR_UND", FR_UND, aliases);

    ADD_FREQ_CONSTANT("FR_ANN", FR_ANN);

    if(init_freq_group(12, 8, FR_ANN,
        month_names, ANN_prefixes, ANN_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    ADD_FREQ_CONSTANT("FR_QTR", FR_QTR);

    if(init_freq_group(12, 8, FR_QTREDEC,
        month_names, QTRE_prefixes, QTRE_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    if(init_freq_group(12, 4, FR_QTRSDEC,
        month_names, QTRS_prefixes, QTRS_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    ADD_FREQ_CONSTANT("FR_WK", FR_WK);

    if(init_freq_group(7, 4, FR_WK,
                    day_names, WK_prefixes, WK_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    if(reverse_dict(freq_dict, freq_dict_rev) == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    return 0;
}


/* take user specified frequency and convert to int representation
   of the frequency */
int check_freq(PyObject *freq_spec) {

    if (PyInt_Check(freq_spec)) {
        return (int)PyInt_AsLong(freq_spec);
    } else if (PyString_Check(freq_spec)) {
        char *freq_str, *freq_str_uc;
        PyObject *freq_val;

        freq_str = PyString_AsString(freq_spec);
        if((freq_str_uc = str_uppercase(freq_str)) == NULL) {return INT_ERR_CODE;}

        freq_val = PyDict_GetItemString(freq_dict_rev, freq_str_uc);

        free(freq_str_uc);

        if (freq_val == NULL) {
            PyErr_SetString(PyExc_ValueError, "invalid frequency specification");
            return INT_ERR_CODE;
        } else {
            int ret_val = (int)PyInt_AsLong(freq_val);
            return ret_val;
        }
    } else if (freq_spec == Py_None) {
        return FR_UND;
    } else {
        int retval = (int)PyInt_AsLong(freq_spec);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError, "invalid frequency specification");
            return INT_ERR_CODE;
        } else { return retval; }
    }

}

static PyObject *
DateObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {

    DateObject *self;

    self = (DateObject*)type->tp_alloc(type, 0);
    if (self != NULL) {
        // initialize attributes that need initializing in here
        self->freq = FR_UND;
        self->value = -1;
    }

    return (PyObject *)self;
}

/* for use in C code */
static DateObject *
DateObject_New(void) {
    PyObject *dummy;
    return (DateObject*)DateObject_new(&DateType, dummy, dummy);
}

#define INIT_ERR(errortype, errmsg) PyErr_SetString(errortype,errmsg);return -1

static int
DateObject_init(DateObject *self, PyObject *args, PyObject *kwds) {

    PyObject *freq=NULL, *value=NULL, *datetime=NULL, *string=NULL;
    char *INSUFFICIENT_MSG = "insufficient parameters to initialize Date";

    int def_info=INT_ERR_CODE;

    int year=def_info, month=def_info, day=def_info, quarter=def_info,
        hour=def_info, minute=def_info, second=def_info;

    int free_dt=0;

    static char *kwlist[] = {"freq", "value", "string",
                             "year", "month", "day", "quarter",
                             "hour", "minute", "second",
                             "datetime", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|OOiiiiiiiO", kwlist,
                                      &freq, &value, &string,
                                      &year, &month, &day, &quarter,
                                      &hour, &minute, &second,
                                      &datetime)) {
        return -1;
    }

    if (PyObject_HasAttrString(freq, "freq")) {
        PyObject *freq_attr = PyObject_GetAttrString(freq, "freq");
        self->freq = PyInt_AS_LONG(freq_attr);
        Py_DECREF(freq_attr);
    } else {
        if((self->freq = check_freq(freq)) == INT_ERR_CODE) return -1;
    }

    if ((value && PyString_Check(value)) || string) {

        PyObject *string_arg = PyTuple_New(1);
        int freq_group = get_freq_group(self->freq);

        free_dt = 1;

        if (!string) {
            string = value;
        }

        PyTuple_SET_ITEM(string_arg, 0, string);
        Py_INCREF(string);

        if (freq_group == FR_HR ||
            freq_group == FR_MIN ||
            freq_group == FR_SEC)
             { datetime = PyEval_CallObject(DateTimeFromString, string_arg); }
        else { datetime = PyEval_CallObject(DateFromString, string_arg); }

        Py_DECREF(string_arg);

        value = NULL;
    }

    if (value && (PyDateTime_Check(value) || PyDate_Check(value))) {
        if (!datetime) {
            datetime = value;
        }
        value = NULL;
    } // datetime = (datetime||value), value = NULL


    if (value) {
        self->value = PyInt_AsLong(value);
    } else {

        int freq_group = get_freq_group(self->freq);

        if (datetime) {
            if (PyDateTime_Check(datetime) || PyDate_Check(datetime)) {
                year=PyDateTime_GET_YEAR(datetime);
                month=PyDateTime_GET_MONTH(datetime);
                quarter=((month-1)/3)+1;
                day=PyDateTime_GET_DAY(datetime);
                hour=PyDateTime_DATE_GET_HOUR(datetime);
                minute=PyDateTime_DATE_GET_MINUTE(datetime);
                second=PyDateTime_DATE_GET_SECOND(datetime);
            } else {
                PyObject *err_msg, *_type;
                _type = PyObject_Type(datetime);
                err_msg = PyString_FromString("Expected datetime object, received: ");
                PyString_ConcatAndDel(&err_msg, PyObject_Str(_type));
                PyErr_SetString(PyExc_TypeError, PyString_AsString(err_msg));
                Py_DECREF(_type);
                Py_DECREF(err_msg);
                return -1;
            }
        }

        if (!datetime) {

            // First, some basic checks.....
            if (year == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
            if (self->freq == FR_BUS ||
               self->freq == FR_DAY ||
               self->freq == FR_WK ||
               self->freq == FR_UND) {
                if (month == def_info || day == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }

                // if FR_BUS, check for week day

            } else if (self->freq == FR_MTH) {
                if (month == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            } else if (freq_group == FR_QTR) {
                if (quarter == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            } else if (self->freq == FR_SEC) {
                if (month == def_info ||
                    day == def_info ||
                    second == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
                if (hour == def_info) {
                    hour = second/3600;
                    minute = (second % 3600)/60;
                    second = second % 60;
                } else if (minute == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            } else if (self->freq == FR_MIN) {
                if (month == def_info ||
                    day == def_info ||
                    minute == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
                if (hour == def_info) {
                    hour = minute/60;
                    minute = minute % 60;
                }
            } else if (self->freq == FR_HR) {
                if (month == def_info ||
                    day == def_info ||
                    hour == def_info) {
                    INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
                }
            }

        }

        if (self->freq == FR_SEC) {
            long absdays, delta;
            absdays = absdate_from_ymd(year, month, day);
            delta = (absdays - HIGHFREQ_ORIG);
            self->value = (int)(delta*86400 + hour*3600 + minute*60 + second + 1);
        } else if (self->freq == FR_MIN) {
            long absdays, delta;
            absdays = absdate_from_ymd(year, month, day);
            delta = (absdays - HIGHFREQ_ORIG);
            self->value = (int)(delta*1440 + hour*60 + minute + 1);
        } else if (self->freq == FR_HR) {
            long absdays, delta;
            if((absdays = absdate_from_ymd(year, month, day)) == INT_ERR_CODE) return -1;
            delta = (absdays - HIGHFREQ_ORIG);
            self->value = (int)(delta*24 + hour + 1);
        } else if (self->freq == FR_DAY) {
            if((self->value = (int)absdate_from_ymd(year, month, day)) == INT_ERR_CODE) return -1;
        } else if (self->freq == FR_UND) {
            if((self->value = (int)absdate_from_ymd(year, month, day)) == INT_ERR_CODE) return -1;
        } else if (self->freq == FR_BUS) {
            long weeks, days;
            if((days = absdate_from_ymd(year, month, day)) == INT_ERR_CODE) return -1;
            weeks = days/7;
            self->value = (int)(days - weeks*2);
        } else if (freq_group == FR_WK) {
            int adj_ordinal, ordinal, day_adj;
            if((ordinal = (int)absdate_from_ymd(year, month, day)) == INT_ERR_CODE) return -1;
            day_adj = (7 - (self->freq - FR_WK)) % 7;
            adj_ordinal = ordinal + ((7 - day_adj) - ordinal % 7) % 7;
            self->value = adj_ordinal/7;
        } else if (self->freq == FR_MTH) {
            self->value = (year-1)*12 + month;
        } else if (freq_group == FR_QTR) {
            if ((self->freq - freq_group) > 12) {
                // quarterly frequency with year determined by ending period
                self->value = year*4 + quarter;
            } else {
                /* quarterly frequency with year determined by ending period
                   or has December year end*/
                self->value = (year-1)*4 + quarter;
            }
        } else if (freq_group == FR_ANN) {
            self->value = year;
        }

    }

    if (free_dt) { Py_DECREF(datetime); }

    return 0;
}

static PyMemberDef DateObject_members[] = {
    {"freq", T_INT, offsetof(DateObject, freq), 0,
     "frequency"},
    {"value", T_INT, offsetof(DateObject, value), 0,
     "integer representation of the Date"},
    {NULL}  /* Sentinel */
};

static char DateObject_toordinal_doc[] =
"Returns the proleptic Gregorian ordinal of the date, as an integer.\n"
"This corresponds to the number of days since Jan., 1st, 1AD.\n\n"
"When the instance has a frequency less than daily, the proleptic date \n"
"is calculated for the last day of the period.\n\n"
"   >>> ts.Date('D', '2001-01-01').toordinal()\n"
"   730486\n"
"   >>> ts.Date('H', '2001-01-01 18:00').toordinal()\n"
"   730486\n"
"   >>> ts.Date('M', '2001-01-01').toordinal()\n"
"   730516\n"
"   >>> # Note that 730516 = 730486 + 31 - 1\n"
"   >>> ts.Date('Y', '2001-01-01').toordinal()\n"
"   730850\n"
"   >>> # Note that 730850 = 730486 + 365 - 1\n";

static PyObject *
DateObject_toordinal(DateObject* self)
{
    if (self->freq == FR_DAY) {
        return PyInt_FromLong(self->value);
    } else {
        long (*toDaily)(long, char, asfreq_info*) = NULL;
        asfreq_info af_info;

        toDaily = get_asfreq_func(self->freq, FR_DAY, 0);
        get_asfreq_info(self->freq, FR_DAY, &af_info);

        return PyInt_FromLong(toDaily(self->value, 'E', &af_info));
    }
}

static char DateObject_asfreq_doc[] =
"   asfreq(freq, relation='END')\n"
"\n"
"   Returns a :class:`Date` object converted to a specified frequency.\n"
"\n"
"   :Parameters:\n"
"\n"
"      **freq** : {string, integer}\n"
"         Frequency to convert the instance to. Accepts any valid frequency\n"
"         specification (string or integer).\n"
"\n"
"      **relation** : {'END', 'START'} (optional)\n"
"         Applies only when converting a :class:`Date` to a higher frequency,\n"
"         or when converting a weekend Date to a business frequency Date.\n"
"         Valid values are 'START' and 'END'.\n"
"         For example, when converting a monthly :class:`Date` to the daily\n"
"         frequency, ``relation='START'`` gives the first day of the month\n"
"         while ``relation='END'`` gives the last day of the month.\n"
"\n"
"   .. warning::\n"
"\n"
"      Some information will be lost when a :class:`Date` is converted to \n"
"      a lower frequency and then back to the original one.\n"
"      For example, if a daily :class:`Date` is converted to monthly and \n"
"      then back to a daily one, the :attr:`day` information is lost::\n"
"\n"
"         >>> D = ts.Date('D', year=2007, month=12, day=15)\n"
"         >>> D.asfreq('M')\n"
"         <M: Dec-2007>\n"
"         >>> D.asfreq('M').asfreq('D', relation='START')\n"
"         <D: 01-Dec-2007>\n"
"         >>> D.asfreq('M').asfreq('D', relation=\"END\")\n"
"         <D: 31-Dec-2007>\n"
"\n";

static PyObject *
DateObject_asfreq(DateObject *self, PyObject *args, PyObject *kwds)
{

    PyObject *freq=NULL;
    char *relation_raw=NULL;
    char *relation_uc;
    char relation;
    int invalid_relation=0;
    int toFreq;
    int result_val;
    DateObject *result = DateObject_New();

    static char *kwlist[] = {"freq", "relation", NULL};

    long (*asfreq_func)(long, char, asfreq_info*) = NULL;
    asfreq_info af_info;

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|s", kwlist,
                                      &freq, &relation_raw)) return NULL;

    if(relation_raw) {
        if (strlen(relation_raw) > 0) {
            if((relation_uc = str_uppercase(relation_raw)) == NULL)
            {return PyErr_NoMemory();}

            // 'BEFORE' and 'AFTER' values for this parameter are deprecated
            if (strcmp(relation_uc, "END") == 0 ||
                strcmp(relation_uc, "E") == 0 ||
                strcmp(relation_uc, "START") == 0 ||
                strcmp(relation_uc, "S") == 0 ||
                strcmp(relation_uc, "BEFORE") == 0 ||
                strcmp(relation_uc, "B") == 0 ||
                strcmp(relation_uc, "AFTER") == 0 ||
                strcmp(relation_uc, "A") == 0) {
                 if(relation_uc[0] == 'E' || relation_uc[0] == 'A') { relation = 'E'; }
                 else { relation = 'S'; }

            } else { invalid_relation=1; }

            free(relation_uc);

        } else {
            invalid_relation=1;
        }

        if (invalid_relation) {
            PyErr_SetString(PyExc_ValueError,"Invalid relation specification");
            return NULL;
        }
    } else {
        relation = 'E';
    }

    if ((toFreq = check_freq(freq)) == INT_ERR_CODE) return NULL;

    if (toFreq == self->freq) {
        result->freq = self->freq;
        result->value = self->value;
        return (PyObject*)result;
    }

    get_asfreq_info(self->freq, toFreq, &af_info);
    asfreq_func = get_asfreq_func(self->freq, toFreq, 0);

    result_val = asfreq_func(self->value, relation, &af_info);

    if (result_val == INT_ERR_CODE) return NULL;

    result->freq = toFreq;
    result->value = result_val;

    return (PyObject*)result;

}

static char DateObject_strfmt_doc[] =
"Deprecated alias for strftime method";

static char DateObject_strftime_doc[] =
"\n"
"   Returns the string representation of the :class:`Date`, \n"
"   depending on the selected :keyword:`format`.\n"
"   :keyword:`format` must be a string containing one or several directives.\n"
"   The method recognizes the same directives as the :func:`time.strftime` \n"
"   function of the standard Python distribution, as well as the specific \n"
"   additional directives ``%f``, ``%F``, ``%q``.\n"
"\n"
"   +-----------+--------------------------------+-------+\n"
"   | Directive | Meaning                        | Notes |\n"
"   +===========+================================+=======+\n"
"   | ``%a``    | Locale's abbreviated weekday   |       |\n"
"   |           | name.                          |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%A``    | Locale's full weekday name.    |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%b``    | Locale's abbreviated month     |       |\n"
"   |           | name.                          |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%B``    | Locale's full month name.      |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%c``    | Locale's appropriate date and  |       |\n"
"   |           | time representation.           |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%d``    | Day of the month as a decimal  |       |\n"
"   |           | number [01,31].                |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%f``    | 'Fiscal' year without a        | \(1)  |\n"
"   |           | century  as a decimal number   |       |\n"
"   |           | [00,99]                        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%F``    | 'Fiscal' year with a century   | \(2)  |\n"
"   |           | as a decimal number            |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%H``    | Hour (24-hour clock) as a      |       |\n"
"   |           | decimal number [00,23].        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%I``    | Hour (12-hour clock) as a      |       |\n"
"   |           | decimal number [01,12].        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%j``    | Day of the year as a decimal   |       |\n"
"   |           | number [001,366].              |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%m``    | Month as a decimal number      |       |\n"
"   |           | [01,12].                       |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%M``    | Minute as a decimal number     |       |\n"
"   |           | [00,59].                       |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%p``    | Locale's equivalent of either  | \(3)  |\n"
"   |           | AM or PM.                      |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%q``    | Quarter as a decimal number    |       |\n"
"   |           | [01,04]                        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%S``    | Second as a decimal number     | \(4)  |\n"
"   |           | [00,61].                       |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%U``    | Week number of the year        | \(5)  |\n"
"   |           | (Sunday as the first day of    |       |\n"
"   |           | the week) as a decimal number  |       |\n"
"   |           | [00,53].  All days in a new    |       |\n"
"   |           | year preceding the first       |       |\n"
"   |           | Sunday are considered to be in |       |\n"
"   |           | week 0.                        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%w``    | Weekday as a decimal number    |       |\n"
"   |           | [0(Sunday),6].                 |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%W``    | Week number of the year        | \(5)  |\n"
"   |           | (Monday as the first day of    |       |\n"
"   |           | the week) as a decimal number  |       |\n"
"   |           | [00,53].  All days in a new    |       |\n"
"   |           | year preceding the first       |       |\n"
"   |           | Monday are considered to be in |       |\n"
"   |           | week 0.                        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%x``    | Locale's appropriate date      |       |\n"
"   |           | representation.                |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%X``    | Locale's appropriate time      |       |\n"
"   |           | representation.                |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%y``    | Year without century as a      |       |\n"
"   |           | decimal number [00,99].        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%Y``    | Year with century as a decimal |       |\n"
"   |           | number.                        |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%Z``    | Time zone name (no characters  |       |\n"
"   |           | if no time zone exists).       |       |\n"
"   +-----------+--------------------------------+-------+\n"
"   | ``%%``    | A literal ``'%'`` character.   |       |\n"
"   +-----------+--------------------------------+-------+\n"
"\n"
"   .. note::\n"
"\n"
"      (1)\n"
"         The ``%f`` directive is the same as ``%y`` if the frequency is \n"
"         not quarterly.\n"
"         Otherwise, it corresponds to the 'fiscal' year, as defined by \n"
"         the :attr:`qyear` attribute.\n"
"\n"
"      (2)\n"
"         The ``%F`` directive is the same as ``%Y`` if the frequency is \n"
"         not quarterly.\n"
"         Otherwise, it corresponds to the 'fiscal' year, as defined by \n"
"         the :attr:`qyear` attribute.\n"
"\n"
"      (3)\n"
"         The ``%p`` directive only affects the output hour field \n"
"         if the ``%I`` directive is used to parse the hour.\n"
"\n"
"      (4)\n"
"         The range really is ``0`` to ``61``; this accounts for leap seconds \n"
"         and the (very rare) double leap seconds.\n"
"\n"
"      (5)\n"
"         The ``%U`` and ``%W`` directives are only used in calculations \n"
"         when the day of the week and the year are specified.\n"
"\n"
"\n"
"   .. rubric::  Examples\n"
"\n"
"   >>> a = ts.Date(freq='q-jul', year=2006, quarter=1)\n"
"   >>> a.strftime('%F-Q%q')\n"
"   '2006-Q1'\n"
"   >>> # Output the last month in the quarter of this date\n"
"   >>> a.strftime('%b-%Y')\n"
"   'Oct-2005'\n"
"   >>> \n"
"   >>> a = ts.Date(freq='d', year=2001, month=1, day=1)\n"
"   >>> a.strftime('%d-%b-%Y')\n"
"   '01-Jan-2006'\n"
"   >>> a.strftime('%b. %d, %Y was a %A')\n"
"   'Jan. 01, 2001 was a Monday'\n";
static PyObject *
DateObject_strftime(DateObject *self, PyObject *args)
{

    char *orig_fmt_str, *fmt_str;
    char *result;

    int num_extra_fmts = 3;

    char extra_fmts[3][2][10] = {{"%q", "^`AB`^"},
                                 {"%f", "^`CD`^"},
                                 {"%F", "^`EF`^"}};

    int extra_fmts_found[3] = {0,0,0};
    int extra_fmts_found_one = 0;
    struct tm c_date;
    struct date_info tempDate;
    long absdate;
    double abstime;
    int i, result_len;
    PyObject *py_result;

    long (*toDaily)(long, char, asfreq_info*) = NULL;
    asfreq_info af_info;

    if (!PyArg_ParseTuple(args, "s:strftime(fmt)", &orig_fmt_str)) return NULL;

    toDaily = get_asfreq_func(self->freq, FR_DAY, 0);
    get_asfreq_info(self->freq, FR_DAY, &af_info);

    absdate = toDaily(self->value, 'E', &af_info);
    abstime = getAbsTime(self->freq, absdate, self->value);

    if(dInfoCalc_SetFromAbsDateTime(&tempDate, absdate, abstime,
                                    GREGORIAN_CALENDAR)) return NULL;

    // populate standard C date struct with info from our date_info struct
    c_date.tm_sec = (int)tempDate.second;
    c_date.tm_min = tempDate.minute;
    c_date.tm_hour = tempDate.hour;
    c_date.tm_mday = tempDate.day;
    c_date.tm_mon = tempDate.month - 1;
    c_date.tm_year = tempDate.year - 1900;
    c_date.tm_wday = (tempDate.day_of_week + 1) % 7;
    c_date.tm_yday = tempDate.day_of_year - 1;
    c_date.tm_isdst = -1;

    result_len = strlen(orig_fmt_str) + 50;
    if ((result = PyArray_malloc(result_len * sizeof(char))) == NULL) {return PyErr_NoMemory();}

    fmt_str = orig_fmt_str;

    // replace any special format characters with their place holder
    for(i=0; i < num_extra_fmts; i++) {
        char *special_loc;
        if ((special_loc = strstr(fmt_str,extra_fmts[i][0])) != NULL) {
            char *tmp_str = fmt_str;
            fmt_str = str_replace(fmt_str, extra_fmts[i][0],
                                           extra_fmts[i][1]);
            /* only free the previous loop value if this is not the first
               special format string found */
            if (extra_fmts_found_one) { free(tmp_str); }

            if (fmt_str == NULL) {return NULL;}

            extra_fmts_found[i] = 1;
            extra_fmts_found_one = 1;
        }
    }

    strftime(result, result_len, fmt_str, &c_date);
    if (extra_fmts_found_one) { free(fmt_str); }

    // replace any place holders with the appropriate value
    for(i=0; i < num_extra_fmts; i++) {
        if (extra_fmts_found[i]) {
            char *tmp_str = result;
            char *extra_str;

            if (strcmp(extra_fmts[i][0], "%q") == 0 ||
                strcmp(extra_fmts[i][0], "%f") == 0 ||
                strcmp(extra_fmts[i][0], "%F") == 0) {

                asfreq_info af_info;
                int qtr_freq, year, quarter, year_len;

                if (get_freq_group(self->freq) == FR_QTR) {
                    qtr_freq = self->freq;
                } else { qtr_freq = FR_QTR; }
                get_asfreq_info(FR_DAY, qtr_freq, &af_info);

                if(DtoQ_yq(absdate, &af_info, &year, &quarter) == INT_ERR_CODE)
                { return NULL; }

                if(strcmp(extra_fmts[i][0], "%q") == 0) {
                    if ((extra_str = PyArray_malloc(2 * sizeof(char))) == NULL) {
                        free(tmp_str);
                        return PyErr_NoMemory();
                    }
                    sprintf(extra_str, "%i", quarter);
                } else {
                    if ((qtr_freq % 1000) > 12) { year -= 1; }

                    if (strcmp(extra_fmts[i][0], "%f") == 0) {
                        year_len = 2;
                        year = year % 100;
                    } else { year_len = 4; }

                    if ((extra_str = PyArray_malloc((year_len+1) * sizeof(char))) == NULL) {
                        free(tmp_str);
                        return PyErr_NoMemory();
                    }

                    if (year_len == 2 && year < 10) {
                        sprintf(extra_str, "0%i", year);
                    } else { sprintf(extra_str, "%i", year); }
                }

            } else {
                PyErr_SetString(PyExc_RuntimeError,"Unrecognized format string");
                return NULL;
            }

            result = str_replace(result, extra_fmts[i][1], extra_str);
            free(tmp_str);
            free(extra_str);
            if (result == NULL) { return NULL; }
        }
    }

    py_result = PyString_FromString(result);
    free(result);

    return py_result;
}

static PyObject *
DateObject___str__(DateObject* self)
{

    int freq_group = get_freq_group(self->freq);
    PyObject *string_arg, *retval;

    string_arg = NULL;
    if (freq_group == FR_UND) {
        retval = PyString_FromFormat("%ld", self->value);
        return retval;}
    else if (freq_group == FR_ANN) { string_arg = Py_BuildValue("(s)", "%Y"); }
    else if (freq_group == FR_QTR) { string_arg = Py_BuildValue("(s)", "%FQ%q"); }
    else if (freq_group == FR_MTH) { string_arg = Py_BuildValue("(s)", "%b-%Y"); }
    else if (freq_group == FR_DAY ||
             freq_group == FR_BUS ||
             freq_group == FR_WK) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y"); }
    else if (freq_group == FR_HR) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:00"); }
    else if (freq_group == FR_MIN) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M"); }
    else if (freq_group == FR_SEC) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M:%S"); }

    if (string_arg == NULL) { return NULL; }

    retval = DateObject_strftime(self, string_arg);
    Py_DECREF(string_arg);

    return retval;
}

static PyObject *
DateObject_freqstr(DateObject *self, void *closure) {
    PyObject *key = PyInt_FromLong(self->freq);
    PyObject *freq_aliases = PyDict_GetItem(freq_dict, key);
    PyObject *main_alias = PyTuple_GET_ITEM(freq_aliases, 0);
    Py_DECREF(key);
    Py_INCREF(main_alias);
    return main_alias;
}


static PyObject *
DateObject___repr__(DateObject* self)
{
    PyObject *py_str_rep, *py_freqstr, *py_repr;
    char *str_rep, *freqstr, *repr;
    int repr_len;

    py_str_rep = DateObject___str__(self);
    if (py_str_rep == NULL) { return NULL; }

    py_freqstr = DateObject_freqstr(self, NULL);

    str_rep = PyString_AsString(py_str_rep);
    freqstr = PyString_AsString(py_freqstr);

    repr_len = strlen(str_rep) + strlen(freqstr) + 6;

    if((repr = PyArray_malloc((repr_len + 1) * sizeof(char))) == NULL)
    { return PyErr_NoMemory(); }

    strcpy(repr, "<");
    strcat(repr, freqstr);
    strcat(repr, " : ");
    strcat(repr, str_rep);
    strcat(repr, ">");

    py_repr = PyString_FromString(repr);

    Py_DECREF(py_str_rep);
    Py_DECREF(py_freqstr);

    free(repr);

    return py_repr;
}

/******************************
   These methods seem rather useless. May or may not implement them.
fromordinal(self, ordinal):
    return Date(self.freq, datetime=dt.datetime.fromordinal(ordinal))
tostring(self):
    return str(self)
toobject(self):
    return self
isvalid(self):
    return True
*******************************/


static DateObject *
DateObject_FromFreqAndValue(int freq, int value) {

    DateObject *result = DateObject_New();

    PyObject *args = PyTuple_New(0);
    PyObject *kw = PyDict_New();
    PyObject *py_freq = PyInt_FromLong(freq);
    PyObject *py_value = PyInt_FromLong(value);

    PyDict_SetItemString(kw, "freq", py_freq);
    PyDict_SetItemString(kw, "value", py_value);

    Py_DECREF(py_freq);
    Py_DECREF(py_value);

    DateObject_init(result, args, kw);

    Py_DECREF(args);
    Py_DECREF(kw);

    return result;
}

static PyObject *
DateObject_date_plus_int(PyObject *date, PyObject *pyint) {
    DateObject *dateobj = (DateObject*)date;

    if (!PyInt_Check(pyint) && !PyObject_HasAttrString(pyint, "__int__")) {
        // invalid type for addition

        char *err_str, *type_str;
        PyObject *type_repr, *obj_type;

        obj_type = PyObject_Type(pyint);
        type_repr = PyObject_Repr(obj_type);
        type_str = PyString_AsString(type_repr);

        if ((err_str = PyArray_malloc(255 * sizeof(char))) == NULL) {
            return PyErr_NoMemory();
        }
        sprintf(err_str, "Cannot add Date and %s", type_str);
        Py_DECREF(obj_type);
        Py_DECREF(type_repr);
        PyErr_SetString(PyExc_TypeError, err_str);
        free(err_str);
        return NULL;
    }

    return (PyObject*)DateObject_FromFreqAndValue(
        dateobj->freq, PyInt_AsLong(pyint) + dateobj->value);
}

static PyObject *
DateObject___add__(PyObject *left, PyObject *right)
{
    if (DateObject_Check(left) && DateObject_Check(right)) {
        PyErr_SetString(PyExc_TypeError, "Cannot add Date to Date");
        return NULL;
    } else if (DateObject_Check(left)) {
        return DateObject_date_plus_int(left, right);
    } else {
        return DateObject_date_plus_int(right, left);
    }
}

static PyObject *
DateObject___subtract__(PyObject *left, PyObject *right)
{
    int result;
    DateObject *dleft;
    if (!DateObject_Check(left)) {
        PyErr_SetString(PyExc_ValueError, "Cannot subtract a Date from a non-Date object.");
        return NULL;
    }

    dleft = (DateObject*)left;

    if (DateObject_Check(right)) {
        DateObject *dright = (DateObject*)right;
        if (dleft->freq != dright->freq) {
            PyErr_SetString(PyExc_ValueError, "Cannot subtract Date objects with different frequencies.");
            return NULL;
        }
        result = dleft->value - dright->value;
        return PyInt_FromLong(result);
    } else {
        result = dleft->value - PyInt_AsLong(right);
        return (PyObject*)DateObject_FromFreqAndValue(dleft->freq, result);
    }
}

static int
DateObject___compare__(DateObject * obj1, DateObject * obj2)
{
    if (obj1->freq != obj2->freq) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot compare Date objects with different frequencies.");
        return -1;
    }

    if (obj1->value < obj2->value) return -1;
    if (obj1->value > obj2->value) return 1;
    if (obj1->value == obj2->value) return 0;
    return -1;
}

static long
DateObject___hash__(DateObject *self)
{
    register int freq_group = get_freq_group(self->freq);

    /* within a given frequency, hash values are guaranteed to be unique
       for different dates. For different frequencies, we make a reasonable
       effort to ensure hash values will be unique, but it is not guaranteed */
    if (freq_group == FR_BUS) {
        return self->value + 10000000;
    } else if (freq_group == FR_WK) {
        return self->value + 100000000;
    } else { return self->value; }
}

static PyObject *
DateObject___int__(DateObject *self)
{
    return PyInt_FromLong(self->value);
}

static PyObject *
DateObject___float__(DateObject *self)
{
    return PyFloat_FromDouble((double)(self->value));
}

static PyObject *
DateObject___long__(DateObject *self)
{
    return PyLong_FromLong(self->value);
}


/***************************************************
           ====== Date Properties ======
****************************************************/

// helper function for date property funcs
static int
DateObject_set_date_info(DateObject *self, struct date_info *dinfo) {
    PyObject *daily_obj = DateObject_toordinal(self);
    long absdate = PyInt_AsLong(daily_obj);

    Py_DECREF(daily_obj);

    if(dInfoCalc_SetFromAbsDate(dinfo, absdate,
                                GREGORIAN_CALENDAR)) return -1;

    return 0;
}

// helper function for date property funcs
static int
DateObject_set_date_info_wtime(DateObject *self, struct date_info *dinfo) {
    PyObject *daily_obj = DateObject_toordinal(self);
    long absdate = PyInt_AsLong(daily_obj);
    double abstime;

    Py_DECREF(daily_obj);

    abstime = getAbsTime(self->freq, absdate, self->value);

    if(dInfoCalc_SetFromAbsDateTime(dinfo, absdate, abstime,
                                    GREGORIAN_CALENDAR)) return -1;

    return 0;
}

static PyObject *
DateObject_year(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.year);
}

static int _DateObject_quarter_year(DateObject *self, int *year, int *quarter) {

    PyObject *daily_obj;
    long absdate;

    asfreq_info af_info;
    int qtr_freq;

    daily_obj = DateObject_toordinal(self);
    absdate = PyInt_AsLong(daily_obj);
    Py_DECREF(daily_obj);

    if (get_freq_group(self->freq) == FR_QTR) {
        qtr_freq = self->freq;
    } else { qtr_freq = FR_QTR; }
    get_asfreq_info(FR_DAY, qtr_freq, &af_info);

    if(DtoQ_yq(absdate, &af_info, year, quarter) == INT_ERR_CODE)
    { return INT_ERR_CODE; }

    if ((qtr_freq % 1000) > 12) { *year -= 1; }

    return 0;
}

static PyObject *
DateObject_qyear(DateObject *self, void *closure) {
    int year, quarter;
    if(_DateObject_quarter_year(self,
            &year, &quarter) == INT_ERR_CODE) { return NULL; }
    return PyInt_FromLong(year);
}

static PyObject *
DateObject_quarter(DateObject *self, void *closure) {
    int year, quarter;
    if(_DateObject_quarter_year(self,
            &year, &quarter) == INT_ERR_CODE) { return NULL; }
    return PyInt_FromLong(quarter);
}

static PyObject *
DateObject_month(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.month);
}

static PyObject *
DateObject_day(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day);
}

static PyObject *
DateObject_weekday(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day_of_week);
}

static PyObject *
DateObject_day_of_week(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day_of_week);
}

static PyObject *
DateObject_day_of_year(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day_of_year);
}

static PyObject *
DateObject_week(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dInfoCalc_ISOWeek(&dinfo));
}

static PyObject *
DateObject_hour(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.hour);
}

static PyObject *
DateObject_minute(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.minute);
}

static PyObject *
DateObject_second(DateObject *self, void *closure) {
    struct date_info dinfo;
    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong((int)dinfo.second);
}

static PyObject *
DateObject_datetime(DateObject *self, void *closure) {
    PyObject *datetime;
    int hour=0, minute=0, second=0;
    int freq_group;
    struct date_info dinfo;

    if(DateObject_set_date_info_wtime(self, &dinfo) == -1) return NULL;
    freq_group = get_freq_group(self->freq);

    switch(freq_group) {
        case FR_HR:
            hour = dinfo.hour;
            break;
        case FR_MIN:
            hour = dinfo.hour;
            minute = dinfo.minute;
            break;
        case FR_SEC:
            hour = dinfo.hour;
            minute = dinfo.minute;
            second = (int)dinfo.second;
            break;
    }

    datetime = PyDateTime_FromDateAndTime(
                dinfo.year, dinfo.month, dinfo.day, hour, minute, second, 0);
    return datetime;
}

static int
DateObject_ReadOnlyErr(DateObject *self, PyObject *value, void *closure) {
   PyErr_SetString(PyExc_AttributeError, "Cannot set read-only property");
   return -1;
}

static PyGetSetDef DateObject_getseters[] = {
    {"year", (getter)DateObject_year, (setter)DateObject_ReadOnlyErr,
            "Returns the year.", NULL},
    {"qyear", (getter)DateObject_qyear, (setter)DateObject_ReadOnlyErr,
            "For quarterly frequency dates, returns the year corresponding to the\n"
            "year end (start) month. When using QTR or QTR-E based quarterly\n"
            "frequencies, this is the fiscal year in a financial context.\n\n"
            "For non-quarterly dates, this simply returns the year of the date.",
            NULL},
    {"quarter", (getter)DateObject_quarter, (setter)DateObject_ReadOnlyErr,
            "Returns the quarter.", NULL},
    {"month", (getter)DateObject_month, (setter)DateObject_ReadOnlyErr,
            "Returns the month.", NULL},
    {"week", (getter)DateObject_week, (setter)DateObject_ReadOnlyErr,
            "Returns the week.", NULL},
    {"day", (getter)DateObject_day, (setter)DateObject_ReadOnlyErr,
            "Returns the day of month.", NULL},
    {"weekday", (getter)DateObject_weekday, (setter)DateObject_ReadOnlyErr,
            "Returns the day of week.", NULL},
    // deprecated alias for weekday property
    {"day_of_week", (getter)DateObject_weekday, (setter)DateObject_ReadOnlyErr,
            "Returns the day of week.", NULL},
    {"day_of_year", (getter)DateObject_day_of_year, (setter)DateObject_ReadOnlyErr,
            "Returns the day of year.", NULL},
    {"second", (getter)DateObject_second, (setter)DateObject_ReadOnlyErr,
            "Returns the second.", NULL},
    {"minute", (getter)DateObject_minute, (setter)DateObject_ReadOnlyErr,
            "Returns the minute.", NULL},
    {"hour", (getter)DateObject_hour, (setter)DateObject_ReadOnlyErr,
            "Returns the hour.", NULL},

    {"freqstr", (getter)DateObject_freqstr, (setter)DateObject_ReadOnlyErr,
            "Returns the string representation of frequency.", NULL},
    {"datetime", (getter)DateObject_datetime, (setter)DateObject_ReadOnlyErr,
            "Returns the Date object converted to standard python datetime object",
            NULL},

    {NULL}  /* Sentinel */
};


static PyNumberMethods DateObject_as_number = {
    (binaryfunc)DateObject___add__,      /* nb_add */
    (binaryfunc)DateObject___subtract__, /* nb_subtract */
    0,                                   /* nb_multiply */
    0,                                   /* nb_divide */
    0,                                   /* nb_remainder */
    0,                                   /* nb_divmod */
    0,                                   /* nb_power */
    0,                                   /* nb_negative */
    0,                                   /* nb_positive */
    0,                                   /* nb_absolute */
    0,                                   /* nb_nonzero */
    0,                                   /* nb_invert */
    0,                                   /* nb_lshift */
    0,                                   /* nb_rshift */
    0,                                   /* nb_and */
    0,                                   /* nb_xor */
    0,                                   /* nb_or */
    0,                                   /* nb_coerce */
    (unaryfunc)DateObject___int__,       /* nb_int */
    (unaryfunc)DateObject___long__,      /* nb_long */
    (unaryfunc)DateObject___float__,     /* nb_float */
    (unaryfunc)0,                        /* nb_oct */
    (unaryfunc)0,                        /* nb_hex */
};

static PyMethodDef DateObject_methods[] = {
    {"toordinal", (PyCFunction)DateObject_toordinal, METH_NOARGS,
     DateObject_toordinal_doc},
    {"strftime", (PyCFunction)DateObject_strftime, METH_VARARGS,
     DateObject_strftime_doc},
    // deprecated alias for strftime
    {"strfmt", (PyCFunction)DateObject_strftime, METH_VARARGS,
     DateObject_strfmt_doc},
    {"asfreq", (PyCFunction)DateObject_asfreq, METH_VARARGS | METH_KEYWORDS,
     DateObject_asfreq_doc},
    {NULL}  /* Sentinel */
};


static PyTypeObject DateType = {
    PyObject_HEAD_INIT(NULL)
    0,                               /* ob_size */
    "timeseries.Date",               /* tp_name */
    sizeof(DateObject),              /* tp_basicsize */
    0,                               /* tp_itemsize */
    (destructor)DateObject_dealloc,  /* tp_dealloc */
    0,                               /* tp_print */
    0,                               /* tp_getattr */
    0,                               /* tp_setattr */
    (cmpfunc)DateObject___compare__, /* tp_compare */
    (reprfunc)DateObject___repr__,   /* tp_repr */
    &DateObject_as_number,           /* tp_as_number */
    0,                               /* tp_as_sequence */
    0,                               /* tp_as_mapping */
    (hashfunc)DateObject___hash__,   /* tp_hash */
    0,                               /* tp_call*/
    (reprfunc)DateObject___str__,    /* tp_str */
    0,                               /* tp_getattro */
    0,                               /* tp_setattro */
    0,                               /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |             /* tp_flags */
    Py_TPFLAGS_CHECKTYPES |
    Py_TPFLAGS_BASETYPE,
    "Defines a Date object, as the combination of a date and a frequency.\n"
    "Several options are available to construct a Date object explicitly:\n\n"
    "- Give appropriate values to the `year`, `month`, `day`, `quarter`, `hours`,\n"
    "  `minutes`, `seconds` arguments.\n\n"
    "  >>> td.Date(freq='Q',year=2004,quarter=3)\n"
    "  >>> td.Date(freq='D',year=2001,month=1,day=1)\n\n"
    "- Use the `string` keyword. This method uses a modified version of the\n"
    "  mx.DateTime parser submodule. More information is available in its\n"
    "  documentation.\n\n"
    "  >>> ts.Date('D', '2007-01-01')\n\n"
    "- Use the `datetime` keyword with an existing datetime.datetime object.\n\n"
    "  >>> td.Date('D', datetime=datetime.datetime.now())",  /* tp_doc */
    0,                               /* tp_traverse */
    0,                               /* tp_clear */
    0,                               /* tp_richcompare */
    0,                               /* tp_weaklistoffset */
    0,                               /* tp_iter */
    0,                               /* tp_iternext */
    DateObject_methods,              /* tp_methods */
    DateObject_members,              /* tp_members */
    DateObject_getseters,            /* tp_getset */
    0,                               /* tp_base */
    0,                               /* tp_dict */
    0,                               /* tp_descr_get */
    0,                               /* tp_descr_set */
    0,                               /* tp_dictoffset */
    (initproc)DateObject_init,       /* tp_init */
    0,                               /* tp_alloc */
    DateObject_new,                  /* tp_new */
};


///////////////////////////////////////////////////////////////////////

PyObject *
c_dates_check_freq(PyObject *self, PyObject *args) {

    PyObject *freq;
    int freq_val;

    if (!PyArg_ParseTuple(args, "O:check_freq(freq)", &freq)) return NULL;
    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) return NULL;

    return PyInt_FromLong(freq_val);
}

PyObject *
c_dates_check_freq_str(PyObject *self, PyObject *args) {

    PyObject *alias_tuple, *result, *freq_key;

    if ((freq_key = c_dates_check_freq(self, args)) == NULL) return NULL;

    alias_tuple = PyDict_GetItem(freq_dict, freq_key);
    result = PyTuple_GET_ITEM(alias_tuple, 0);

    Py_INCREF(result);

    Py_DECREF(freq_key);

    return result;
}

PyObject *
c_dates_get_freq_group(PyObject *self, PyObject *args) {

    PyObject *freq;
    int freq_val;

    if (!PyArg_ParseTuple(args, "O:get_freq_group(freq)", &freq)) return NULL;
    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) return NULL;

    return PyInt_FromLong(get_freq_group(freq_val));
}

PyObject *
c_dates_now(PyObject *self, PyObject *args) {

    PyObject *freq, *init_args, *init_kwargs;

#ifdef WIN32
    __time64_t rawtime;
#else
    time_t rawtime;
#endif
    struct tm *timeinfo;
    int freq_val;

    DateObject *secondly_date;

    if (!PyArg_ParseTuple(args, "O:now(freq)", &freq)) return NULL;

    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) return NULL;
#ifdef WIN32
    _time64(&rawtime);
#else
    time(&rawtime);
#endif


#ifdef WIN32
    timeinfo = _localtime64(&rawtime);
#else
    timeinfo = localtime(&rawtime);
#endif

    init_args = PyTuple_New(0);
    init_kwargs = PyDict_New();

    DICT_SETINT_STRKEY(init_kwargs, "freq", FR_SEC);
    DICT_SETINT_STRKEY(init_kwargs, "year", timeinfo->tm_year+1900);
    DICT_SETINT_STRKEY(init_kwargs, "month", timeinfo->tm_mon+1);
    DICT_SETINT_STRKEY(init_kwargs, "day", timeinfo->tm_mday);
    DICT_SETINT_STRKEY(init_kwargs, "hour", timeinfo->tm_hour);
    DICT_SETINT_STRKEY(init_kwargs, "minute", timeinfo->tm_min);
    DICT_SETINT_STRKEY(init_kwargs, "second", timeinfo->tm_sec);

    secondly_date = DateObject_New();
    DateObject_init(secondly_date, init_args, init_kwargs);

    Py_DECREF(init_args);
    Py_DECREF(init_kwargs);

    if (freq_val != FR_SEC) {
        DateObject *result = DateObject_New();

        long (*asfreq_func)(long, char, asfreq_info*) = NULL;
        asfreq_info af_info;

        int date_val;

        get_asfreq_info(FR_SEC, freq_val, &af_info);
        asfreq_func = get_asfreq_func(FR_SEC, freq_val, 0);

        date_val = asfreq_func(secondly_date->value, 'S', &af_info);

        Py_DECREF(secondly_date);

        result->freq = freq_val;
        result->value = date_val;

        return (PyObject*)result;

    } else { return (PyObject*)secondly_date; }
}


PyObject *
DateArray_asfreq(PyObject *self, PyObject *args)
{
    PyArrayObject *fromDates, *toDates;
    PyArrayIterObject *iterFrom, *iterTo;
    PyObject *fromDateObj, *toDateObj;
    char *relation;
    int fromFreq, toFreq;
    long fromDate, toDate;
    long (*asfreq_main)(long, char, asfreq_info*) = NULL;
    asfreq_info af_info;

    if (!PyArg_ParseTuple(args,
                "Oiis:asfreq(fromDates, fromfreq, tofreq, relation)",
                &fromDates, &fromFreq, &toFreq, &relation)) return NULL;

    get_asfreq_info(fromFreq, toFreq, &af_info);

    asfreq_main = get_asfreq_func(fromFreq, toFreq, 0);

    toDates = (PyArrayObject *)PyArray_Copy(fromDates);

    iterFrom = (PyArrayIterObject *)PyArray_IterNew((PyObject *)fromDates);
    if (iterFrom == NULL) return NULL;

    iterTo = (PyArrayIterObject *)PyArray_IterNew((PyObject *)toDates);
    if (iterTo == NULL) return NULL;

    while (iterFrom->index < iterFrom->size) {

        fromDateObj = PyArray_GETITEM(fromDates, iterFrom->dataptr);
        fromDate = PyInt_AsLong(fromDateObj);
        CHECK_ASFREQ(toDate = asfreq_main(fromDate, relation[0], &af_info));
        toDateObj = PyInt_FromLong(toDate);

        PyArray_SETITEM(toDates, iterTo->dataptr, toDateObj);

        Py_DECREF(fromDateObj);
        Py_DECREF(toDateObj);

        PyArray_ITER_NEXT(iterFrom);
        PyArray_ITER_NEXT(iterTo);
    }

    Py_DECREF(iterFrom);
    Py_DECREF(iterTo);

    return (PyObject *)toDates;

}

/**************************************************************
** The following functions are used by DateArray_getDateInfo **
** to determine how many consecutive periods will have the   **
** same result                                               **
**************************************************************/

// also used for qyear
static int __skip_periods_year(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_QTR:
            return 4;
        case FR_MTH:
            return 12;
        case FR_WK:
            return 51;
        case FR_BUS:
            return 260;
        case FR_DAY:
            return 365;
        case FR_HR:
            return 365*24;
        case FR_MIN:
            return 365*24*60;
        case FR_SEC:
            return 365*24*60*60;
        default:
            return 1;
    }
}

static int __skip_periods_quarter(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_MTH:
            return 3;
        case FR_WK:
            return 12;
        case FR_BUS:
            return 64;
        case FR_DAY:
            return 90;
        case FR_HR:
            return 90*24;
        case FR_MIN:
            return 90*24*60;
        case FR_SEC:
            return 90*24*60*60;
        default:
            return 1;
    }
}

static int __skip_periods_month(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_WK:
            return 3;
        case FR_BUS:
            return 20;
        case FR_DAY:
            return 28;
        case FR_HR:
            return 28*24;
        case FR_MIN:
            return 28*24*60;
        case FR_SEC:
            return 28*24*60*60;
        default:
            return 1;
    }
}

// also used for day_of_year, day_of_week
static int __skip_periods_day(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_HR:
            return 24;
        case FR_MIN:
            return 24*60;
        case FR_SEC:
            return 24*60*60;
        default:
            return 1;
    }
}

static int __skip_periods_week(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_BUS:
            return 5;
        case FR_DAY:
            return 7;
        case FR_HR:
            return 7*24;
        case FR_MIN:
            return 7*24*60;
        case FR_SEC:
            return 7*24*60*60;
        default:
            return 1;
    }
}

static int __skip_periods_hour(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_MIN:
            return 60;
        case FR_SEC:
            return 60*60;
        default:
            return 1;
    }
}

static int __skip_periods_minute(int freq) {

    int freq_group = get_freq_group(freq);

    switch(freq_group)
    {
        case FR_SEC:
            return 60;
        default:
            return 1;
    }
}

PyObject *
DateArray_getDateInfo(PyObject *self, PyObject *args)
{
    int freq, is_full, skip_periods, counter=1, val_changed=0;
    char *info;

    PyObject *prev_val=NULL;
    PyArrayObject *array, *newArray;
    PyArrayIterObject *iterSource, *iterResult;

    PyObject* (*getDateInfo)(DateObject*, void*) = NULL;

    if (!PyArg_ParseTuple(args, "Oisi:getDateInfo(array, freq, info, is_full)",
                                &array, &freq, &info, &is_full)) return NULL;
    newArray = (PyArrayObject *)PyArray_Copy(array);

    iterSource = (PyArrayIterObject *)PyArray_IterNew((PyObject *)array);
    iterResult = (PyArrayIterObject *)PyArray_IterNew((PyObject *)newArray);


    switch(*info)
    {
        case 'Y': //year
            getDateInfo = &DateObject_year;
            skip_periods = __skip_periods_year(freq);
            break;
        case 'F': //"fiscal" year
            getDateInfo = &DateObject_qyear;
            skip_periods = __skip_periods_year(freq);
            break;
        case 'Q': //quarter
            getDateInfo = &DateObject_quarter;
            skip_periods = __skip_periods_quarter(freq);
            break;
        case 'M': //month
            getDateInfo = &DateObject_month;
            skip_periods = __skip_periods_month(freq);
            break;
        case 'D': //day
            getDateInfo = &DateObject_day;
            skip_periods = __skip_periods_day(freq);
            break;
        case 'R': //day of year
            getDateInfo = &DateObject_day_of_year;
            skip_periods = __skip_periods_day(freq);
            break;
        case 'W': //day of week
            getDateInfo = &DateObject_day_of_week;
            skip_periods = __skip_periods_day(freq);
            break;
        case 'I': //week of year
            getDateInfo = &DateObject_week;
            skip_periods = __skip_periods_week(freq);
            break;
        case 'H': //hour
            getDateInfo = &DateObject_hour;
            skip_periods = __skip_periods_hour(freq);
            break;
        case 'T': //minute
            getDateInfo = &DateObject_minute;
            skip_periods = __skip_periods_minute(freq);
            break;
        case 'S': //second
            getDateInfo = &DateObject_second;
            skip_periods = 1;
            break;
        default:
            return NULL;
    }

    {
        DateObject *curr_date;
        PyObject *val, *dInfo;

        while (iterSource->index < iterSource->size) {

            if ((val_changed == 0) ||
                (is_full == 0) ||
                (prev_val == NULL) ||
                (counter >= skip_periods)) {

                   val = PyArray_GETITEM(array, iterSource->dataptr);
                   curr_date = DateObject_FromFreqAndValue(freq, PyInt_AsLong(val));
                   dInfo = getDateInfo(curr_date, NULL);

                   if ((prev_val != NULL) &&
                       (PyInt_AsLong(prev_val) != PyInt_AsLong(dInfo))) {
                       val_changed = 1;
                       counter = 0;
                   }

                   Py_DECREF(val);
                   Py_DECREF(curr_date);

                   if (prev_val != NULL) {
                       Py_DECREF(prev_val);
                   }

                   prev_val = dInfo;
            }

            PyArray_SETITEM(newArray, iterResult->dataptr, dInfo);

            PyArray_ITER_NEXT(iterSource);
            PyArray_ITER_NEXT(iterResult);

            counter += 1;
        }
    }

    if (prev_val != NULL) {
        Py_DECREF(prev_val);
    }
    Py_DECREF(iterSource);
    Py_DECREF(iterResult);

    return (PyObject *) newArray;
}


void import_c_dates(PyObject *m)
{

    if (PyType_Ready(&DateType) < 0) return;

    DateCalc_Error =
        PyErr_NewException("c_dates.DateCalc_Error", NULL, NULL);
    DateCalc_RangeError =
        PyErr_NewException("c_dates.DateCalc_RangeError", NULL, NULL);

    import_array();
    PyDateTime_IMPORT;

    Py_INCREF(&DateType);
    PyModule_AddObject(m, "Date", (PyObject *)(&DateType));

    if(build_freq_dict() == INT_ERR_CODE) {
        PyErr_SetString(                    \
            PyExc_ImportError,              \
            "initialization of module timeseries.c_dates failed");
        return;
    };

    PyModule_AddObject(m, "freq_dict", freq_dict);
    PyModule_AddObject(m, "freq_dict_rev", freq_dict_rev);
    PyModule_AddObject(m, "freq_constants", freq_constants);

    PyModule_AddObject(m, "DateCalc_Error", DateCalc_Error);
    PyModule_AddObject(m, "DateCalc_RangeError", DateCalc_RangeError);

}
