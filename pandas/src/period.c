#include "period.h"
#include "limits.h"
#include "numpy/ndarraytypes.h"

/*
 * Borrowed and derived code from scikits.timeseries that we will expose via
 * Cython to pandas. This primarily concerns period representation and
 * frequency conversion routines.
 */

/* see end of file for stuff pandas uses (search for 'pandas') */

/* ------------------------------------------------------------------
 * Code derived from scikits.timeseries
 * ------------------------------------------------------------------*/

static asfreq_info NULL_AF_INFO;

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

/* Return 1/0 iff year points to a leap year in calendar. */
static int dInfoCalc_Leapyear(long year, int calendar)
{
    if (calendar == GREGORIAN_CALENDAR) {
        return (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));
    } else {
        return (year % 4 == 0);
    }
}

/* Return the day of the week for the given absolute date. */
static int dInfoCalc_DayOfWeek(long absdate)
{
    int day_of_week;

    if (absdate >= 1) {
        day_of_week = (absdate - 1) % 7;
    } else {
        day_of_week = 6 - ((-absdate) % 7);
    }
    return day_of_week;
}

static int monthToQuarter(int month) { return ((month-1)/3)+1; }

/* Return the year offset, that is the absolute date of the day
   31.12.(year-1) in the given calendar.

   Note:
   For the Julian calendar we shift the absdate (which is measured
   using the Gregorian Epoch) value by two days because the Epoch
   (0001-01-01) in the Julian calendar lies 2 days before the Epoch in
   the Gregorian calendar. */
static int dInfoCalc_YearOffset(long year, int calendar)
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
    Py_Error(PyExc_ValueError, "unknown calendar");
 onError:
    return INT_ERR_CODE;
}

/* Set the instance's value using the given date and time. calendar may be set
 * to the flags: GREGORIAN_CALENDAR, JULIAN_CALENDAR to indicate the calendar
 * to be used. */

static int dInfoCalc_SetFromDateAndTime(struct date_info *dinfo,
        int year, int month, int day, int hour, int minute, double second,
        int calendar)
{

    /* Calculate the absolute date */
    {
        int leap;
        long yearoffset,absdate;

        /* Range check */
        Py_AssertWithArg(year > -(INT_MAX / 366) && year < (INT_MAX / 366),
                 PyExc_ValueError,
                 "year out of range: %i",
                 year);

        /* Is it a leap year ? */
        leap = dInfoCalc_Leapyear(year, calendar);

        /* Negative month values indicate months relative to the years end */
        if (month < 0) month += 13;
        Py_AssertWithArg(month >= 1 && month <= 12,
                 PyExc_ValueError,
                 "month out of range (1-12): %i",
                 month);

        /* Negative values indicate days relative to the months end */
        if (day < 0) day += days_in_month[leap][month - 1] + 1;
        Py_AssertWithArg(day >= 1 && day <= days_in_month[leap][month - 1],
                 PyExc_ValueError,
                 "day out of range: %i",
                 day);

        yearoffset = dInfoCalc_YearOffset(year, calendar);
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
                PyExc_ValueError,
                "hour out of range (0-23): %i",
                hour);
        Py_AssertWithArg(minute >= 0 && minute <= 59,
                PyExc_ValueError,
                "minute out of range (0-59): %i",
                minute);
        Py_AssertWithArg(second >= (double)0.0 &&
                (second < (double)60.0 ||
                (hour == 23 && minute == 59 &&
                second < (double)61.0)),
                PyExc_ValueError,
                "second out of range (0.0 - <60.0; <61.0 for 23:59): %f",
                second);

        dinfo->abstime = (double)(hour*3600 + minute*60) + second;

        dinfo->hour = hour;
        dinfo->minute = minute;
        dinfo->second = second;
    }
    return 0;

 onError:
    return INT_ERR_CODE;
}

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
        Py_Error(PyExc_ValueError, "unknown calendar");
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
    return INT_ERR_CODE;
}

///////////////////////////////////////////////

// frequency specifc conversion routines
// each function must take an integer fromDate and
// a char relation ('S' or 'E' for 'START' or 'END')
///////////////////////////////////////////////////////////////////////

// helpers for frequency conversion routines //

static long DtoB_weekday(long fromDate) {
    return (((fromDate) / 7) * 5) + (fromDate)%7;
}

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
    if (dInfoCalc_SetFromDateAndTime(&tempDate, y, m, d, 0, 0, 0, GREGORIAN_CALENDAR)) {
        return INT_ERR_CODE;
    }
    return tempDate.absdate;
}

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

    if (DtoQ_yq(fromDate, af_info, &year, &quarter) == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    return (long)((year - 1) * 4 + quarter);
}

static long asfreq_DtoM(long fromDate, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;
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
    if (dInfoCalc_SetFromAbsDate(&dinfo, fromDate, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;

    if (dinfo.day_of_week > 4) {
        return INT_ERR_CODE;
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
    } else { return INT_ERR_CODE; }
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

static long nofunc(long fromDate, char relation, asfreq_info *af_info) { return INT_ERR_CODE; }
static long no_op(long fromDate, char relation, asfreq_info *af_info) { return fromDate; }

// end of frequency specific conversion routines

static int get_freq_group(int freq) { return (freq/1000)*1000; }

static int calc_a_year_end(int freq, int group) {
    int result = (freq - group) % 12;
    if (result == 0) {return 12;}
    else {return result;}
}

static int calc_week_end(int freq, int group) {
    return freq - group;
}

static void get_asfreq_info(int fromFreq, int toFreq, asfreq_info *af_info) {
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


static freq_conv_func get_asfreq_func(int fromFreq, int toFreq, int forConvert)
{
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
                case FR_MTH: return &no_op;
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
                case FR_BUS: return &no_op;
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
                case FR_HR: return &no_op;
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
                case FR_MIN: return &no_op;
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
                case FR_SEC: return &no_op;
                default: return &nofunc;
            }
        default: return &nofunc;
    }
}

double getAbsTime(int freq, long dailyDate, long originalDate) {

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
		  return 0; // 24*60*60 - 1;
    }

    startOfDay = asfreq_DtoHIGHFREQ(dailyDate, 'S', periodsPerDay);
    return (24*60*60)*((double)(originalDate - startOfDay))/((double)periodsPerDay);
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
             PyExc_ValueError,
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
    return INT_ERR_CODE;
}

/* ------------------------------------------------------------------
 * New pandas API-helper code, to expose to cython
 * ------------------------------------------------------------------*/

long asfreq(long period_ordinal, int freq1, int freq2, char relation)
{
    freq_conv_func func = get_asfreq_func(freq1, freq2, 0);

    asfreq_info finfo;
    get_asfreq_info(freq1, freq2, &finfo);

    long val = (*func)(period_ordinal, relation, &finfo);

    if (val == INT_ERR_CODE) {
        Py_Error(PyExc_ValueError, "Unable to convert to desired frequency.");
		goto onError;
	}
    return val;
onError:
    return INT_ERR_CODE;
}

/* generate an ordinal in period space */
long get_period_ordinal(int year, int month, int day,
                      int hour, int minute, int second,
                      int freq)
{
    int freq_group = get_freq_group(freq);
    int quarter=((month-1)/3)+1;

    if (freq == FR_SEC) {
        long absdays, delta;
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - HIGHFREQ_ORIG);
        return (long)(delta*86400 + hour*3600 + minute*60 + second + 1);
    }

    if (freq == FR_MIN) {
        long absdays, delta;
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - HIGHFREQ_ORIG);
        return (long)(delta*1440 + hour*60 + minute + 1);
    }

    if (freq == FR_HR) {
        long absdays, delta;
        if ((absdays = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        delta = (absdays - HIGHFREQ_ORIG);
        return (long)(delta*24 + hour + 1);
    }

    if (freq == FR_DAY)
    {
        return (long)absdate_from_ymd(year, month, day);
    }

    if (freq == FR_UND)
    {
        return (long)absdate_from_ymd(year, month, day);
    }

    if (freq == FR_BUS)
    {
        long weeks, days;
        if((days = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        weeks = days/7;
        return (long)(days - weeks*2);
    }

    if (freq_group == FR_WK)
    {
        long adj_ordinal, ordinal, day_adj;
        if((ordinal = (long)absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        day_adj = (7 - (freq - FR_WK)) % 7;
        adj_ordinal = ordinal + ((7 - day_adj) - ordinal % 7) % 7;
        return adj_ordinal/7;
    }

    if (freq == FR_MTH)
    {
        return (year-1)*12 + month;
    }

    if (freq_group == FR_QTR)
    {
        return (year-1)*4 + quarter;
    }

    if (freq_group == FR_ANN)
    {
        return year;
    }

    Py_Error(PyExc_RuntimeError, "Unable to generate frequency ordinal");

onError:
    return INT_ERR_CODE;
}

/*
    Returns the proleptic Gregorian ordinal of the date, as an integer.
    This corresponds to the number of days since Jan., 1st, 1AD.
    When the instance has a frequency less than daily, the proleptic date
    is calculated for the last day of the period.
*/

long get_python_ordinal(long period_ordinal, int freq)
{
    if (freq == FR_DAY)
        return period_ordinal;

    long (*toDaily)(long, char, asfreq_info*) = NULL;
    asfreq_info af_info;

    toDaily = get_asfreq_func(freq, FR_DAY, 0);
    get_asfreq_info(freq, FR_DAY, &af_info);

    return toDaily(period_ordinal, 'E', &af_info);
}

char *str_replace(const char *s, const char *old, const char *new) {
    char *ret;
    int i, count = 0;
    size_t newlen = strlen(new);
    size_t oldlen = strlen(old);

    for (i = 0; s[i] != '\0'; i++) {
        if (strstr(&s[i], old) == &s[i]) {
           count++;
           i += oldlen - 1;
        }
    }

    ret = PyArray_malloc(i + 1 + count * (newlen - oldlen));
    if (ret == NULL) {return (char *)PyErr_NoMemory();}

    i = 0;
    while (*s) {
        if (strstr(s, old) == s) {
            strcpy(&ret[i], new);
            i += newlen;
            s += oldlen;
        } else {
            ret[i++] = *s++;
        }
    }
    ret[i] = '\0';

    return ret;
}

// function to generate a nice string representation of the period
// object, originally from DateObject_strftime

char *skts_strftime(long value, int freq, PyObject *args)
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

    if (!PyArg_ParseTuple(args, "s:strftime(fmt)", &orig_fmt_str))
        return NULL;

    toDaily = get_asfreq_func(freq, FR_DAY, 0);
    get_asfreq_info(freq, FR_DAY, &af_info);

    absdate = toDaily(value, 'E', &af_info);
    abstime = getAbsTime(freq, absdate, value);

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
    if ((result = PyArray_malloc(result_len * sizeof(char))) == NULL) {
        return (char*)PyErr_NoMemory();
    }

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

                if (get_freq_group(freq) == FR_QTR) {
                    qtr_freq = freq;
                } else { qtr_freq = FR_QTR; }
                get_asfreq_info(FR_DAY, qtr_freq, &af_info);

                if(DtoQ_yq(absdate, &af_info, &year, &quarter) == INT_ERR_CODE)
                { return NULL; }

                if(strcmp(extra_fmts[i][0], "%q") == 0) {
                    if ((extra_str = PyArray_malloc(2 * sizeof(char))) == NULL) {
                        free(tmp_str);
                        return (char *)PyErr_NoMemory();
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
                        return (char *)PyErr_NoMemory();
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

    return result;
}

char *period_to_string(long value, int freq)
{
    int freq_group = get_freq_group(freq);
    PyObject *string_arg;
    char *retval;

    string_arg = NULL;
    if (freq_group == FR_UND) {
        int digits = log10(value) + 1;
        if ((retval = PyArray_malloc(digits * sizeof(char))) == NULL) {
            return (char *)PyErr_NoMemory();
        }
        sprintf(retval, "%ld", value);
        return retval;
    }
    else if (freq_group == FR_ANN) { string_arg = Py_BuildValue("(s)", "%Y"); }
    else if (freq_group == FR_QTR) { string_arg = Py_BuildValue("(s)", "%FQ%q"); }
    else if (freq_group == FR_MTH) { string_arg = Py_BuildValue("(s)", "%b-%Y"); }
    else if (freq_group == FR_DAY ||
             freq_group == FR_BUS ||
             freq_group == FR_WK) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y"); }
    else if (freq_group == FR_HR) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:00"); }
    else if (freq_group == FR_MIN) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M"); }
    else if (freq_group == FR_SEC) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M:%S"); }

    if (string_arg == NULL) { return (char *)NULL; }

    retval = skts_strftime(value, freq, string_arg);
    Py_DECREF(string_arg);

    return retval;
}

char *period_to_string2(long value, int freq, char *fmt)
{
    PyObject *string_arg;
    char *retval;
    string_arg = Py_BuildValue("(s)", fmt);
    if (string_arg == NULL) { return (char *)NULL; }
    retval = skts_strftime(value, freq, string_arg);
    Py_DECREF(string_arg);
    return retval;
}

static int _quarter_year(long ordinal, int freq, int *year, int *quarter) {
    asfreq_info af_info;
    int qtr_freq;

    ordinal = get_python_ordinal(ordinal, freq);

    if (get_freq_group(freq) == FR_QTR)
        qtr_freq = freq;
    else
        qtr_freq = FR_QTR;

    get_asfreq_info(FR_DAY, qtr_freq, &af_info);

    if (DtoQ_yq(ordinal, &af_info, year, quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;

    if ((qtr_freq % 1000) > 12)
        *year -= 1;

    return 0;
}

static int _ISOWeek(struct date_info *dinfo)
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

int get_date_info(long ordinal, int freq, struct date_info *dinfo)
{
    long absdate = get_python_ordinal(ordinal, freq);
    double abstime = getAbsTime(freq, absdate, ordinal);

    if(dInfoCalc_SetFromAbsDateTime(dinfo, absdate, abstime, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;

    return 0;
}

int pyear(long ordinal, int freq) {
    struct date_info dinfo;
    get_date_info(ordinal, freq, &dinfo);
    return dinfo.year;
}

int pqyear(long ordinal, int freq) {
    int year, quarter;
    if( _quarter_year(ordinal, freq, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return year;
}

int pquarter(long ordinal, int freq) {
    int year, quarter;
    if(_quarter_year(ordinal, freq, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return quarter;
}

int pmonth(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.month;
}

int pday(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day;
}

int pweekday(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_week;
}

int pday_of_week(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_week;
}

int pday_of_year(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_year;
}

int pweek(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return _ISOWeek(&dinfo);
}

int phour(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.hour;
}

int pminute(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.minute;
}

int psecond(long ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return (int)dinfo.second;
}
