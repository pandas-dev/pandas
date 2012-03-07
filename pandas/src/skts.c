#include "skts.h"
#include "limits.h"

/* 
 * Borrowed and derived code from scikits.timeseries that we will expose via
 * Cython to pandas. This primarily concerns interval frequency conversion
 * routines.
 */

/* see end of file for stuff pandas uses */

/* ------------------------------------------------------------------
 * Code derived from skts
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
    return -1;
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
    return -1;
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
    return -1;
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

    if (DtoQ_yq(fromDate, af_info, &year, &quarter) == INT_ERR_CODE)
    { 
        return INT_ERR_CODE; 
    }

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

/* ------------------------------------------------------------------
 * New pandas API-helper code, to expose to cython 
 * ------------------------------------------------------------------*/

int frequency_conversion(long dtordinal, int freq1, int freq2, char relation)
{
    freq_conv_func func = get_asfreq_func(freq1, freq2, 1);

    asfreq_info finfo;
    get_asfreq_info(freq1, freq2, &finfo);

    return (*func)(dtordinal, relation, &finfo);
}

/* generate an ordinal in skts space */
long get_skts_ordinal(int year, int month, int day,
                     int hour, int minute, int second,
                     int freq)
{
    int freq_group = get_freq_group(freq);
    int quarter=((month-1)/3)+1;

    if (freq == FR_SEC) {
        long absdays, delta;
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - HIGHFREQ_ORIG);
        return (int)(delta*86400 + hour*3600 + minute*60 + second + 1);
    }

    if (freq == FR_MIN) {
        long absdays, delta;
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - HIGHFREQ_ORIG);
        return (int)(delta*1440 + hour*60 + minute + 1);
    }

    if (freq == FR_HR) {
        long absdays, delta;
        if ((absdays = absdate_from_ymd(year, month, day)) == INT_ERR_CODE) 
        {
            goto onError;
        }
        delta = (absdays - HIGHFREQ_ORIG);
        return (int)(delta*24 + hour + 1);
    }

    if (freq == FR_DAY)
    {
        return (int)absdate_from_ymd(year, month, day);
    }

    if (freq == FR_UND)
    {
        return (int)absdate_from_ymd(year, month, day);
    }

    if (freq == FR_BUS)
    {
        long weeks, days;
        if((days = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        weeks = days/7;
        return (int)(days - weeks*2);
    }

    if (freq_group == FR_WK)
    {
        int adj_ordinal, ordinal, day_adj;
        if((ordinal = (int)absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
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
        if ((freq - freq_group) > 12) {
            // quarterly frequency with year determined by ending period
            return year*4 + quarter;
        } else {
            /* quarterly frequency with year determined by ending period
            or has December year end*/
            return (year-1)*4 + quarter;
        }
    }

    if (freq_group == FR_ANN)
    {
        return year;
    }

onError:
    Py_Error(PyExc_Exception, "Unable to generate frequency ordinal");
    return INT_ERR_CODE;
}

/*
    Returns the proleptic Gregorian ordinal of the date, as an integer.
    This corresponds to the number of days since Jan., 1st, 1AD.
    When the instance has a frequency less than daily, the proleptic date 
    is calculated for the last day of the period.
*/

long get_python_ordinal(long skts_ordinal, int freq)
{
    if (freq == FR_DAY) 
        return skts_ordinal;

    long (*toDaily)(long, char, asfreq_info*) = NULL;
    asfreq_info af_info;

    toDaily = get_asfreq_func(freq, FR_DAY, 0);
    get_asfreq_info(freq, FR_DAY, &af_info);

    return toDaily(skts_ordinal, 'E', &af_info);
}

