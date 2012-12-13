#include "period.h"


/*
 * Borrowed and derived code from scikits.timeseries that we will expose via
 * Cython to pandas. This primarily concerns period representation and
 * frequency conversion routines.
 */

/* see end of file for stuff pandas uses (search for 'pandas') */

/* ------------------------------------------------------------------
 * Code derived from scikits.timeseries
 * ------------------------------------------------------------------*/


static i8 mod_compat(i8 x, i8 m) {
    i8 result = x % m;

    if (result < 0)
        result += m;

    return result;
}

static i8 floordiv(i8 x, i8 divisor) {
    i8 x_div_d = x / divisor;

    if (x < 0 && mod_compat(x, divisor))
        --x_div_d;

    return x_div_d;
}

static asfreq_info NULL_AF_INFO;

/* Table with day offsets for each month (0-based, without and with leap) */
static i8 month_offset[2][13] = {
    { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 },
    { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 }
};

/* Table of number of days in a month (0-based, without and with leap) */
static i8 days_in_month[2][12] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

/* Return 1/0 iff year points to a leap year in calendar. */
static i8 dInfoCalc_Leapyear(i8 year, i8 calendar)
{
    i8 ymod4_is0 = year % 4 == 0;

    if (calendar == GREGORIAN)
        return ymod4_is0 && (year % 100 != 0 || year % 400 == 0);

    return ymod4_is0;
}

/* Return the day of the week for the given absolute date. */
static i8 dInfoCalc_DayOfWeek(i8 absdate)
{
    return absdate >= 1 ? (absdate - 1) % 7 : 6 - (-absdate % 7);
}

static i8 monthToQuarter(i8 month) { return (month - 1) / 3 + 1; }

/* Return the year offset, that is the absolute date of the day
   31.12.(year-1) in the given calendar.

   Note:
   For the Julian calendar we shift the absdate (which is measured
   using the Gregorian Epoch) value by two days because the Epoch
   (0001-01-01) in the Julian calendar lies 2 days before the Epoch in
   the Gregorian calendar. */
static i8 dInfoCalc_YearOffset(i8 year, i8 calendar)
{
    --year;

    if (calendar == GREGORIAN)
        if (year >= 0 || -1 / 4 == -1)
            return year * 365 + year / 4 - year / 100 + year / 400;
        else
            return year * 365 + (year - 3) / 4 - (year - 99) / 100 + (year - 399) / 400;
    else if (calendar == JULIAN)
        if (year >= 0 || -1 / 4 == -1)
            return year * 365 + year / 4 - 2;
        else
            return year * 365 + (year - 3) / 4 - 2;
    else
        Py_Error(PyExc_ValueError, "unknown calendar");
onError:
    return INT_ERR_CODE;
}

/* Set the instance's value using the given date and time. calendar may be set
 * to the flags: GREGORIAN, JULIAN to indicate the calendar
 * to be used. */

static i8 dInfoCalc_SetFromDateAndTime(struct date_info *dinfo, i8 year,
                                       i8 month, i8 day, i8 hour, i8 minute,
                                       i8 second, i8 microsecond, i8 calendar)
{

    /* Calculate the absolute date */
    {
        i8 leap, absdate, yearoffset;

        /* Range check */
        Py_AssertWithArg(year > -(INT_MAX / 366) && year < (INT_MAX / 366),
                         PyExc_ValueError,
                         "year out of range: %li",
                         year);

        /* Is it a leap year ? */
        leap = dInfoCalc_Leapyear(year, calendar);

        /* Negative month values indicate months relative to the years end */
        if (month < 0)
            month += 13;

        Py_AssertWithArg(month >= 1 && month <= 12,
                         PyExc_ValueError,
                         "month out of range (1-12): %li",
                         month);

        /* Negative values indicate days relative to the months end */
        if (day < 0)
            day += days_in_month[leap][month - 1] + 1;

        Py_AssertWithArg(day >= 1 && day <= days_in_month[leap][month - 1],
                         PyExc_ValueError,
                         "day out of range: %li",
                         day);

        yearoffset = dInfoCalc_YearOffset(year, calendar);

        if (PyErr_Occurred())
            goto onError;

        absdate = day + month_offset[leap][month - 1] + yearoffset;

        dinfo->absdate = absdate;

        dinfo->year = year;
        dinfo->month = month;
        dinfo->quarter = (month - 1) / 3 + 1;
        dinfo->day = day;

        dinfo->day_of_week = dInfoCalc_DayOfWeek(absdate);
        dinfo->day_of_year = absdate - yearoffset;

        dinfo->calendar = calendar;
    }

    /* Calculate the absolute time */
    {
        Py_AssertWithArg(hour >= 0 && hour <= 23,
                         PyExc_ValueError,
                         "hour out of range (0-23): %li",
                         hour);
        Py_AssertWithArg(minute >= 0 && minute <= 59,
                         PyExc_ValueError,
                         "minute out of range (0-59): %li",
                         minute);
        Py_AssertWithArg(second >= 0 &&
                         (second < 60L || (hour == 23 && minute == 59 &&
                                           second < 61)),
                         PyExc_ValueError,
                         "second out of range (0 - <60L; <61 for 23:59): %li",
                         second);
        Py_AssertWithArg(microsecond >= 0 && (microsecond < 1000000 ||
                                              (hour == 23 && minute == 59 &&
                                               second < 61 &&
                                               microsecond < 1000001)),
                         PyExc_ValueError,
                         "microsecond out of range (0 - <100000; <100001 for 23:59:59): %li",
                         microsecond);

        dinfo->abstime = hour * US_PER_HOUR + minute * US_PER_MINUTE +
            second * US_PER_SECOND + microsecond;

        dinfo->hour = hour;
        dinfo->minute = minute;
        dinfo->second = second;
        dinfo->microsecond = microsecond;
    }

    return 0;

onError:
    return INT_ERR_CODE;
}

/* Sets the date part of the date_info struct using the indicated
   calendar.

   XXX This could also be done using some i8eger arithmetics rather
   than with this iterative approach... */
static
i8 dInfoCalc_SetFromAbsDate(register struct date_info *dinfo,
                            i8 absdate, i8 calendar)
{
    register i8 year;
    i8 yearoffset, leap, dayoffset;
    i8 *monthoffset = NULL;

    /* Approximate year */
    switch (calendar) {
    case GREGORIAN:
        year = absdate / 365.2425;
        break;
    case JULIAN:
        year = absdate / 365.25;
        break;
    default:
        Py_Error(PyExc_ValueError, "unknown calendar");
        break;
    }

    if (absdate > 0)
        ++year;

    /* Apply corrections to reach the correct year */
    while (1) {

        /* Calculate the year offset */
        yearoffset = dInfoCalc_YearOffset(year, calendar);
        if (PyErr_Occurred())
            goto onError;

        /* Backward correction: absdate must be greater than the
           yearoffset */
        if (yearoffset >= absdate) {
            --year;
            continue;
        }

        dayoffset = absdate - yearoffset;
        leap = dInfoCalc_Leapyear(year,calendar);

        /* Forward correction: non leap years only have 365 days */
        if (dayoffset > 365 && !leap) {
            ++year;
            continue;
        }

        break;
    }

    dinfo->year = year;
    dinfo->calendar = calendar;

    /* Now iterate to find the month */
    monthoffset = month_offset[leap];
    {
        register i8 month;

        for (month = 1; month < 13; ++month)
            if (monthoffset[month] >= dayoffset)
                break;

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
// each function must take an i8eger fromDate and
// a char relation ('S' or 'E' for 'START' or 'END')
///////////////////////////////////////////////////////////////////////

// helpers for frequency conversion routines //

static i8 DtoB_weekday(i8 absdate) {
    return (absdate / 7) * 5 + absdate % 7 - BDAY_OFFSET;
}

static i8 DtoB_WeekendToMonday(i8 absdate, i8 day_of_week) {
    if (day_of_week > 4) {
        //change to Monday after weekend
        absdate += 7 - day_of_week;
    }
    return DtoB_weekday(absdate);
}

static i8 DtoB_WeekendToFriday(i8 absdate, i8 day_of_week) {
    if (day_of_week > 4) {
        //change to friday before weekend
        absdate -= day_of_week - 4;
    }
    return DtoB_weekday(absdate);
}

static i8 absdate_from_ymd(i8 y, i8 m, i8 d) {
    struct date_info td;

    if (dInfoCalc_SetFromDateAndTime(&td, y, m, d, 0, 0, 0, 0, GREGORIAN))
        return INT_ERR_CODE;

    return td.absdate;
}

//************ FROM DAILY ***************

static i8 asfreq_DtoA(i8 ordinal, const char* relation, asfreq_info *af_info) {

    struct date_info dinfo;

    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN))
        return INT_ERR_CODE;

    if (dinfo.month > af_info->to_a_year_end)
        return dinfo.year + 1 - BASE_YEAR;
	else
        return dinfo.year - BASE_YEAR;
}

static i8 DtoQ_yq(i8 ordinal, asfreq_info *af_info, i8 *year, i8 *quarter) {
    struct date_info dinfo;

    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN))
        return INT_ERR_CODE;

    if (af_info->to_q_year_end != 12) {
        dinfo.month -= af_info->to_q_year_end;

        if (dinfo.month <= 0)
            dinfo.month += 12;
        else
            ++dinfo.year;

        dinfo.quarter = monthToQuarter(dinfo.month);
    }

    *year = dinfo.year;
    *quarter = dinfo.quarter;

    return 0;
}


static i8 asfreq_DtoQ(i8 ordinal, const char* relation, asfreq_info *af_info) {

    i8 year, quarter;

    if (DtoQ_yq(ordinal, af_info, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;

    return (year - BASE_YEAR) * 4 + quarter - 1;
}

static i8 asfreq_DtoM(i8 ordinal, const char* relation, asfreq_info *af_info) {

    struct date_info dinfo;

    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN))
        return INT_ERR_CODE;

    return (dinfo.year - BASE_YEAR) * 12 + dinfo.month - 1;
}

static i8 asfreq_DtoW(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return (ordinal + ORD_OFFSET - (1 + af_info->to_week_end)) / 7
        + 1 - WEEK_OFFSET;
}

static i8 asfreq_DtoB(i8 ordinal, const char* relation, asfreq_info *af_info) {
    struct date_info dinfo;

    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN))
        return INT_ERR_CODE;

    if (!strcmp(relation, "S"))
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
    else
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
}

// needed for getDateInfo function
static i8 asfreq_DtoD(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return ordinal;
}

static i8 asfreq_DtoHIGHFREQ(i8 ordinal, const char* relation, i8 per_day) {
    if (!strcmp(relation, "S"))
        return ordinal * per_day;
    else
        return (ordinal + 1L) * per_day - 1L;
}

static i8 asfreq_DtoH(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoHIGHFREQ(ordinal, relation, 24L);
}

static i8 asfreq_DtoT(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoHIGHFREQ(ordinal, relation, 24 * 60L);
}

static i8 asfreq_DtoS(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoHIGHFREQ(ordinal, relation, 24 * 60L * 60L);
}

static i8 asfreq_DtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoHIGHFREQ(ordinal, relation, US_PER_DAY);
}

//************ FROM SECONDLY ***************

static i8 asfreq_StoD(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return ordinal / (60L * 60L * 24L);
}

static i8 asfreq_StoA(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_StoQ(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_StoM(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoM(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation,
                     &NULL_AF_INFO); }

static i8 asfreq_StoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation,
                     af_info); }

static i8 asfreq_StoB(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoB(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation,
                     &NULL_AF_INFO); }


static i8 asfreq_StoT(i8 ordinal, const char* relation, asfreq_info *af_info) {
	return ordinal / 60L;
}

static i8 asfreq_StoH(i8 ordinal, const char* relation, asfreq_info *af_info) {
	return ordinal / (60L * 60L);
}

//************ FROM MINUTELY ***************

static i8 asfreq_TtoD(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return ordinal / (60L * 24L); }

static i8 asfreq_TtoA(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoA(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static i8 asfreq_TtoQ(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoQ(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static i8 asfreq_TtoM(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoM(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static i8 asfreq_TtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static i8 asfreq_TtoB(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoB(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static i8 asfreq_TtoH(i8 ordinal, const char* relation, asfreq_info *af_info) {
	return ordinal / 60L;
}

static i8 asfreq_TtoS(i8 ordinal, const char* relation, asfreq_info *af_info) {
    i8 out = ordinal * 60L;

    if (strcmp(relation, "S"))
        out += 59;

    return out;
}

//************ FROM HOURLY ***************

static i8 asfreq_HtoD(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return ordinal / 24; }
static i8 asfreq_HtoA(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoA(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static i8 asfreq_HtoQ(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoQ(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static i8 asfreq_HtoM(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoM(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static i8 asfreq_HtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static i8 asfreq_HtoB(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoB(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

// calculation works out the same as TtoS, so we just call that function for HtoT
static i8 asfreq_HtoT(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_TtoS(ordinal, relation, &NULL_AF_INFO);
}

static i8 asfreq_HtoS(i8 ordinal, const char* relation, asfreq_info *af_info) {
    i8 is_S = !strcmp(relation, "S");
    return is_S ? ordinal * 60 * 60 : (ordinal + 1) * 60 * 60 - 1;
}

//************ FROM BUSINESS ***************

static i8 asfreq_BtoD(i8 ordinal, const char* relation, asfreq_info *af_info)
{
    i8 ord = ordinal;
    ord += BDAY_OFFSET;
    return (ord - 1) / 5 * 7 + mod_compat(ord - 1, 5) + 1 - ORD_OFFSET;
}

static i8 asfreq_BtoA(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoA(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static i8 asfreq_BtoQ(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoQ(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static i8 asfreq_BtoM(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoM(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static i8 asfreq_BtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static i8 asfreq_BtoH(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoH(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static i8 asfreq_BtoT(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoT(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static i8 asfreq_BtoS(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoS(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM WEEKLY ***************

static i8 asfreq_WtoD(i8 ordinal, const char* relation, asfreq_info *af_info) {
    i8 k = 0, ord = ordinal;
	ord += WEEK_OFFSET;
    ord *= 7;

    if (!strcmp(relation, "S"))
        k = -6;

    return ord + k + af_info->from_week_end - ORD_OFFSET;
}

static i8 asfreq_WtoA(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_WtoD(ordinal, "E", af_info), relation, af_info);
}

static i8 asfreq_WtoQ(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_WtoD(ordinal, "E", af_info), relation, af_info);
}

static i8 asfreq_WtoM(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_WtoD(ordinal, "E", af_info), relation,
                       &NULL_AF_INFO);
}

static i8 asfreq_WtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_WtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_WtoB(i8 ordinal, const char* relation, asfreq_info *af_info) {

    struct date_info dinfo;
    i8 wtod = asfreq_WtoD(ordinal, relation, af_info) + ORD_OFFSET;

    if (dInfoCalc_SetFromAbsDate(&dinfo, wtod, GREGORIAN))
        return INT_ERR_CODE;

    i8 (*f)(i8 absdate, i8 day_of_week) = NULL;
    f = !strcmp(relation, "S") ? DtoB_WeekendToMonday : DtoB_WeekendToFriday;
    return f(dinfo.absdate, dinfo.day_of_week);
}

static i8 asfreq_WtoH(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoH(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static i8 asfreq_WtoT(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoT(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static i8 asfreq_WtoS(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoS(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }

//************ FROM MONTHLY ***************
static void MtoD_ym(i8 ordinal, i8 *y, i8 *m) {
    *y = floordiv(ordinal, 12) + BASE_YEAR;
    *m = mod_compat(ordinal, 12) + 1;
}


static i8 asfreq_MtoD(i8 ordinal, const char* relation, asfreq_info *af_info) {

    i8 absdate, y, m;

    if (!strcmp(relation, "S")) {
        MtoD_ym(ordinal, &y, &m);

        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE)
            return INT_ERR_CODE;

        return absdate - ORD_OFFSET;
    } else {
        MtoD_ym(ordinal + 1, &y, &m);

        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE)
            return INT_ERR_CODE;

        return absdate - 1 - ORD_OFFSET;
    }
}

static i8 asfreq_MtoA(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_MtoD(ordinal, "E", &NULL_AF_INFO), relation, af_info); }

static i8 asfreq_MtoQ(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_MtoD(ordinal, "E", &NULL_AF_INFO), relation, af_info); }

static i8 asfreq_MtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static i8 asfreq_MtoB(i8 ordinal, const char* relation, asfreq_info *af_info) {

    struct date_info dinfo;
    i8 mtod = asfreq_MtoD(ordinal, relation, &NULL_AF_INFO) + ORD_OFFSET;

    if (dInfoCalc_SetFromAbsDate(&dinfo, mtod, GREGORIAN))
        return INT_ERR_CODE;

    i8 (*f)(i8 absdate, i8 day_of_week);
    f = !strcmp(relation, "S") ? DtoB_WeekendToMonday : DtoB_WeekendToFriday;

    return f(dinfo.absdate, dinfo.day_of_week);
}

static i8 asfreq_MtoH(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoH(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static i8 asfreq_MtoT(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoT(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static i8 asfreq_MtoS(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoS(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM QUARTERLY ***************

static void QtoD_ym(i8 ordinal, i8 *y, i8 *m, asfreq_info *af_info) {
    *y = floordiv(ordinal, 4) + BASE_YEAR;
    *m = mod_compat(ordinal, 4) * 3 + 1;

    if (af_info->from_q_year_end != 12) {
        *m += af_info->from_q_year_end;

        if (*m > 12)
            *m -= 12;
        else
            *y -= 1;
    }
}

static i8 asfreq_QtoD(i8 ordinal, const char* relation, asfreq_info *af_info) {

    i8 absdate, y, m;

    if (!strcmp(relation, "S")) {
        QtoD_ym(ordinal, &y, &m, af_info);

        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE)
            return INT_ERR_CODE;

        return absdate - ORD_OFFSET;
    } else {
        QtoD_ym(ordinal + 1, &y, &m, af_info);

        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE)
            return INT_ERR_CODE;

        return absdate - 1 - ORD_OFFSET;
    }
}

static i8 asfreq_QtoQ(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoQ(asfreq_QtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_QtoA(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_QtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_QtoM(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }

static i8 asfreq_QtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_QtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_QtoB(i8 ordinal, const char* relation, asfreq_info *af_info) {

    struct date_info dinfo;
    i8 qtod = asfreq_QtoD(ordinal, relation, af_info) + ORD_OFFSET;

    if (dInfoCalc_SetFromAbsDate(&dinfo, qtod, GREGORIAN))
        return INT_ERR_CODE;

    i8 (*f)(i8 abdate, i8 day_of_week) = NULL;
    f = !strcmp(relation, "S") ? DtoB_WeekendToMonday : DtoB_WeekendToFriday;

    return f(dinfo.absdate, dinfo.day_of_week);
}


static i8 asfreq_QtoH(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoH(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static i8 asfreq_QtoT(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoT(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static i8 asfreq_QtoS(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoS(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }


//************ FROM ANNUAL ***************

static i8 asfreq_AtoD(i8 ordinal, const char* relation, asfreq_info *af_info) {
    i8 absdate, final_adj, year;
    i8 month = af_info->from_a_year_end % 12, ord = ordinal;

	// start from 1970
	ord += BASE_YEAR;

    if (month == 0)
        month = 1;
    else
        ++month;


    if (!strcmp(relation, "S")) {
        year = af_info->from_a_year_end == 12 ? ord : ord - 1;
        final_adj = 0;
    } else {
        if (af_info->from_a_year_end == 12)
            year = ord + 1;
        else
            year = ord;

        final_adj = -1;
    }

    absdate = absdate_from_ymd(year, month, 1);

    if (absdate == INT_ERR_CODE)
        return INT_ERR_CODE;

    return absdate + final_adj - ORD_OFFSET;
}

static i8 asfreq_AtoA(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoA(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_AtoQ(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoQ(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_AtoM(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoM(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_AtoW(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoW(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static i8 asfreq_AtoB(i8 ordinal, const char* relation, asfreq_info *af_info) {
    struct date_info dinfo;

    i8 atob = asfreq_AtoD(ordinal, relation, af_info) + ORD_OFFSET;

    if (dInfoCalc_SetFromAbsDate(&dinfo, atob, GREGORIAN))
        return INT_ERR_CODE;

    i8 (*f)(i8 date, i8 day_of_week) = NULL;
    f = !strcmp(relation, "S") ? DtoB_WeekendToMonday : DtoB_WeekendToFriday;

    return f(dinfo.absdate, dinfo.day_of_week);
}

static i8 asfreq_AtoH(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoH(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static i8 asfreq_AtoT(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoT(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static i8 asfreq_AtoS(i8 ordinal, const char* relation, asfreq_info *af_info)
{ return asfreq_DtoS(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }

static i8 nofunc(i8 ordinal, const char* relation, asfreq_info *af_info) { return INT_ERR_CODE; }
static i8 no_op(i8 ordinal, const char* relation, asfreq_info *af_info) { return ordinal; }

// end of frequency specific conversion routines

static i8 get_freq_group(i8 freq) { return (freq / 1000) * 1000; }

static i8 calc_a_year_end(i8 freq, i8 group) {
    i8 result = (freq - group) % 12;
    return !result ? 12 : result;
}

static i8 calc_week_end(i8 freq, i8 group) {
    return freq - group;
}

void get_asfreq_info(i8 fromFreq, i8 toFreq, asfreq_info *af_info) {
    i8 fromGroup = get_freq_group(fromFreq);
    i8 toGroup = get_freq_group(toFreq);

    switch(fromGroup)
    {
    case FR_WK:
        af_info->from_week_end = calc_week_end(fromFreq, fromGroup);
        break;
    case FR_ANN:
        af_info->from_a_year_end = calc_a_year_end(fromFreq, fromGroup);
        break;
    case FR_QTR:
        af_info->from_q_year_end = calc_a_year_end(fromFreq, fromGroup);
        break;
    }

    switch(toGroup)
    {
    case FR_WK:
        af_info->to_week_end = calc_week_end(toFreq, toGroup);
        break;
    case FR_ANN:
        af_info->to_a_year_end = calc_a_year_end(toFreq, toGroup);
        break;
    case FR_QTR:
        af_info->to_q_year_end = calc_a_year_end(toFreq, toGroup);
        break;
    }
}


static i8 asfreq_UtoS(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return ordinal / US_PER_SECOND;
}

static i8 asfreq_UtoD(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_StoD(asfreq_UtoS(ordinal, relation, &NULL_AF_INFO), relation,
                       &NULL_AF_INFO);
}

static i8 asfreq_UtoA(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_UtoQ(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_UtoM(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_UtoW(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoW(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_UtoB(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoB(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       af_info);
}

static i8 asfreq_UtoH(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoH(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO),
                       relation, &NULL_AF_INFO);

}

static i8 asfreq_UtoT(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoT(asfreq_UtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       &NULL_AF_INFO);
}


static i8 asfreq_StoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_StoD(ordinal, relation, af_info), relation, &NULL_AF_INFO);
}

static i8 asfreq_AtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO);
}

static i8 asfreq_QtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO);
}

static i8 asfreq_MtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_MtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO);
}

static i8 asfreq_WtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO);
}

static i8 asfreq_BtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_BtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO);
}

static i8 asfreq_HtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    /* return ordinal * US_PER_HOUR; */
    return asfreq_DtoU(asfreq_HtoD(ordinal, relation, af_info), relation,
                       af_info);
}

static i8 asfreq_TtoU(i8 ordinal, const char* relation, asfreq_info *af_info) {
    return asfreq_DtoU(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation,
                       &NULL_AF_INFO);
}


freq_conv_func get_asfreq_func(i8 fromFreq, i8 toFreq)
{
    i8 fromGroup = get_freq_group(fromFreq);
    i8 toGroup = get_freq_group(toFreq);

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
        case FR_USEC: return &asfreq_AtoU;
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
        case FR_USEC: return &asfreq_QtoU;
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
        case FR_USEC: return &asfreq_MtoU;
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
        case FR_USEC: return &asfreq_WtoU;
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
        case FR_USEC: return &asfreq_BtoU;
        default: return &nofunc;
        }

    case FR_DAY:
        switch(toGroup)
        {
        case FR_ANN: return &asfreq_DtoA;
        case FR_QTR: return &asfreq_DtoQ;
        case FR_MTH: return &asfreq_DtoM;
        case FR_WK: return &asfreq_DtoW;
        case FR_BUS: return &asfreq_DtoB;
        case FR_DAY: return &asfreq_DtoD;
        case FR_HR: return &asfreq_DtoH;
        case FR_MIN: return &asfreq_DtoT;
        case FR_SEC: return &asfreq_DtoS;
        case FR_USEC: return &asfreq_DtoU;
        default: return &nofunc;
        }

    case FR_HR:
        switch(toGroup)
        {
        case FR_ANN: return &asfreq_HtoA;
        case FR_QTR: return &asfreq_HtoQ;
        case FR_MTH: return &asfreq_HtoM;
        case FR_WK: return &asfreq_HtoW;
        case FR_BUS: return &asfreq_HtoB;
        case FR_DAY: return &asfreq_HtoD;
        case FR_HR: return &no_op;
        case FR_MIN: return &asfreq_HtoT;
        case FR_SEC: return &asfreq_HtoS;
        case FR_USEC: return &asfreq_HtoU;
        default: return &nofunc;
        }

    case FR_MIN:
        switch(toGroup)
        {
        case FR_ANN: return &asfreq_TtoA;
        case FR_QTR: return &asfreq_TtoQ;
        case FR_MTH: return &asfreq_TtoM;
        case FR_WK: return &asfreq_TtoW;
        case FR_BUS: return &asfreq_TtoB;
        case FR_DAY: return &asfreq_TtoD;
        case FR_HR: return &asfreq_TtoH;
        case FR_MIN: return &no_op;
        case FR_SEC: return &asfreq_TtoS;
        case FR_USEC: return &asfreq_TtoU;
        default: return &nofunc;
        }

    case FR_SEC:
        switch(toGroup)
        {
        case FR_ANN: return &asfreq_StoA;
        case FR_QTR: return &asfreq_StoQ;
        case FR_MTH: return &asfreq_StoM;
        case FR_WK: return &asfreq_StoW;
        case FR_BUS: return &asfreq_StoB;
        case FR_DAY: return &asfreq_StoD;
        case FR_HR: return &asfreq_StoH;
        case FR_MIN: return &asfreq_StoT;
        case FR_SEC: return &no_op;
        case FR_USEC: return &asfreq_StoU;
        default: return &nofunc;
        }

    case FR_USEC:
        switch(toGroup)
        {
        case FR_ANN: return &asfreq_UtoA;
        case FR_QTR: return &asfreq_UtoQ;
        case FR_MTH: return &asfreq_UtoM;
        case FR_WK: return &asfreq_UtoW;
        case FR_BUS: return &asfreq_UtoB;
        case FR_DAY: return &asfreq_UtoD;
        case FR_HR: return &asfreq_UtoH;
        case FR_MIN: return &asfreq_UtoT;
        case FR_SEC: return &asfreq_UtoS;
        case FR_USEC: return &no_op;
        default: return &nofunc;
        }

    default: return &nofunc;
    }
}


static i8 get_abs_time(i8 freq, i8 daily_ord, i8 ordinal) {

    i8 start_ord, per_day, unit;

    switch (freq)
    {
    case FR_HR:
        per_day = 24L;
        unit = US_PER_HOUR;
        break;

    case FR_MIN:
        per_day = 24L * 60L;
        unit = US_PER_MINUTE;
        break;

    case FR_SEC:
        per_day = 24L * 60L * 60L;
        unit = US_PER_SECOND;
        break;

    case FR_USEC:
        per_day = US_PER_DAY;
        unit = 1L;
        break;

    default:
        return 0L; // 24*60*60 - 1;
    }

    start_ord = asfreq_DtoHIGHFREQ(daily_ord, "S", per_day);
	return unit * (ordinal - start_ord);
}

/* Sets the time part of the DateTime object. */
static
i8 dInfoCalc_SetFromAbsTime(struct date_info *dinfo, i8 abstime)
{
    i8 hour = abstime / US_PER_HOUR;
    i8 minute = (abstime % US_PER_HOUR) / US_PER_MINUTE;
    i8 second = (abstime - (hour * US_PER_HOUR + minute * US_PER_MINUTE)) /
        US_PER_SECOND;
    i8 microsecond = abstime - (hour * US_PER_HOUR + minute * US_PER_MINUTE +
                                second * US_PER_SECOND);

    dinfo->hour = hour;
    dinfo->minute = minute;
    dinfo->second = second;
    dinfo->microsecond = microsecond;

    dinfo->abstime = abstime;

    return 0;
}

/* Set the instance's value using the given date and time. calendar
   may be set to the flags: GREGORIAN, JULIAN to
   indicate the calendar to be used. */
static
i8 dInfoCalc_SetFromAbsDateTime(struct date_info *dinfo, i8 absdate,
                                i8 abstime, i8 calendar)
{

    /* Bounds check */
    Py_AssertWithArg(abstime >= 0 && abstime <= US_PER_DAY,
                     PyExc_ValueError,
                     "abstime out of range (0, 86400000000): %li",
                     abstime);

    /* Calculate the date */
    if (dInfoCalc_SetFromAbsDate(dinfo, absdate, calendar))
        goto onError;

    /* Calculate the time */
    if (dInfoCalc_SetFromAbsTime(dinfo, abstime))
        goto onError;

    return 0;

onError:
    return INT_ERR_CODE;
}

/* ------------------------------------------------------------------
 * New pandas API-helper code, to expose to cython
 * ------------------------------------------------------------------*/

i8 asfreq(i8 period_ordinal, i8 freq1, i8 freq2, const char* relation)
{
	freq_conv_func func = get_asfreq_func(freq1, freq2);

    asfreq_info finfo;
    get_asfreq_info(freq1, freq2, &finfo);

    i8 val = func(period_ordinal, relation, &finfo);

    if (val == INT_ERR_CODE) {
		goto onError;
	}

    return val;

onError:
    return INT_ERR_CODE;
}


/* generate an ordinal in period space */
i8 get_period_ordinal(i8 year, i8 month, i8 day, i8 hour, i8 minute, i8 second,
                      i8 microsecond, i8 freq)
{
    i8 absdays, delta, weeks, days, ordinal, day_adj, fmonth, mdiff;
    i8 freq_group = get_freq_group(freq);

    if (freq == FR_USEC) {
        absdays = absdate_from_ymd(year, month, day);
        delta = absdays - ORD_OFFSET;
        return delta * US_PER_DAY + hour * US_PER_HOUR +
            minute * US_PER_MINUTE + second * US_PER_SECOND + microsecond;
    }

    if (freq == FR_SEC) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        return delta*86400 + hour*3600 + minute*60 + second;
    }

    if (freq == FR_MIN) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        return delta*1440 + hour*60 + minute;
    }

    if (freq == FR_HR) {
        if ((absdays = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        delta = (absdays - ORD_OFFSET);
        return delta*24 + hour;
    }

    if (freq == FR_DAY)
    {
        return absdate_from_ymd(year, month, day) - ORD_OFFSET;
    }

    if (freq == FR_UND)
    {
        return absdate_from_ymd(year, month, day) - ORD_OFFSET;
    }

    if (freq == FR_BUS)
    {
        if ((days = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        weeks = days / 7;
        return (days - weeks * 2) - BDAY_OFFSET;
    }

    if (freq_group == FR_WK)
    {
        if ((ordinal = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        day_adj = freq - FR_WK;
        return (ordinal - (1 + day_adj)) / 7 + 1 - WEEK_OFFSET;
    }

    if (freq == FR_MTH)
    {
        return (year - BASE_YEAR) * 12 + month - 1;
    }

    if (freq_group == FR_QTR)
    {
        fmonth = freq - FR_QTR;
        if (fmonth == 0) fmonth = 12;

        mdiff = month - fmonth;
        if (mdiff < 0) mdiff += 12;
        if (month >= fmonth) mdiff += 12;

        return (year - BASE_YEAR) * 4 + (mdiff - 1) / 3;
    }

    if (freq_group == FR_ANN)
    {
        fmonth = freq - FR_ANN;
        if (fmonth == 0) fmonth = 12;
        if (month <= fmonth) {
            return year - BASE_YEAR;
        }
        else {
            return year - BASE_YEAR + 1;
        }
    }

    Py_Error(PyExc_RuntimeError, "Unable to generate frequency ordinal");

onError:
    return INT_ERR_CODE;
}

/*
  Returns the proleptic Gregorian ordinal of the date, as an i8eger.
  This corresponds to the number of days since Jan., 1st, 1AD.
  When the instance has a frequency less than daily, the proleptic date
  is calculated for the last day of the period.
*/

i8 get_python_ordinal(i8 period_ordinal, i8 freq)
{
    if (freq == FR_DAY)
        return period_ordinal + ORD_OFFSET;

    i8 (*toDaily)(i8, const char*, asfreq_info*) = get_asfreq_func(freq, FR_DAY);
    asfreq_info af_info;

    get_asfreq_info(freq, FR_DAY, &af_info);

    return toDaily(period_ordinal, "E", &af_info) + ORD_OFFSET;
}

char *str_replace(const char *s, const char *old, const char *new)
{
    char *ret = NULL;
    i8 i, count = 0;
    size_t newlen = strlen(new);
    size_t oldlen = strlen(old);

    for (i = 0; s[i] != '\0'; i++) {
        if (strstr(&s[i], old) == &s[i]) {
            ++count;
            i += oldlen - 1;
        }
    }

    ret = PyArray_malloc(i + 1 + count * (newlen - oldlen));
    if (ret == NULL)
        return (char*) PyErr_NoMemory();


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

char* c_strftime(struct date_info *tmp, char *fmt)
{
    struct tm c_date;
    struct date_info dinfo = *tmp;

    size_t result_len = strlen(fmt) + 50L;

    c_date.tm_sec = dinfo.second;
    c_date.tm_min = dinfo.minute;
    c_date.tm_hour = dinfo.hour;
    c_date.tm_mday = dinfo.day;
    c_date.tm_mon = dinfo.month - 1L;
    c_date.tm_year = dinfo.year - 1900L;
    c_date.tm_wday = (dinfo.day_of_week + 1L) % 7L;
    c_date.tm_yday = dinfo.day_of_year - 1L;
    c_date.tm_isdst = -1L;

    // this will be freed in cython
    char* result = (char*) malloc(result_len * sizeof(char));

    strftime(result, result_len, fmt, &c_date);

    // check to see if the format is sub second
    i8 is_subsec = strstr(fmt, ".") != NULL;

    // workaround the fact that strftime doesn't do sub second printing
    char* result2 = NULL;

    if (is_subsec) {
        result2 = (char*) malloc(result_len * sizeof(char));
        snprintf(result2, result_len, result, dinfo.microsecond);
        free(result);
        result = result2;
    }

    return result;
}


i8 get_yq(i8 ordinal, i8 freq, i8 *quarter, i8 *year)
{
    asfreq_info af_info;
    i8 qtr_freq, daily_ord;
    i8 (*toDaily)(i8, const char*, asfreq_info*) = get_asfreq_func(freq, FR_DAY);

    get_asfreq_info(freq, FR_DAY, &af_info);

    daily_ord = toDaily(ordinal, "E", &af_info);

    qtr_freq = get_freq_group(freq) == FR_QTR ? freq : FR_QTR;

    get_asfreq_info(FR_DAY, qtr_freq, &af_info);

    if (DtoQ_yq(daily_ord, &af_info, year, quarter) == INT_ERR_CODE)
        return -1;

    return 0;
}


static i8 _quarter_year(i8 ordinal, i8 freq, i8 *year, i8 *quarter)
{
    asfreq_info af_info;

    i8 ord = get_python_ordinal(ordinal, freq) - ORD_OFFSET;
    i8 qtr_freq = get_freq_group(freq) == FR_QTR ? freq : FR_QTR;

    get_asfreq_info(FR_DAY, qtr_freq, &af_info);

    if (DtoQ_yq(ord, &af_info, year, quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;

    if (qtr_freq % 1000 > 12)
        *year -= 1;

    return 0;
}


static i8 _ISOWeek(struct date_info *dinfo)
{
    i8 week;

    /* Estimate */
    week = (dinfo->day_of_year - 1) - dinfo->day_of_week + 3;

    if (week >= 0) {
        week /= 7;
        ++week;
    }


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
        if (31 - dinfo->day + dinfo->day_of_week < 3) {
            week = 1;
        }
    }

    return week;
}


i8 get_date_info(i8 ordinal, i8 freq, struct date_info *dinfo)
{
    i8 absdate = get_python_ordinal(ordinal, freq);
    i8 abstime = get_abs_time(freq, absdate - ORD_OFFSET, ordinal);

	if (abstime < 0) {
        // add a day's worth of us to time
		abstime += US_PER_DAY;

        // subtract a day from the date
		--absdate;
	}

    if (dInfoCalc_SetFromAbsDateTime(dinfo, absdate, abstime, GREGORIAN))
        return INT_ERR_CODE;

    return 0;
}

i8 pyear(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    get_date_info(ordinal, freq, &dinfo);
    return dinfo.year;
}

i8 pqyear(i8 ordinal, i8 freq)
{
    i8 year, quarter;
    if (_quarter_year(ordinal, freq, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return year;
}

i8 pquarter(i8 ordinal, i8 freq)
{
    i8 year, quarter;
    if (_quarter_year(ordinal, freq, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return quarter;
}

i8 pmonth(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.month;
}

i8 pday(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day;
}

i8 pweekday(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_week;
}

i8 pday_of_week(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_week;
}

i8 pday_of_year(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_year;
}

i8 pweek(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return _ISOWeek(&dinfo);
}

i8 phour(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.hour;
}

i8 pminute(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.minute;
}

i8 psecond(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.second;
}

i8 pmicrosecond(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.microsecond;
}

i8 pnanosecond(i8 ordinal, i8 freq)
{
    struct date_info dinfo;
    if (get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.nanosecond;
}
