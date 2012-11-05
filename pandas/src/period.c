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


static int mod_compat(int x, int m) {
  int result = x % m;
  if (result < 0) return result + m;
  return result;
}

static int floordiv(int x, int divisor) {
    if (x < 0) {
        if (mod_compat(x, divisor)) {
            return x / divisor - 1;
        }
        else return x / divisor;
    } else {
        return x / divisor;
    }
}

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
static int dInfoCalc_Leapyear(npy_int64 year, int calendar)
{
    if (calendar == GREGORIAN_CALENDAR) {
        return (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));
    } else {
        return (year % 4 == 0);
    }
}

/* Return the day of the week for the given absolute date. */
static int dInfoCalc_DayOfWeek(npy_int64 absdate)
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
static int dInfoCalc_YearOffset(npy_int64 year, int calendar)
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
		npy_int64 absdate;
        int yearoffset;

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
							 npy_int64 absdate, int calendar)
{
    register npy_int64 year;
    npy_int64 yearoffset;
    int leap,dayoffset;
    int *monthoffset;

    /* Approximate year */
    if (calendar == GREGORIAN_CALENDAR) {
        year = (npy_int64)(((double)absdate) / 365.2425);
    } else if (calendar == JULIAN_CALENDAR) {
        year = (npy_int64)(((double)absdate) / 365.25);
    } else {
        Py_Error(PyExc_ValueError, "unknown calendar");
    }

    if (absdate > 0) year++;

    /* Apply corrections to reach the correct year */
    while (1) {
        /* Calculate the year offset */
        yearoffset = dInfoCalc_YearOffset(year, calendar);
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

static npy_int64 DtoB_weekday(npy_int64 absdate) {
    return (((absdate) / 7) * 5) + (absdate) % 7 - BDAY_OFFSET;
}

static npy_int64 DtoB_WeekendToMonday(npy_int64 absdate, int day_of_week) {
    if (day_of_week > 4) {
        //change to Monday after weekend
        absdate += (7 - day_of_week);
    }
    return DtoB_weekday(absdate);
}

static npy_int64 DtoB_WeekendToFriday(npy_int64 absdate, int day_of_week) {
    if (day_of_week > 4) {
        //change to friday before weekend
        absdate -= (day_of_week - 4);
    }
    return DtoB_weekday(absdate);
}

static npy_int64 absdate_from_ymd(int y, int m, int d) {
    struct date_info tempDate;
    if (dInfoCalc_SetFromDateAndTime(&tempDate, y, m, d, 0, 0, 0, GREGORIAN_CALENDAR)) {
        return INT_ERR_CODE;
    }
    return tempDate.absdate;
}

//************ FROM DAILY ***************

static npy_int64 asfreq_DtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
								 GREGORIAN_CALENDAR)) return INT_ERR_CODE;
    if (dinfo.month > af_info->to_a_year_end) {
	  return (npy_int64)(dinfo.year + 1 - BASE_YEAR);
	}
    else {
	  return (npy_int64)(dinfo.year - BASE_YEAR);
	}
}

static npy_int64 DtoQ_yq(npy_int64 ordinal, asfreq_info *af_info,
					   int *year, int *quarter) {
    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
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


static npy_int64 asfreq_DtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    int year, quarter;

    if (DtoQ_yq(ordinal, af_info, &year, &quarter) == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    return (npy_int64)((year - BASE_YEAR) * 4 + quarter - 1);
}

static npy_int64 asfreq_DtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;
    return (npy_int64)((dinfo.year - BASE_YEAR) * 12 + dinfo.month - 1);
}

static npy_int64 asfreq_DtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return (ordinal + ORD_OFFSET - (1 + af_info->to_week_end))/7 + 1 - WEEK_OFFSET;
}

static npy_int64 asfreq_DtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
								 GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') {
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
    } else {
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
    }
}

// needed for getDateInfo function
static npy_int64 asfreq_DtoD(npy_int64 ordinal, char relation, asfreq_info *af_info) { return ordinal; }

static npy_int64 asfreq_DtoHIGHFREQ(npy_int64 ordinal, char relation, npy_int64 per_day) {
	if (relation == 'S') {
	  return ordinal * per_day;
	}
	else {
	  return (ordinal+ 1) * per_day - 1;
	}
}

static npy_int64 asfreq_DtoH(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(ordinal, relation, 24); }
static npy_int64 asfreq_DtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(ordinal, relation, 24*60); }
static npy_int64 asfreq_DtoS(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoHIGHFREQ(ordinal, relation, 24*60*60); }

//************ FROM SECONDLY ***************

static npy_int64 asfreq_StoD(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return (ordinal)/(60*60*24); }

static npy_int64 asfreq_StoA(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_StoQ(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_StoM(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_StoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_StoB(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_StoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }


static npy_int64 asfreq_StoT(npy_int64 ordinal, char relation, asfreq_info *af_info) {
	return ordinal / 60;
}

static npy_int64 asfreq_StoH(npy_int64 ordinal, char relation, asfreq_info *af_info) {
	return ordinal / (60*60);
}

//************ FROM MINUTELY ***************

static npy_int64 asfreq_TtoD(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return (ordinal)/(60*24); }

static npy_int64 asfreq_TtoA(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static npy_int64 asfreq_TtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static npy_int64 asfreq_TtoM(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_TtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static npy_int64 asfreq_TtoB(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_TtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_TtoH(npy_int64 ordinal, char relation, asfreq_info *af_info) {
	return ordinal / 60;
}

static npy_int64 asfreq_TtoS(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    if (relation == 'S') {
		return ordinal*60; }
    else                 {
		return ordinal*60 + 59;
	}
}

//************ FROM HOURLY ***************

static npy_int64 asfreq_HtoD(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return ordinal / 24; }
static npy_int64 asfreq_HtoA(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static npy_int64 asfreq_HtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static npy_int64 asfreq_HtoM(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_HtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }
static npy_int64 asfreq_HtoB(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoB(asfreq_HtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

// calculation works out the same as TtoS, so we just call that function for HtoT
static npy_int64 asfreq_HtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_TtoS(ordinal, relation, &NULL_AF_INFO); }

static npy_int64 asfreq_HtoS(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    if (relation == 'S') {
		return ordinal*60*60;
	}
    else {
		return (ordinal + 1)*60*60 - 1;
	}
}

//************ FROM BUSINESS ***************

static npy_int64 asfreq_BtoD(npy_int64 ordinal, char relation, asfreq_info *af_info)
    {
		ordinal += BDAY_OFFSET;
		return (((ordinal - 1) / 5) * 7 +
				mod_compat(ordinal - 1, 5) + 1 - ORD_OFFSET);
	}

static npy_int64 asfreq_BtoA(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_BtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_BtoM(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_BtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_BtoH(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_BtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_BtoS(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_BtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM WEEKLY ***************

static npy_int64 asfreq_WtoD(npy_int64 ordinal, char relation, asfreq_info *af_info) {
	ordinal += WEEK_OFFSET;
    if (relation == 'S') {
	  return ordinal * 7 - 6 + af_info->from_week_end - ORD_OFFSET;
	}
    else {
	  return ordinal * 7 + af_info->from_week_end - ORD_OFFSET;
	}
}

static npy_int64 asfreq_WtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_WtoD(ordinal, 'E', af_info), relation, af_info); }
static npy_int64 asfreq_WtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_WtoD(ordinal, 'E', af_info), relation, af_info); }
static npy_int64 asfreq_WtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_WtoD(ordinal, 'E', af_info), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_WtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_WtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_WtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
								 asfreq_WtoD(ordinal, relation, af_info) + ORD_OFFSET,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') {
		return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
	}
    else {
		return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
	}
}

static npy_int64 asfreq_WtoH(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_WtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_WtoS(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_WtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }

//************ FROM MONTHLY ***************
static void MtoD_ym(npy_int64 ordinal, int *y, int *m) {
    *y = floordiv(ordinal, 12) + BASE_YEAR;
    *m = mod_compat(ordinal, 12) + 1;
}


static npy_int64 asfreq_MtoD(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    npy_int64 absdate;
    int y, m;

    if (relation == 'S') {
        MtoD_ym(ordinal, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate - ORD_OFFSET;
    } else {
        MtoD_ym(ordinal + 1, &y, &m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate - 1 - ORD_OFFSET;
    }
}

static npy_int64 asfreq_MtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_MtoD(ordinal, 'E', &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_MtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoQ(asfreq_MtoD(ordinal, 'E', &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_MtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, af_info); }

static npy_int64 asfreq_MtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
								 asfreq_MtoD(ordinal, relation, &NULL_AF_INFO) + ORD_OFFSET,
								 GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static npy_int64 asfreq_MtoH(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_MtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_MtoS(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_MtoD(ordinal, relation, &NULL_AF_INFO), relation, &NULL_AF_INFO); }

//************ FROM QUARTERLY ***************

static void QtoD_ym(npy_int64 ordinal, int *y, int *m, asfreq_info *af_info) {
    *y = floordiv(ordinal, 4) + BASE_YEAR;
    *m = mod_compat(ordinal, 4) * 3 + 1;

    if (af_info->from_q_year_end != 12) {
        *m += af_info->from_q_year_end;
        if (*m > 12) { *m -= 12; }
        else { *y -= 1; }
    }
}

static npy_int64 asfreq_QtoD(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    npy_int64 absdate;
    int y, m;

    if (relation == 'S') {
        QtoD_ym(ordinal, &y, &m, af_info);
		// printf("ordinal: %d, year: %d, month: %d\n", (int) ordinal, y, m);
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate - ORD_OFFSET;
    } else {
        QtoD_ym(ordinal+1, &y, &m, af_info);
		/* printf("ordinal: %d, year: %d, month: %d\n", (int) ordinal, y, m); */
        if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
        return absdate - 1 - ORD_OFFSET;
    }
}

static npy_int64 asfreq_QtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_QtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_QtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoA(asfreq_QtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_QtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return asfreq_DtoM(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }

static npy_int64 asfreq_QtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_QtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_QtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
								 asfreq_QtoD(ordinal, relation, af_info) + ORD_OFFSET,
								 GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}


static npy_int64 asfreq_QtoH(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_QtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_QtoS(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_QtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }


//************ FROM ANNUAL ***************

static npy_int64 asfreq_AtoD(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    npy_int64 absdate, final_adj;
	int year;
    int month = (af_info->from_a_year_end) % 12;

	// start from 1970
	ordinal += BASE_YEAR;

    if (month == 0) { month = 1; }
    else { month += 1; }

    if (relation == 'S') {
        if (af_info->from_a_year_end == 12) {year = ordinal;}
        else {year = ordinal - 1;}
        final_adj = 0;
    } else {
        if (af_info->from_a_year_end == 12) {year = ordinal+1;}
        else {year = ordinal;}
        final_adj = -1;
    }
    absdate = absdate_from_ymd(year, month, 1);
    if (absdate  == INT_ERR_CODE) {
	  return INT_ERR_CODE;
	}
    return absdate + final_adj - ORD_OFFSET;
}

static npy_int64 asfreq_AtoA(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoA(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_AtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoQ(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_AtoM(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoM(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_AtoW(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoW(asfreq_AtoD(ordinal, relation, af_info), relation, af_info); }

static npy_int64 asfreq_AtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
								 asfreq_AtoD(ordinal, relation, af_info) + ORD_OFFSET,
                    GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static npy_int64 asfreq_AtoH(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoH(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_AtoT(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoT(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }
static npy_int64 asfreq_AtoS(npy_int64 ordinal, char relation, asfreq_info *af_info)
    { return asfreq_DtoS(asfreq_AtoD(ordinal, relation, af_info), relation, &NULL_AF_INFO); }

static npy_int64 nofunc(npy_int64 ordinal, char relation, asfreq_info *af_info) { return INT_ERR_CODE; }
static npy_int64 no_op(npy_int64 ordinal, char relation, asfreq_info *af_info) { return ordinal; }

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


freq_conv_func get_asfreq_func(int fromFreq, int toFreq)
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
                case FR_BUS: return &asfreq_DtoB;
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
                case FR_BUS: return &asfreq_HtoB;
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
                case FR_BUS: return &asfreq_TtoB;
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
                case FR_BUS: return &asfreq_StoB;
                case FR_DAY: return &asfreq_StoD;
                case FR_HR: return &asfreq_StoH;
                case FR_MIN: return &asfreq_StoT;
                case FR_SEC: return &no_op;
                default: return &nofunc;
            }
        default: return &nofunc;
    }
}

double get_abs_time(int freq, npy_int64 daily_ord, npy_int64 ordinal) {

    npy_int64 start_ord, per_day, unit;
    switch(freq)
    {
        case FR_HR:
            per_day = 24;
			unit = 60 * 60;
            break;
        case FR_MIN:
            per_day = 24*60;
			unit = 60;
            break;
        case FR_SEC:
            per_day = 24*60*60;
			unit = 1;
            break;
        default:
		  return 0; // 24*60*60 - 1;
    }

    start_ord = asfreq_DtoHIGHFREQ(daily_ord, 'S', per_day);
	/* printf("start_ord: %d\n", start_ord); */
	return (double) ( unit * (ordinal - start_ord));
	/* if (ordinal >= 0) { */
	/* } */
	/* else { */
	/* 	return (double) (unit * mod_compat(ordinal - start_ord, per_day)); */
	/* } */
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
                  npy_int64 absdate,
                  double abstime,
                  int calendar)
{

    /* Bounds check */
    Py_AssertWithArg(abstime >= 0.0 && abstime <= SECONDS_PER_DAY,
             PyExc_ValueError,
             "abstime out of range (0.0 - 86400.0): %f",
             abstime);

    /* Calculate the date */
    if (dInfoCalc_SetFromAbsDate(dinfo, absdate, calendar)) goto onError;

    /* Calculate the time */
    if (dInfoCalc_SetFromAbsTime(dinfo, abstime)) goto onError;

    return 0;
 onError:
    return INT_ERR_CODE;
}

/* ------------------------------------------------------------------
 * New pandas API-helper code, to expose to cython
 * ------------------------------------------------------------------*/

npy_int64 asfreq(npy_int64 period_ordinal, int freq1, int freq2, char relation)
{
    npy_int64 val;
    freq_conv_func func;
    asfreq_info finfo;

	func = get_asfreq_func(freq1, freq2);
    get_asfreq_info(freq1, freq2, &finfo);

    val = (*func)(period_ordinal, relation, &finfo);

    if (val == INT_ERR_CODE) {
        // Py_Error(PyExc_ValueError, "Unable to convert to desired frequency.");
		goto onError;
	}
    return val;
onError:
    return INT_ERR_CODE;
}


/* generate an ordinal in period space */
npy_int64 get_period_ordinal(int year, int month, int day,
                      int hour, int minute, int second,
                      int freq)
{
    npy_int64 absdays, delta;
    npy_int64 weeks, days;
    npy_int64 ordinal, day_adj;
    int freq_group, fmonth, mdiff;
    freq_group = get_freq_group(freq);

    if (freq == FR_SEC) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        return (npy_int64)(delta*86400 + hour*3600 + minute*60 + second);
    }

    if (freq == FR_MIN) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        return (npy_int64)(delta*1440 + hour*60 + minute);
    }

    if (freq == FR_HR) {
        if ((absdays = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        delta = (absdays - ORD_OFFSET);
        return (npy_int64)(delta*24 + hour);
    }

    if (freq == FR_DAY)
    {
	  return (npy_int64) (absdate_from_ymd(year, month, day) - ORD_OFFSET);
    }

    if (freq == FR_UND)
    {
	  return (npy_int64) (absdate_from_ymd(year, month, day) - ORD_OFFSET);
    }

    if (freq == FR_BUS)
    {
        if((days = absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
        {
            goto onError;
        }
        weeks = days / 7;
        return (npy_int64)(days - weeks * 2) - BDAY_OFFSET;
    }

    if (freq_group == FR_WK)
    {
        if((ordinal = (npy_int64)absdate_from_ymd(year, month, day)) == INT_ERR_CODE)
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
    Returns the proleptic Gregorian ordinal of the date, as an integer.
    This corresponds to the number of days since Jan., 1st, 1AD.
    When the instance has a frequency less than daily, the proleptic date
    is calculated for the last day of the period.
*/

npy_int64 get_python_ordinal(npy_int64 period_ordinal, int freq)
{
    asfreq_info af_info;
    npy_int64 (*toDaily)(npy_int64, char, asfreq_info*);

    if (freq == FR_DAY)
        return period_ordinal + ORD_OFFSET;

    toDaily = get_asfreq_func(freq, FR_DAY);
    get_asfreq_info(freq, FR_DAY, &af_info);
    return toDaily(period_ordinal, 'E', &af_info) + ORD_OFFSET;
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

char* c_strftime(struct date_info *tmp, char *fmt) {
    struct tm c_date;
    char* result;
    struct date_info dinfo = *tmp;
    int result_len = strlen(fmt) + 50;

    c_date.tm_sec = (int)dinfo.second;
    c_date.tm_min = dinfo.minute;
    c_date.tm_hour = dinfo.hour;
    c_date.tm_mday = dinfo.day;
    c_date.tm_mon = dinfo.month - 1;
    c_date.tm_year = dinfo.year - 1900;
    c_date.tm_wday = (dinfo.day_of_week + 1) % 7;
    c_date.tm_yday = dinfo.day_of_year - 1;
    c_date.tm_isdst = -1;

    result = malloc(result_len * sizeof(char));

    strftime(result, result_len, fmt, &c_date);

    return result;
}

int get_yq(npy_int64 ordinal, int freq, int *quarter, int *year) {
    asfreq_info af_info;
    int qtr_freq;
    npy_int64 daily_ord;
    npy_int64 (*toDaily)(npy_int64, char, asfreq_info*) = NULL;

    toDaily = get_asfreq_func(freq, FR_DAY);
    get_asfreq_info(freq, FR_DAY, &af_info);

    daily_ord = toDaily(ordinal, 'E', &af_info);

    if (get_freq_group(freq) == FR_QTR) {
        qtr_freq = freq;
    } else { qtr_freq = FR_QTR; }
    get_asfreq_info(FR_DAY, qtr_freq, &af_info);

    if(DtoQ_yq(daily_ord, &af_info, year, quarter) == INT_ERR_CODE)
        return -1;

    return 0;
}





static int _quarter_year(npy_int64 ordinal, int freq, int *year, int *quarter) {
    asfreq_info af_info;
    int qtr_freq;

    ordinal = get_python_ordinal(ordinal, freq) - ORD_OFFSET;

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

int get_date_info(npy_int64 ordinal, int freq, struct date_info *dinfo)
{
    npy_int64 absdate = get_python_ordinal(ordinal, freq);
	/* printf("freq: %d, absdate: %d\n", freq, (int) absdate); */
    double abstime = get_abs_time(freq, absdate - ORD_OFFSET, ordinal);
	if (abstime < 0) {
		abstime += 86400;
		absdate -= 1;
	}

    if(dInfoCalc_SetFromAbsDateTime(dinfo, absdate,
									abstime, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;

    return 0;
}

int pyear(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    get_date_info(ordinal, freq, &dinfo);
    return dinfo.year;
}

int pqyear(npy_int64 ordinal, int freq) {
    int year, quarter;
    if( _quarter_year(ordinal, freq, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return year;
}

int pquarter(npy_int64 ordinal, int freq) {
    int year, quarter;
    if(_quarter_year(ordinal, freq, &year, &quarter) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return quarter;
}

int pmonth(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.month;
}

int pday(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day;
}

int pweekday(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_week;
}

int pday_of_week(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_week;
}

int pday_of_year(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.day_of_year;
}

int pweek(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return _ISOWeek(&dinfo);
}

int phour(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.hour;
}

int pminute(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return dinfo.minute;
}

int psecond(npy_int64 ordinal, int freq) {
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    return (int)dinfo.second;
}
