#include "period_helper.h"


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
        if (yearoffset == INT_ERR_CODE) goto onError;

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
        if (yearoffset == INT_ERR_CODE) goto onError;

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

static int daytime_conversion_factors[][2] = {
    { FR_DAY, 1 },
    { FR_HR,  24 },
    { FR_MIN, 60 },
    { FR_SEC, 60 },
    { FR_MS,  1000 },
    { FR_US,  1000 },
    { FR_NS,  1000 },
    { 0, 0 }
};

static npy_int64** daytime_conversion_factor_matrix = NULL;

PANDAS_INLINE int max_value(int a, int b) {
    return a > b ? a : b;
}

PANDAS_INLINE int min_value(int a, int b) {
    return a < b ? a : b;
}

PANDAS_INLINE int get_freq_group(int freq) {
    return (freq/1000)*1000;
}

PANDAS_INLINE int get_freq_group_index(int freq) {
    return freq/1000;
}

static int calc_conversion_factors_matrix_size(void) {
    int matrix_size = 0;
    int index;
    for (index=0;; index++) {
        int period_value = get_freq_group_index(daytime_conversion_factors[index][0]);
        if (period_value == 0) {
            break;
        }
        matrix_size = max_value(matrix_size, period_value);
    }
    return matrix_size + 1;
}

static void alloc_conversion_factors_matrix(int matrix_size) {
    int row_index;
    int column_index;
	daytime_conversion_factor_matrix = malloc(matrix_size * sizeof(**daytime_conversion_factor_matrix));
    for (row_index = 0; row_index < matrix_size; row_index++) {
        daytime_conversion_factor_matrix[row_index] = malloc(matrix_size * sizeof(**daytime_conversion_factor_matrix));
        for (column_index = 0; column_index < matrix_size; column_index++) {
            daytime_conversion_factor_matrix[row_index][column_index] = 0;
        }
    }
}

static npy_int64 calculate_conversion_factor(int start_value, int end_value) {
    npy_int64 conversion_factor = 0;
    int index;
    for (index=0;; index++) {
        int freq_group = daytime_conversion_factors[index][0];

        if (freq_group == 0) {
            conversion_factor = 0;
            break;
        }

        if (freq_group == start_value) {
            conversion_factor = 1;
        } else {
            conversion_factor *= daytime_conversion_factors[index][1];
        }

        if (freq_group == end_value) {
            break;
        }
    }
    return conversion_factor;
}

static void populate_conversion_factors_matrix(void) {
    int row_index_index;
	int row_value, row_index;
    int column_index_index;
	int column_value, column_index;

	for (row_index_index = 0;; row_index_index++) {
        row_value = daytime_conversion_factors[row_index_index][0];
        if (row_value == 0) {
            break;
        }
        row_index = get_freq_group_index(row_value);
        for (column_index_index = row_index_index;; column_index_index++) {
            column_value = daytime_conversion_factors[column_index_index][0];
            if (column_value == 0) {
                break;
            }
            column_index = get_freq_group_index(column_value);

            daytime_conversion_factor_matrix[row_index][column_index] = calculate_conversion_factor(row_value, column_value);
        }
    }
}

void initialize_daytime_conversion_factor_matrix() {
    if (daytime_conversion_factor_matrix == NULL) {
        int matrix_size = calc_conversion_factors_matrix_size();
        alloc_conversion_factors_matrix(matrix_size);
        populate_conversion_factors_matrix();
    }
}

PANDAS_INLINE npy_int64 get_daytime_conversion_factor(int from_index, int to_index)
{
    return daytime_conversion_factor_matrix[min_value(from_index, to_index)][max_value(from_index, to_index)];
}

PANDAS_INLINE npy_int64 upsample_daytime(npy_int64 ordinal, asfreq_info *af_info, int atEnd)
{
    if (atEnd) {
        return (ordinal + 1) * af_info->intraday_conversion_factor - 1;
    } else {
        return ordinal * af_info->intraday_conversion_factor;
    }
}

PANDAS_INLINE npy_int64 downsample_daytime(npy_int64 ordinal, asfreq_info *af_info, int atEnd)
{
    return ordinal / (af_info->intraday_conversion_factor);
}

PANDAS_INLINE npy_int64 transform_via_day(npy_int64 ordinal, char relation, asfreq_info *af_info, freq_conv_func first_func, freq_conv_func second_func) {
    //printf("transform_via_day(%ld, %ld, %d)\n", ordinal, af_info->intraday_conversion_factor, af_info->intraday_conversion_upsample);
	npy_int64 result;

    result = (*first_func)(ordinal, relation, af_info);
    result = (*second_func)(result, relation, af_info);

    return result;
}

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

static npy_int64 asfreq_DTtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    struct date_info dinfo;
    ordinal = downsample_daytime(ordinal, af_info, 0);
    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;
    if (dinfo.month > af_info->to_a_year_end) {
        return (npy_int64)(dinfo.year + 1 - BASE_YEAR);
    }
    else {
        return (npy_int64)(dinfo.year - BASE_YEAR);
    }
}

static npy_int64 DtoQ_yq(npy_int64 ordinal, asfreq_info *af_info, int *year, int *quarter) {
    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN_CALENDAR)) 
        return INT_ERR_CODE;
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

static npy_int64 asfreq_DTtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    int year, quarter;

    ordinal = downsample_daytime(ordinal, af_info, 0);

    if (DtoQ_yq(ordinal, af_info, &year, &quarter) == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    return (npy_int64)((year - BASE_YEAR) * 4 + quarter - 1);
}

static npy_int64 asfreq_DTtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    struct date_info dinfo;

    ordinal = downsample_daytime(ordinal, af_info, 0);

    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;
    return (npy_int64)((dinfo.year - BASE_YEAR) * 12 + dinfo.month - 1);
}

static npy_int64 asfreq_DTtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    ordinal = downsample_daytime(ordinal, af_info, 0);
    return (ordinal + ORD_OFFSET - (1 + af_info->to_week_end))/7 + 1 - WEEK_OFFSET;
}

static npy_int64 asfreq_DTtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    struct date_info dinfo;

	ordinal = downsample_daytime(ordinal, af_info, 0);

    if (dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET, GREGORIAN_CALENDAR))
        return INT_ERR_CODE;

    if (relation == 'S') {
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
    } else {
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
    }
}

// all intra day calculations are now done within one function
static npy_int64 asfreq_DownsampleWithinDay(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return downsample_daytime(ordinal, af_info, relation == 'E');
}

static npy_int64 asfreq_UpsampleWithinDay(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return upsample_daytime(ordinal, af_info, relation == 'E');
}
//************ FROM BUSINESS ***************

static npy_int64 asfreq_BtoDT(npy_int64 ordinal, char relation, asfreq_info *af_info)
{
    ordinal += BDAY_OFFSET;
    ordinal = (((ordinal - 1) / 5) * 7 +
            mod_compat(ordinal - 1, 5) + 1 - ORD_OFFSET);

    return upsample_daytime(ordinal, af_info, relation != 'S');
}

static npy_int64 asfreq_BtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_BtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_BtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_BtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT, asfreq_DTtoW);
}

//************ FROM WEEKLY ***************

static npy_int64 asfreq_WtoDT(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    ordinal += WEEK_OFFSET;
    if (relation != 'S') {
        ordinal += 1;
    }

    ordinal = ordinal * 7 - 6 + af_info->from_week_end - ORD_OFFSET;

    if (relation != 'S') {
        ordinal -= 1;
    }

    return upsample_daytime(ordinal, af_info, relation != 'S');
}

static npy_int64 asfreq_WtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_WtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_WtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_WtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_WtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
                asfreq_WtoDT(ordinal, relation, af_info) + ORD_OFFSET,
                GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') {
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week);
    }
    else {
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week);
    }
}

//************ FROM MONTHLY ***************
static void MtoD_ym(npy_int64 ordinal, int *y, int *m) {
    *y = floordiv(ordinal, 12) + BASE_YEAR;
    *m = mod_compat(ordinal, 12) + 1;
}


static npy_int64 asfreq_MtoDT(npy_int64 ordinal, char relation, asfreq_info* af_info) {
    npy_int64 absdate;
    int y, m;

    if (relation == 'E') {
      ordinal += 1;
    }
    MtoD_ym(ordinal, &y, &m);
    if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;
    ordinal = absdate - ORD_OFFSET;

    if (relation == 'E') {
      ordinal -= 1;
    }

    return upsample_daytime(ordinal, af_info, relation != 'S');
}

static npy_int64 asfreq_MtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_MtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_MtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_MtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_MtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_MtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_MtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    struct date_info dinfo;
    
    if (dInfoCalc_SetFromAbsDate(&dinfo,
                asfreq_MtoDT(ordinal, relation, af_info) + ORD_OFFSET,
                GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

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

static npy_int64 asfreq_QtoDT(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    npy_int64 absdate;
    int y, m;

    if (relation == 'E') {
      ordinal += 1;
    }

    QtoD_ym(ordinal, &y, &m, af_info);

    if ((absdate = absdate_from_ymd(y, m, 1)) == INT_ERR_CODE) return INT_ERR_CODE;

    if (relation == 'E') {
      absdate -= 1;
    }

    return upsample_daytime(absdate - ORD_OFFSET, af_info, relation != 'S');
}

static npy_int64 asfreq_QtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_QtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_QtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_QtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_QtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
                asfreq_QtoDT(ordinal, relation, af_info) + ORD_OFFSET,
                GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}


//************ FROM ANNUAL ***************

static npy_int64 asfreq_AtoDT(npy_int64 year, char relation, asfreq_info *af_info) {
    npy_int64 absdate;
    int month = (af_info->from_a_year_end) % 12;

    // start from 1970
    year += BASE_YEAR;

    month += 1;

    if (af_info->from_a_year_end != 12) {
      year -= 1;
    }

    if (relation == 'E') {
      year += 1;
    }

    absdate = absdate_from_ymd(year, month, 1);

    if (absdate  == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    if (relation == 'E') {
      absdate -= 1;
    }

    return upsample_daytime(absdate - ORD_OFFSET, af_info, relation != 'S');
}

static npy_int64 asfreq_AtoA(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_AtoQ(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_AtoM(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_AtoW(npy_int64 ordinal, char relation, asfreq_info *af_info) {
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_AtoB(npy_int64 ordinal, char relation, asfreq_info *af_info) {

    struct date_info dinfo;
    if (dInfoCalc_SetFromAbsDate(&dinfo,
                asfreq_AtoDT(ordinal, relation, af_info) + ORD_OFFSET,
                GREGORIAN_CALENDAR)) return INT_ERR_CODE;

    if (relation == 'S') { return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week); }
    else                 { return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week); }
}

static npy_int64 nofunc(npy_int64 ordinal, char relation, asfreq_info *af_info) { return INT_ERR_CODE; }
static npy_int64 no_op(npy_int64 ordinal, char relation, asfreq_info *af_info) { return ordinal; }

// end of frequency specific conversion routines

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

	af_info->intraday_conversion_factor =
	    get_daytime_conversion_factor(
	        get_freq_group_index(max_value(fromGroup, FR_DAY)),
	        get_freq_group_index(max_value(toGroup, FR_DAY))
	    );

    //printf("get_asfreq_info(%d, %d) %ld, %d\n", fromFreq, toFreq, af_info->intraday_conversion_factor, af_info->intraday_conversion_upsample);

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
                case FR_DAY: 
                case FR_HR: 
                case FR_MIN: 
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                             return &asfreq_AtoDT;

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
                case FR_DAY: 
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                             return &asfreq_QtoDT;
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
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                             return &asfreq_MtoDT;
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
                case FR_DAY: 
                case FR_HR: 
                case FR_MIN: 
                case FR_SEC: 
                case FR_MS:
                case FR_US:
                case FR_NS:
                             return &asfreq_WtoDT;
                default: return &nofunc;
            }

        case FR_BUS:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_BtoA;
                case FR_QTR: return &asfreq_BtoQ;
                case FR_MTH: return &asfreq_BtoM;
                case FR_WK: return &asfreq_BtoW;
                case FR_BUS: return &no_op;
                case FR_DAY: 
                case FR_HR: 
                case FR_MIN: 
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                             return &asfreq_BtoDT;
                default: return &nofunc;
            }

        case FR_DAY:
        case FR_HR:
        case FR_MIN:
        case FR_SEC:
        case FR_MS:
        case FR_US:
        case FR_NS:
            switch(toGroup)
            {
                case FR_ANN: return &asfreq_DTtoA;
                case FR_QTR: return &asfreq_DTtoQ;
                case FR_MTH: return &asfreq_DTtoM;
                case FR_WK: return &asfreq_DTtoW;
                case FR_BUS: return &asfreq_DTtoB;
                case FR_DAY: 
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    if (fromGroup > toGroup) {
                        return &asfreq_DownsampleWithinDay;
                    } else {
                        return &asfreq_UpsampleWithinDay;
                    }
                default: return &nofunc;
            }

        default: return &nofunc;
    }
}

double get_abs_time(int freq, npy_int64 date_ordinal, npy_int64 ordinal) {
    //printf("get_abs_time %d %lld %lld\n", freq, date_ordinal, ordinal);

	int freq_index, day_index, base_index;
	npy_int64 per_day, start_ord;
	double unit, result;

    if (freq <= FR_DAY) {
      return 0;
    }

    freq_index = get_freq_group_index(freq);
    day_index = get_freq_group_index(FR_DAY);
    base_index = get_freq_group_index(FR_SEC);

    //printf("  indices: day %d, freq %d, base %d\n", day_index, freq_index, base_index);

    per_day = get_daytime_conversion_factor(day_index, freq_index);
    unit = get_daytime_conversion_factor(freq_index, base_index);

    //printf("  per_day: %lld, unit: %f\n", per_day, unit);

    if (base_index < freq_index) {
      unit = 1 / unit;
      //printf("  corrected unit: %f\n", unit);
    }

    start_ord = date_ordinal * per_day;
    //printf("start_ord: %lld\n", start_ord);
    result = (double) ( unit * (ordinal - start_ord));
    //printf("  result: %f\n", result);
    return result;
}

/* Sets the time part of the DateTime object. */
static int dInfoCalc_SetFromAbsTime(struct date_info *dinfo,
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
static int dInfoCalc_SetFromAbsDateTime(struct date_info *dinfo,
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

    //printf("\n%x %d %d %ld %ld\n", func, freq1, freq2, finfo.intraday_conversion_factor, -finfo.intraday_conversion_factor);

    val = (*func)(period_ordinal, relation, &finfo);

    if (val == INT_ERR_CODE) {
        //Py_Error(PyExc_ValueError, "Unable to convert to desired frequency.");
        goto onError;
    }
    return val;
onError:
    return INT_ERR_CODE;
}


/* generate an ordinal in period space */
npy_int64 get_period_ordinal(int year, int month, int day,
        int hour, int minute, int second, int microseconds, int picoseconds,
        int freq)
{
    npy_int64 absdays, delta, seconds;
    npy_int64 weeks, days;
    npy_int64 ordinal, day_adj;
    int freq_group, fmonth, mdiff;
    freq_group = get_freq_group(freq);

    if (freq == FR_SEC || freq == FR_MS || freq == FR_US || freq == FR_NS) {

        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        seconds = (npy_int64)(delta * 86400 + hour * 3600 + minute * 60 + second);

        switch(freq) {
          case FR_MS:
            return seconds * 1000 + microseconds / 1000;

          case FR_US:
            return seconds * 1000000 + microseconds;

          case FR_NS:
            return seconds * 1000000000 + microseconds * 1000 + picoseconds / 1000;
        }

        return seconds;
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
        // calculate the current week assuming sunday as last day of a week
        weeks = (days - BASE_WEEK_TO_DAY_OFFSET) / DAYS_PER_WEEK;
        // calculate the current weekday (in range 1 .. 7)
        delta = (days - BASE_WEEK_TO_DAY_OFFSET) % DAYS_PER_WEEK + 1;
        // return the number of business days in full weeks plus the business days in the last - possible partial - week
        return (npy_int64)(weeks * BUSINESS_DAYS_PER_WEEK)
            + (delta <= BUSINESS_DAYS_PER_WEEK
                ? delta
                : BUSINESS_DAYS_PER_WEEK + 1)
             - BDAY_OFFSET;
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
	freq_conv_func toDaily = NULL;

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
    double abstime = get_abs_time(freq, absdate - ORD_OFFSET, ordinal);

    while (abstime < 0) {
        abstime += 86400;
        absdate -= 1;
    }
    while (abstime >= 86400) {
        abstime -= 86400;
	absdate += 1;
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

int pdays_in_month(npy_int64 ordinal, int freq) {
    int days;
    struct date_info dinfo;
    if(get_date_info(ordinal, freq, &dinfo) == INT_ERR_CODE)
        return INT_ERR_CODE;
    
    days = days_in_month[dInfoCalc_Leapyear(dinfo.year, dinfo.calendar)][dinfo.month-1];
    return days;
}
