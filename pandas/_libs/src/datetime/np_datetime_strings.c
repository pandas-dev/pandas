/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Written by Mark Wiebe (mwwiebe@gmail.com)
Copyright (c) 2011 by Enthought, Inc.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

See NUMPY_LICENSE.txt for the license.

This file implements string parsing and creation for NumPy datetime.

*/

#define PY_SSIZE_T_CLEAN
#define NO_IMPORT

#include <Python.h>

#include <time.h>

#include <numpy/arrayobject.h>
#include "numpy/arrayscalars.h"

#include "np_datetime.h"
#include "np_datetime_strings.h"


/* Platform-specific time_t typedef */
typedef time_t NPY_TIME_T;

/* We *do* want these symbols, but for Cython, not for C.
   Fine in Mac OSX, but Linux complains.

static void _suppress_unused_variable_warning(void) {
    int x = days_per_month_table[0][0];
    x = x;

    int y = _month_offset[0][0];
    y = y;

    char *z = _datetime_strings[0];
    z = z;
} */

/* Exported as DATETIMEUNITS in multiarraymodule.c */
static char *_datetime_strings[PANDAS_DATETIME_NUMUNITS] = {
    "Y", "M", "W", "D", "h", "m", "s", "ms", "us", "ns", "ps", "fs", "as",
};
/*
 * Wraps `localtime` functionality for multiple platforms. This
 * converts a time value to a time structure in the local timezone.
 *
 * Returns 0 on success, -1 on failure.
 */
static int get_localtime(NPY_TIME_T *ts, struct tm *tms) {
    char *func_name = "<unknown>";
#if defined(_WIN32)
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
    if (localtime_s(tms, ts) != 0) {
        func_name = "localtime_s";
        goto fail;
    }
#elif defined(__GNUC__) && defined(NPY_MINGW_USE_CUSTOM_MSVCR)
    if (_localtime64_s(tms, ts) != 0) {
        func_name = "_localtime64_s";
        goto fail;
    }
#else
    struct tm *tms_tmp;
    localtime_r(ts, tms_tmp);
    if (tms_tmp == NULL) {
        func_name = "localtime";
        goto fail;
    }
    memcpy(tms, tms_tmp, sizeof(struct tm));
#endif
#else
    if (localtime_r(ts, tms) == NULL) {
        func_name = "localtime_r";
        goto fail;
    }
#endif

    return 0;

fail:
    PyErr_Format(PyExc_OSError,
                 "Failed to use '%s' to convert "
                 "to a local time",
                 func_name);
    return -1;
}


/*
 * Converts a datetimestruct in UTC to a datetimestruct in local time,
 * also returning the timezone offset applied.
 *
 * Returns 0 on success, -1 on failure.
 */
static int convert_datetimestruct_utc_to_local(
    pandas_datetimestruct *out_dts_local, const pandas_datetimestruct *dts_utc,
    int *out_timezone_offset) {
    NPY_TIME_T rawtime = 0, localrawtime;
    struct tm tm_;
    npy_int64 year_correction = 0;

    /* Make a copy of the input 'dts' to modify */
    *out_dts_local = *dts_utc;

    /* HACK: Use a year < 2038 for later years for small time_t */
    if (sizeof(NPY_TIME_T) == 4 && out_dts_local->year >= 2038) {
        if (is_leapyear(out_dts_local->year)) {
            /* 2036 is a leap year */
            year_correction = out_dts_local->year - 2036;
            out_dts_local->year -= year_correction;
        } else {
            /* 2037 is not a leap year */
            year_correction = out_dts_local->year - 2037;
            out_dts_local->year -= year_correction;
        }
    }

    /*
     * Convert everything in 'dts' to a time_t, to minutes precision.
     * This is POSIX time, which skips leap-seconds, but because
     * we drop the seconds value from the pandas_datetimestruct, everything
     * is ok for this operation.
     */
    rawtime = (time_t)get_datetimestruct_days(out_dts_local) * 24 * 60 * 60;
    rawtime += dts_utc->hour * 60 * 60;
    rawtime += dts_utc->min * 60;

    /* localtime converts a 'time_t' into a local 'struct tm' */
    if (get_localtime(&rawtime, &tm_) < 0) {
        return -1;
    }

    /* Copy back all the values except seconds */
    out_dts_local->min = tm_.tm_min;
    out_dts_local->hour = tm_.tm_hour;
    out_dts_local->day = tm_.tm_mday;
    out_dts_local->month = tm_.tm_mon + 1;
    out_dts_local->year = tm_.tm_year + 1900;

    /* Extract the timezone offset that was applied */
    rawtime /= 60;
    localrawtime = (time_t)get_datetimestruct_days(out_dts_local) * 24 * 60;
    localrawtime += out_dts_local->hour * 60;
    localrawtime += out_dts_local->min;

    *out_timezone_offset = localrawtime - rawtime;

    /* Reapply the year 2038 year correction HACK */
    out_dts_local->year += year_correction;

    return 0;
}


/*
 * Parses (almost) standard ISO 8601 date strings. The differences are:
 *
 * + Only seconds may have a decimal point, with up to 18 digits after it
 *   (maximum attoseconds precision).
 * + Either a 'T' as in ISO 8601 or a ' ' may be used to separate
 *   the date and the time. Both are treated equivalently.
 * + Doesn't (yet) handle the "YYYY-DDD" or "YYYY-Www" formats.
 * + Doesn't handle leap seconds (seconds value has 60 in these cases).
 * + Doesn't handle 24:00:00 as synonym for midnight (00:00:00) tomorrow
 * + Accepts special values "NaT" (not a time), "Today", (current
 *   day according to local time) and "Now" (current time in UTC).
 * + ':' separator between hours, minutes, and seconds is optional. When
 *   omitted, each component must be 2 digits if it appears. (GH-10041)
 *
 * 'str' must be a NULL-terminated string, and 'len' must be its length.
 * 'unit' should contain -1 if the unit is unknown, or the unit
 *      which will be used if it is.
 *
 * 'out' gets filled with the parsed date-time.
 * 'out_local' gets set to 1 if the parsed time contains timezone,
 *      to 0 otherwise.
 * 'out_tzoffset' gets set to timezone offset by minutes
 *      if the parsed time was in local time,
 *      to 0 otherwise. The values 'now' and 'today' don't get counted
 *      as local, and neither do UTC +/-#### timezone offsets, because
 *      they aren't using the computer's local timezone offset.
 * 'out_bestunit' gives a suggested unit based on the amount of
 *      resolution provided in the string, or -1 for NaT.
 * 'out_special' gets set to 1 if the parsed time was 'today',
 *      'now', or ''/'NaT'. For 'today', the unit recommended is
 *      'D', for 'now', the unit recommended is 's', and for 'NaT'
 *      the unit recommended is 'Y'.
 *
 * Returns 0 on success, -1 on failure.
 */
int parse_iso_8601_datetime(char *str, int len, PANDAS_DATETIMEUNIT unit,
                            pandas_datetimestruct *out,
                            int *out_local, int *out_tzoffset,
                            PANDAS_DATETIMEUNIT *out_bestunit,
                            npy_bool *out_special) {
    int year_leap = 0;
    int i, numdigits;
    char *substr, sublen;
    PANDAS_DATETIMEUNIT bestunit;

    /* If year-month-day are separated by a valid separator,
     * months/days without leading zeroes will be parsed
     * (though not iso8601). If the components aren't separated,
     * 4 (YYYY) or 8 (YYYYMMDD) digits are expected. 6 digits are
     * forbidden here (but parsed as YYMMDD elsewhere).
    */
    int has_ymd_sep = 0;
    char ymd_sep = '\0';
    char valid_ymd_sep[] = {'-', '.', '/', '\\', ' '};
    int valid_ymd_sep_len = sizeof(valid_ymd_sep);

    /* hour-minute-second may or may not separated by ':'. If not, then
     * each component must be 2 digits. */
    int has_hms_sep = 0;
    int hour_was_2_digits = 0;

    /* Initialize the output to all zeros */
    memset(out, 0, sizeof(pandas_datetimestruct));
    out->month = 1;
    out->day = 1;

    /*
     * The string "today" means take today's date in local time, and
     * convert it to a date representation. This date representation, if
     * forced into a time unit, will be at midnight UTC.
     * This is perhaps a little weird, but done so that the
     * 'datetime64[D]' type produces the date you expect, rather than
     * switching to an adjacent day depending on the current time and your
     * timezone.
     */
    if (len == 5 && tolower(str[0]) == 't' && tolower(str[1]) == 'o' &&
        tolower(str[2]) == 'd' && tolower(str[3]) == 'a' &&
        tolower(str[4]) == 'y') {
        NPY_TIME_T rawtime = 0;
        struct tm tm_;

        time(&rawtime);
        if (get_localtime(&rawtime, &tm_) < 0) {
            return -1;
        }
        out->year = tm_.tm_year + 1900;
        out->month = tm_.tm_mon + 1;
        out->day = tm_.tm_mday;

        bestunit = PANDAS_FR_D;

        /*
         * Indicate that this was a special value, and
         * is a date (unit 'D').
         */
        if (out_local != NULL) {
            *out_local = 0;
        }
        if (out_bestunit != NULL) {
            *out_bestunit = bestunit;
        }
        if (out_special != NULL) {
            *out_special = 1;
        }

        return 0;
    }

    /* The string "now" resolves to the current UTC time */
    if (len == 3 && tolower(str[0]) == 'n' && tolower(str[1]) == 'o' &&
        tolower(str[2]) == 'w') {
        NPY_TIME_T rawtime = 0;
        pandas_datetime_metadata meta;

        time(&rawtime);

        /* Set up a dummy metadata for the conversion */
        meta.base = PANDAS_FR_s;
        meta.num = 1;

        bestunit = PANDAS_FR_s;

        /*
         * Indicate that this was a special value, and
         * use 's' because the time() function has resolution
         * seconds.
         */
        if (out_local != NULL) {
            *out_local = 0;
        }
        if (out_bestunit != NULL) {
            *out_bestunit = bestunit;
        }
        if (out_special != NULL) {
            *out_special = 1;
        }

        return convert_datetime_to_datetimestruct(&meta, rawtime, out);
    }

    /* Anything else isn't a special value */
    if (out_special != NULL) {
        *out_special = 0;
    }

    substr = str;
    sublen = len;

    /* Skip leading whitespace */
    while (sublen > 0 && isspace(*substr)) {
        ++substr;
        --sublen;
    }

    /* Leading '-' sign for negative year */
    if (*substr == '-') {
        ++substr;
        --sublen;
    }

    if (sublen == 0) {
        goto parse_error;
    }

    /* PARSE THE YEAR (4 digits) */
    out->year = 0;
    if (sublen >= 4 && isdigit(substr[0]) && isdigit(substr[1]) &&
        isdigit(substr[2]) && isdigit(substr[3])) {
        out->year = 1000 * (substr[0] - '0') + 100 * (substr[1] - '0') +
                    10 * (substr[2] - '0') + (substr[3] - '0');

        substr += 4;
        sublen -= 4;
    }

    /* Negate the year if necessary */
    if (str[0] == '-') {
        out->year = -out->year;
    }
    /* Check whether it's a leap-year */
    year_leap = is_leapyear(out->year);

    /* Next character must be a separator, start of month, or end of string */
    if (sublen == 0) {
        if (out_local != NULL) {
            *out_local = 0;
        }
        bestunit = PANDAS_FR_Y;
        goto finish;
    }

    if (!isdigit(*substr)) {
        for (i = 0; i < valid_ymd_sep_len; ++i) {
            if (*substr == valid_ymd_sep[i]) {
                break;
            }
        }
        if (i == valid_ymd_sep_len) {
            goto parse_error;
        }
        has_ymd_sep = 1;
        ymd_sep = valid_ymd_sep[i];
        ++substr;
        --sublen;
        /* Cannot have trailing separator */
        if (sublen == 0 || !isdigit(*substr)) {
            goto parse_error;
        }
    }

    /* PARSE THE MONTH */
    /* First digit required */
    out->month = (*substr - '0');
    ++substr;
    --sublen;
    /* Second digit optional if there was a separator */
    if (isdigit(*substr)) {
        out->month = 10 * out->month + (*substr - '0');
        ++substr;
        --sublen;
    } else if (!has_ymd_sep) {
        goto parse_error;
    }
    if (out->month < 1 || out->month > 12) {
        PyErr_Format(PyExc_ValueError,
                     "Month out of range in datetime string \"%s\"", str);
        goto error;
    }

    /* Next character must be the separator, start of day, or end of string */
    if (sublen == 0) {
        /* Forbid YYYYMM. Parsed instead as YYMMDD by someone else. */
        if (!has_ymd_sep) {
            goto parse_error;
        }
        if (out_local != NULL) {
            *out_local = 0;
        }
        bestunit = PANDAS_FR_M;
        goto finish;
    }

    if (has_ymd_sep) {
        /* Must have separator, but cannot be trailing */
        if (*substr != ymd_sep || sublen == 1) {
            goto parse_error;
        }
        ++substr;
        --sublen;
    }

    /* PARSE THE DAY */
    /* First digit required */
    if (!isdigit(*substr)) {
        goto parse_error;
    }
    out->day = (*substr - '0');
    ++substr;
    --sublen;
    /* Second digit optional if there was a separator */
    if (isdigit(*substr)) {
        out->day = 10 * out->day + (*substr - '0');
        ++substr;
        --sublen;
    } else if (!has_ymd_sep) {
        goto parse_error;
    }
    if (out->day < 1 ||
        out->day > days_per_month_table[year_leap][out->month - 1]) {
        PyErr_Format(PyExc_ValueError,
                     "Day out of range in datetime string \"%s\"", str);
        goto error;
    }

    /* Next character must be a 'T', ' ', or end of string */
    if (sublen == 0) {
        if (out_local != NULL) {
            *out_local = 0;
        }
        bestunit = PANDAS_FR_D;
        goto finish;
    }

    if ((*substr != 'T' && *substr != ' ') || sublen == 1) {
        goto parse_error;
    }
    ++substr;
    --sublen;

    /* PARSE THE HOURS */
    /* First digit required */
    if (!isdigit(*substr)) {
        goto parse_error;
    }
    out->hour = (*substr - '0');
    ++substr;
    --sublen;
    /* Second digit optional */
    if (isdigit(*substr)) {
        hour_was_2_digits = 1;
        out->hour = 10 * out->hour + (*substr - '0');
        ++substr;
        --sublen;
        if (out->hour >= 24) {
            PyErr_Format(PyExc_ValueError,
                         "Hours out of range in datetime string \"%s\"", str);
            goto error;
        }
    }

    /* Next character must be a ':' or the end of the string */
    if (sublen == 0) {
        if (!hour_was_2_digits) {
            goto parse_error;
        }
        bestunit = PANDAS_FR_h;
        goto finish;
    }

    if (*substr == ':') {
        has_hms_sep = 1;
        ++substr;
        --sublen;
        /* Cannot have a trailing separator */
        if (sublen == 0 || !isdigit(*substr)) {
            goto parse_error;
        }
    } else if (!isdigit(*substr)) {
        if (!hour_was_2_digits) {
            goto parse_error;
        }
        bestunit = PANDAS_FR_h;
        goto parse_timezone;
    }

    /* PARSE THE MINUTES */
    /* First digit required */
    out->min = (*substr - '0');
    ++substr;
    --sublen;
    /* Second digit optional if there was a separator */
    if (isdigit(*substr)) {
        out->min = 10 * out->min + (*substr - '0');
        ++substr;
        --sublen;
        if (out->min >= 60) {
            PyErr_Format(PyExc_ValueError,
                         "Minutes out of range in datetime string \"%s\"", str);
            goto error;
        }
    } else if (!has_hms_sep) {
        goto parse_error;
    }

    if (sublen == 0) {
        bestunit = PANDAS_FR_m;
        goto finish;
    }

    /* If we make it through this condition block, then the next
     * character is a digit. */
    if (has_hms_sep && *substr == ':') {
        ++substr;
        --sublen;
        /* Cannot have a trailing ':' */
        if (sublen == 0 || !isdigit(*substr)) {
            goto parse_error;
        }
    } else if (!has_hms_sep && isdigit(*substr)) {
    } else {
        bestunit = PANDAS_FR_m;
        goto parse_timezone;
    }

    /* PARSE THE SECONDS */
    /* First digit required */
    out->sec = (*substr - '0');
    ++substr;
    --sublen;
    /* Second digit optional if there was a separator */
    if (isdigit(*substr)) {
        out->sec = 10 * out->sec + (*substr - '0');
        ++substr;
        --sublen;
        if (out->sec >= 60) {
            PyErr_Format(PyExc_ValueError,
                         "Seconds out of range in datetime string \"%s\"", str);
            goto error;
        }
    } else if (!has_hms_sep) {
        goto parse_error;
    }

    /* Next character may be a '.' indicating fractional seconds */
    if (sublen > 0 && *substr == '.') {
        ++substr;
        --sublen;
    } else {
        bestunit = PANDAS_FR_s;
        goto parse_timezone;
    }

    /* PARSE THE MICROSECONDS (0 to 6 digits) */
    numdigits = 0;
    for (i = 0; i < 6; ++i) {
        out->us *= 10;
        if (sublen > 0 && isdigit(*substr)) {
            out->us += (*substr - '0');
            ++substr;
            --sublen;
            ++numdigits;
        }
    }

    if (sublen == 0 || !isdigit(*substr)) {
        if (numdigits > 3) {
            bestunit = PANDAS_FR_us;
        } else {
            bestunit = PANDAS_FR_ms;
        }
        goto parse_timezone;
    }

    /* PARSE THE PICOSECONDS (0 to 6 digits) */
    numdigits = 0;
    for (i = 0; i < 6; ++i) {
        out->ps *= 10;
        if (sublen > 0 && isdigit(*substr)) {
            out->ps += (*substr - '0');
            ++substr;
            --sublen;
            ++numdigits;
        }
    }

    if (sublen == 0 || !isdigit(*substr)) {
        if (numdigits > 3) {
            bestunit = PANDAS_FR_ps;
        } else {
            bestunit = PANDAS_FR_ns;
        }
        goto parse_timezone;
    }

    /* PARSE THE ATTOSECONDS (0 to 6 digits) */
    numdigits = 0;
    for (i = 0; i < 6; ++i) {
        out->as *= 10;
        if (sublen > 0 && isdigit(*substr)) {
            out->as += (*substr - '0');
            ++substr;
            --sublen;
            ++numdigits;
        }
    }

    if (numdigits > 3) {
        bestunit = PANDAS_FR_as;
    } else {
        bestunit = PANDAS_FR_fs;
    }

parse_timezone:
    /* trim any whitepsace between time/timeezone */
    while (sublen > 0 && isspace(*substr)) {
        ++substr;
        --sublen;
    }

    if (sublen == 0) {
        // Unlike NumPy, treating no time zone as naive
        goto finish;
    }

    /* UTC specifier */
    if (*substr == 'Z') {
        /* "Z" should be equivalent to tz offset "+00:00" */
        if (out_local != NULL) {
            *out_local = 1;
        }

        if (out_tzoffset != NULL) {
            *out_tzoffset = 0;
        }

        if (sublen == 1) {
            goto finish;
        } else {
            ++substr;
            --sublen;
        }
    } else if (*substr == '-' || *substr == '+') {
        /* Time zone offset */
        int offset_neg = 0, offset_hour = 0, offset_minute = 0;

        /*
         * Since "local" means local with respect to the current
         * machine, we say this is non-local.
         */

        if (*substr == '-') {
            offset_neg = 1;
        }
        ++substr;
        --sublen;

        /* The hours offset */
        if (sublen >= 2 && isdigit(substr[0]) && isdigit(substr[1])) {
            offset_hour = 10 * (substr[0] - '0') + (substr[1] - '0');
            substr += 2;
            sublen -= 2;
            if (offset_hour >= 24) {
                PyErr_Format(PyExc_ValueError,
                             "Timezone hours offset out of range "
                             "in datetime string \"%s\"",
                             str);
                goto error;
            }
        } else if (sublen >= 1 && isdigit(substr[0])) {
            offset_hour = substr[0] - '0';
            ++substr;
            --sublen;
        } else {
            goto parse_error;
        }

        /* The minutes offset is optional */
        if (sublen > 0) {
            /* Optional ':' */
            if (*substr == ':') {
                ++substr;
                --sublen;
            }

            /* The minutes offset (at the end of the string) */
            if (sublen >= 2 && isdigit(substr[0]) && isdigit(substr[1])) {
                offset_minute = 10 * (substr[0] - '0') + (substr[1] - '0');
                substr += 2;
                sublen -= 2;
                if (offset_minute >= 60) {
                    PyErr_Format(PyExc_ValueError,
                                 "Timezone minutes offset out of range "
                                 "in datetime string \"%s\"",
                                 str);
                    goto error;
                }
            } else if (sublen >= 1 && isdigit(substr[0])) {
                offset_minute = substr[0] - '0';
                ++substr;
                --sublen;
            } else {
                goto parse_error;
            }
        }

        /* Apply the time zone offset */
        if (offset_neg) {
            offset_hour = -offset_hour;
            offset_minute = -offset_minute;
        }
        if (out_local != NULL) {
            *out_local = 1;
            // Unlike NumPy, do not change internal value to local time
            *out_tzoffset = 60 * offset_hour + offset_minute;
        }
    }

    /* Skip trailing whitespace */
    while (sublen > 0 && isspace(*substr)) {
        ++substr;
        --sublen;
    }

    if (sublen != 0) {
        goto parse_error;
    }

finish:
    if (out_bestunit != NULL) {
        *out_bestunit = bestunit;
    }

    return 0;

parse_error:
    PyErr_Format(PyExc_ValueError,
                 "Error parsing datetime string \"%s\" at position %d", str,
                 (int)(substr - str));
    return -1;

error:
    return -1;
}

/*
 * Provides a string length to use for converting datetime
 * objects with the given local and unit settings.
 */
int get_datetime_iso_8601_strlen(int local, PANDAS_DATETIMEUNIT base) {
    int len = 0;

    switch (base) {
        /* Generic units can only be used to represent NaT */
        /*case PANDAS_FR_GENERIC:*/
        /*    return 4;*/
        case PANDAS_FR_as:
            len += 3; /* "###" */
        case PANDAS_FR_fs:
            len += 3; /* "###" */
        case PANDAS_FR_ps:
            len += 3; /* "###" */
        case PANDAS_FR_ns:
            len += 3; /* "###" */
        case PANDAS_FR_us:
            len += 3; /* "###" */
        case PANDAS_FR_ms:
            len += 4; /* ".###" */
        case PANDAS_FR_s:
            len += 3; /* ":##" */
        case PANDAS_FR_m:
            len += 3; /* ":##" */
        case PANDAS_FR_h:
            len += 3; /* "T##" */
        case PANDAS_FR_D:
        case PANDAS_FR_W:
            len += 3; /* "-##" */
        case PANDAS_FR_M:
            len += 3; /* "-##" */
        case PANDAS_FR_Y:
            len += 21; /* 64-bit year */
            break;
        default:
            len += 3; /* handle the now defunct NPY_FR_B */
            break;
    }

    if (base >= PANDAS_FR_h) {
        if (local) {
            len += 5; /* "+####" or "-####" */
        } else {
            len += 1; /* "Z" */
        }
    }

    len += 1; /* NULL terminator */

    return len;
}


/*
 * Converts an pandas_datetimestruct to an (almost) ISO 8601
 * NULL-terminated string. If the string fits in the space exactly,
 * it leaves out the NULL terminator and returns success.
 *
 * The differences from ISO 8601 are the 'NaT' string, and
 * the number of year digits is >= 4 instead of strictly 4.
 *
 * If 'local' is non-zero, it produces a string in local time with
 * a +-#### timezone offset, otherwise it uses timezone Z (UTC).
 *
 * 'base' restricts the output to that unit. Set 'base' to
 * -1 to auto-detect a base after which all the values are zero.
 *
 *  'tzoffset' is used if 'local' is enabled, and 'tzoffset' is
 *  set to a value other than -1. This is a manual override for
 *  the local time zone to use, as an offset in minutes.
 *
 *  Returns 0 on success, -1 on failure (for example if the output
 *  string was too short).
 */
int make_iso_8601_datetime(pandas_datetimestruct *dts, char *outstr, int outlen,
                           int local, PANDAS_DATETIMEUNIT base, int tzoffset) {
    pandas_datetimestruct dts_local;
    int timezone_offset = 0;

    char *substr = outstr, sublen = outlen;
    int tmplen;

    /* Only do local time within a reasonable year range */
    if ((dts->year <= 1800 || dts->year >= 10000) && tzoffset == -1) {
        local = 0;
    }

    /*
     * Print weeks with the same precision as days.
     *
     * TODO: Could print weeks with YYYY-Www format if the week
     *       epoch is a Monday.
     */
    if (base == PANDAS_FR_W) {
        base = PANDAS_FR_D;
    }

    /* Use the C API to convert from UTC to local time */
    if (local && tzoffset == -1) {
        if (convert_datetimestruct_utc_to_local(&dts_local, dts,
                                                &timezone_offset) < 0) {
            return -1;
        }

        /* Set dts to point to our local time instead of the UTC time */
        dts = &dts_local;
    } else if (local) {
        // Use the manually provided tzoffset.
        // Make a copy of the pandas_datetimestruct we can modify.
        dts_local = *dts;
        dts = &dts_local;

        /* Set and apply the required timezone offset */
        timezone_offset = tzoffset;
        add_minutes_to_datetimestruct(dts, timezone_offset);
    }

/* YEAR */
/*
 * Can't use PyOS_snprintf, because it always produces a '\0'
 * character at the end, and NumPy string types are permitted
 * to have data all the way to the end of the buffer.
 */
#ifdef _WIN32
    tmplen = _snprintf(substr, sublen, "%04" NPY_INT64_FMT, dts->year);
#else
    tmplen = snprintf(substr, sublen, "%04" NPY_INT64_FMT, dts->year);
#endif
    /* If it ran out of space or there isn't space for the NULL terminator */
    if (tmplen < 0 || tmplen > sublen) {
        goto string_too_short;
    }
    substr += tmplen;
    sublen -= tmplen;

    /* Stop if the unit is years */
    if (base == PANDAS_FR_Y) {
        if (sublen > 0) {
            *substr = '\0';
        }
        return 0;
    }

    /* MONTH */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = '-';
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->month / 10) + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->month % 10) + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is months */
    if (base == PANDAS_FR_M) {
        if (sublen > 0) {
            *substr = '\0';
        }
        return 0;
    }

    /* DAY */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = '-';
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->day / 10) + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->day % 10) + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is days */
    if (base == PANDAS_FR_D) {
        if (sublen > 0) {
            *substr = '\0';
        }
        return 0;
    }

    /* HOUR */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = 'T';
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->hour / 10) + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->hour % 10) + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is hours */
    if (base == PANDAS_FR_h) {
        goto add_time_zone;
    }

    /* MINUTE */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = ':';
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->min / 10) + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->min % 10) + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is minutes */
    if (base == PANDAS_FR_m) {
        goto add_time_zone;
    }

    /* SECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = ':';
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->sec / 10) + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->sec % 10) + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is seconds */
    if (base == PANDAS_FR_s) {
        goto add_time_zone;
    }

    /* MILLISECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = '.';
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->us / 100000) % 10 + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->us / 10000) % 10 + '0');
    if (sublen < 4) {
        goto string_too_short;
    }
    substr[3] = (char)((dts->us / 1000) % 10 + '0');
    substr += 4;
    sublen -= 4;

    /* Stop if the unit is milliseconds */
    if (base == PANDAS_FR_ms) {
        goto add_time_zone;
    }

    /* MICROSECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = (char)((dts->us / 100) % 10 + '0');
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->us / 10) % 10 + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)(dts->us % 10 + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is microseconds */
    if (base == PANDAS_FR_us) {
        goto add_time_zone;
    }

    /* NANOSECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = (char)((dts->ps / 100000) % 10 + '0');
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->ps / 10000) % 10 + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->ps / 1000) % 10 + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is nanoseconds */
    if (base == PANDAS_FR_ns) {
        goto add_time_zone;
    }

    /* PICOSECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = (char)((dts->ps / 100) % 10 + '0');
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->ps / 10) % 10 + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)(dts->ps % 10 + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is picoseconds */
    if (base == PANDAS_FR_ps) {
        goto add_time_zone;
    }

    /* FEMTOSECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = (char)((dts->as / 100000) % 10 + '0');
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->as / 10000) % 10 + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)((dts->as / 1000) % 10 + '0');
    substr += 3;
    sublen -= 3;

    /* Stop if the unit is femtoseconds */
    if (base == PANDAS_FR_fs) {
        goto add_time_zone;
    }

    /* ATTOSECOND */
    if (sublen < 1) {
        goto string_too_short;
    }
    substr[0] = (char)((dts->as / 100) % 10 + '0');
    if (sublen < 2) {
        goto string_too_short;
    }
    substr[1] = (char)((dts->as / 10) % 10 + '0');
    if (sublen < 3) {
        goto string_too_short;
    }
    substr[2] = (char)(dts->as % 10 + '0');
    substr += 3;
    sublen -= 3;

add_time_zone:
    if (local) {
        /* Add the +/- sign */
        if (sublen < 1) {
            goto string_too_short;
        }
        if (timezone_offset < 0) {
            substr[0] = '-';
            timezone_offset = -timezone_offset;
        } else {
            substr[0] = '+';
        }
        substr += 1;
        sublen -= 1;

        /* Add the timezone offset */
        if (sublen < 1) {
            goto string_too_short;
        }
        substr[0] = (char)((timezone_offset / (10 * 60)) % 10 + '0');
        if (sublen < 2) {
            goto string_too_short;
        }
        substr[1] = (char)((timezone_offset / 60) % 10 + '0');
        if (sublen < 3) {
            goto string_too_short;
        }
        substr[2] = (char)(((timezone_offset % 60) / 10) % 10 + '0');
        if (sublen < 4) {
            goto string_too_short;
        }
        substr[3] = (char)((timezone_offset % 60) % 10 + '0');
        substr += 4;
        sublen -= 4;
    } else {
        /* UTC "Zulu" time */
        if (sublen < 1) {
            goto string_too_short;
        }
        substr[0] = 'Z';
        substr += 1;
        sublen -= 1;
    }

    /* Add a NULL terminator, and return */
    if (sublen > 0) {
        substr[0] = '\0';
    }

    return 0;

string_too_short:
    PyErr_Format(PyExc_RuntimeError,
                 "The string provided for NumPy ISO datetime formatting "
                 "was too short, with length %d",
                 outlen);
    return -1;
}
