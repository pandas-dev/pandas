/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

This file is derived from NumPy 1.7. See NUMPY_LICENSE.txt

*/

// Licence at LICENSES/NUMPY_LICENSE

#define NO_IMPORT

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif // NPY_NO_DEPRECATED_API

#include <Python.h>

#include "pandas/vendored/numpy/datetime/np_datetime.h"
#include <numpy/ndarraytypes.h>
#include <numpy/npy_common.h>

#if defined(_WIN32)
#ifndef ENABLE_INTSAFE_SIGNED_FUNCTIONS
#define ENABLE_INTSAFE_SIGNED_FUNCTIONS
#endif
#include <intsafe.h>
#define checked_int64_add(a, b, res) LongLongAdd(a, b, res)
#define checked_int64_sub(a, b, res) LongLongSub(a, b, res)
#define checked_int64_mul(a, b, res) LongLongMult(a, b, res)
#else
#if defined __has_builtin
#if __has_builtin(__builtin_add_overflow)
#define checked_int64_add(a, b, res) __builtin_add_overflow(a, b, res)
#define checked_int64_sub(a, b, res) __builtin_sub_overflow(a, b, res)
#define checked_int64_mul(a, b, res) __builtin_mul_overflow(a, b, res)
#else
_Static_assert(0,
               "Overflow checking not detected; please try a newer compiler");
#endif
// __has_builtin was added in gcc 10, but our muslinux_1_1 build environment
// only has gcc-9.3, so fall back to __GNUC__ macro as long as we have that
#elif __GNUC__ > 7
#define checked_int64_add(a, b, res) __builtin_add_overflow(a, b, res)
#define checked_int64_sub(a, b, res) __builtin_sub_overflow(a, b, res)
#define checked_int64_mul(a, b, res) __builtin_mul_overflow(a, b, res)
#else
_Static_assert(0, "__has_builtin not detected; please try a newer compiler");
#endif
#endif

#define PD_CHECK_OVERFLOW(FUNC)                                                \
  do {                                                                         \
    if ((FUNC) != 0) {                                                         \
      PyGILState_STATE gstate = PyGILState_Ensure();                           \
      PyErr_SetString(PyExc_OverflowError,                                     \
                      "Overflow occurred in npy_datetimestruct_to_datetime");  \
      PyGILState_Release(gstate);                                              \
      return -1;                                                               \
    }                                                                          \
  } while (0)

const int days_per_month_table[2][12] = {
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

/*
 * Returns 1 if the given year is a leap year, 0 otherwise.
 */
int is_leapyear(npy_int64 year) {
  return (year & 0x3) == 0 && /* year % 4 == 0 */
         ((year % 100) != 0 || (year % 400) == 0);
}

/*
 * Adjusts a datetimestruct based on a minutes offset. Assumes
 * the current values are valid.g
 */
void add_minutes_to_datetimestruct(npy_datetimestruct *dts, int minutes) {
  int isleap;

  /* MINUTES */
  dts->min += minutes;
  while (dts->min < 0) {
    dts->min += 60;
    dts->hour--;
  }
  while (dts->min >= 60) {
    dts->min -= 60;
    dts->hour++;
  }

  /* HOURS */
  while (dts->hour < 0) {
    dts->hour += 24;
    dts->day--;
  }
  while (dts->hour >= 24) {
    dts->hour -= 24;
    dts->day++;
  }

  /* DAYS */
  if (dts->day < 1) {
    dts->month--;
    if (dts->month < 1) {
      dts->year--;
      dts->month = 12;
    }
    isleap = is_leapyear(dts->year);
    dts->day += days_per_month_table[isleap][dts->month - 1];
  } else if (dts->day > 28) {
    isleap = is_leapyear(dts->year);
    if (dts->day > days_per_month_table[isleap][dts->month - 1]) {
      dts->day -= days_per_month_table[isleap][dts->month - 1];
      dts->month++;
      if (dts->month > 12) {
        dts->year++;
        dts->month = 1;
      }
    }
  }
}

/*
 * Calculates the days offset from the 1970 epoch.
 */
npy_int64 get_datetimestruct_days(const npy_datetimestruct *dts) {
  int i, month;
  npy_int64 year, days = 0;
  const int *month_lengths;

  year = dts->year - 1970;
  days = year * 365;

  /* Adjust for leap years */
  if (days >= 0) {
    /*
     * 1968 is the closest leap year before 1970.
     * Exclude the current year, so add 1.
     */
    year += 1;
    /* Add one day for each 4 years */
    days += year / 4;
    /* 1900 is the closest previous year divisible by 100 */
    year += 68;
    /* Subtract one day for each 100 years */
    days -= year / 100;
    /* 1600 is the closest previous year divisible by 400 */
    year += 300;
    /* Add one day for each 400 years */
    days += year / 400;
  } else {
    /*
     * 1972 is the closest later year after 1970.
     * Include the current year, so subtract 2.
     */
    year -= 2;
    /* Subtract one day for each 4 years */
    days += year / 4;
    /* 2000 is the closest later year divisible by 100 */
    year -= 28;
    /* Add one day for each 100 years */
    days -= year / 100;
    /* 2000 is also the closest later year divisible by 400 */
    /* Subtract one day for each 400 years */
    days += year / 400;
  }

  month_lengths = days_per_month_table[is_leapyear(dts->year)];
  month = dts->month - 1;

  /* Add the months */
  for (i = 0; i < month; ++i) {
    days += month_lengths[i];
  }

  /* Add the days */
  days += dts->day - 1;

  return days;
}

/*
 * Modifies '*days_' to be the day offset within the year,
 * and returns the year.
 */
static npy_int64 days_to_yearsdays(npy_int64 *days_) {
  const npy_int64 days_per_400years = (400 * 365 + 100 - 4 + 1);
  /* Adjust so it's relative to the year 2000 (divisible by 400) */
  npy_int64 days = (*days_) - (365 * 30 + 7);
  npy_int64 year;

  /* Break down the 400 year cycle to get the year and day within the year */
  if (days >= 0) {
    year = 400 * (days / days_per_400years);
    days = days % days_per_400years;
  } else {
    year = 400 * ((days - (days_per_400years - 1)) / days_per_400years);
    days = days % days_per_400years;
    if (days < 0) {
      days += days_per_400years;
    }
  }

  /* Work out the year/day within the 400 year cycle */
  if (days >= 366) {
    year += 100 * ((days - 1) / (100 * 365 + 25 - 1));
    days = (days - 1) % (100 * 365 + 25 - 1);
    if (days >= 365) {
      year += 4 * ((days + 1) / (4 * 365 + 1));
      days = (days + 1) % (4 * 365 + 1);
      if (days >= 366) {
        year += (days - 1) / 365;
        days = (days - 1) % 365;
      }
    }
  }

  *days_ = days;
  return year + 2000;
}

/*
 * Fills in the year, month, day in 'dts' based on the days
 * offset from 1970.
 */
static void set_datetimestruct_days(npy_int64 days, npy_datetimestruct *dts) {
  const int *month_lengths;
  int i;

  dts->year = days_to_yearsdays(&days);
  month_lengths = days_per_month_table[is_leapyear(dts->year)];

  for (i = 0; i < 12; ++i) {
    if (days < month_lengths[i]) {
      dts->month = i + 1;
      dts->day = (npy_int32)days + 1;
      return;
    } else {
      days -= month_lengths[i];
    }
  }
}

/*
 * Compares two npy_datetimestruct objects chronologically
 */
int cmp_npy_datetimestruct(const npy_datetimestruct *a,
                           const npy_datetimestruct *b) {
  if (a->year > b->year) {
    return 1;
  } else if (a->year < b->year) {
    return -1;
  }

  if (a->month > b->month) {
    return 1;
  } else if (a->month < b->month) {
    return -1;
  }

  if (a->day > b->day) {
    return 1;
  } else if (a->day < b->day) {
    return -1;
  }

  if (a->hour > b->hour) {
    return 1;
  } else if (a->hour < b->hour) {
    return -1;
  }

  if (a->min > b->min) {
    return 1;
  } else if (a->min < b->min) {
    return -1;
  }

  if (a->sec > b->sec) {
    return 1;
  } else if (a->sec < b->sec) {
    return -1;
  }

  if (a->us > b->us) {
    return 1;
  } else if (a->us < b->us) {
    return -1;
  }

  if (a->ps > b->ps) {
    return 1;
  } else if (a->ps < b->ps) {
    return -1;
  }

  if (a->as > b->as) {
    return 1;
  } else if (a->as < b->as) {
    return -1;
  }

  return 0;
}
/*
 * Returns the offset from utc of the timezone as a timedelta.
 * The caller is responsible for ensuring that the tzinfo
 * attribute exists on the datetime object.
 *
 * If the passed object is timezone naive, Py_None is returned.
 * If extraction of the offset fails, NULL is returned.
 *
 * NOTE: This function is not vendored from numpy.
 */
PyObject *extract_utc_offset(PyObject *obj) {
  PyObject *tmp = PyObject_GetAttrString(obj, "tzinfo");
  if (tmp == NULL) {
    return NULL;
  }
  if (tmp != Py_None) {
    PyObject *offset = PyObject_CallMethod(tmp, "utcoffset", "O", obj);
    if (offset == NULL) {
      Py_DECREF(tmp);
      return NULL;
    }
    return offset;
  }
  return tmp;
}

static inline int scaleYearToEpoch(int64_t year, int64_t *result) {
  return checked_int64_sub(year, 1970, result);
}

static inline int scaleYearsToMonths(int64_t years, int64_t *result) {
  return checked_int64_mul(years, 12, result);
}

static inline int scaleDaysToWeeks(int64_t days, int64_t *result) {
  if (days >= 0) {
    *result = days / 7;
    return 0;
  } else {
    int res;
    int64_t checked_days;
    if ((res = checked_int64_sub(days, 6, &checked_days))) {
      return res;
    }

    *result = checked_days / 7;
    return 0;
  }
}

static inline int scaleDaysToHours(int64_t days, int64_t *result) {
  return checked_int64_mul(days, 24, result);
}

static inline int scaleHoursToMinutes(int64_t hours, int64_t *result) {
  return checked_int64_mul(hours, 60, result);
}

static inline int scaleMinutesToSeconds(int64_t minutes, int64_t *result) {
  return checked_int64_mul(minutes, 60, result);
}

static inline int scaleSecondsToMilliseconds(int64_t seconds, int64_t *result) {
  return checked_int64_mul(seconds, 1000, result);
}

static inline int scaleSecondsToMicroseconds(int64_t seconds, int64_t *result) {
  return checked_int64_mul(seconds, 1000000, result);
}

static inline int scaleMicrosecondsToNanoseconds(int64_t microseconds,
                                                 int64_t *result) {
  return checked_int64_mul(microseconds, 1000, result);
}

static inline int scaleMicrosecondsToPicoseconds(int64_t microseconds,
                                                 int64_t *result) {
  return checked_int64_mul(microseconds, 1000000, result);
}

static inline int64_t scalePicosecondsToFemtoseconds(int64_t picoseconds,
                                                     int64_t *result) {
  return checked_int64_mul(picoseconds, 1000, result);
}

static inline int64_t scalePicosecondsToAttoseconds(int64_t picoseconds,
                                                    int64_t *result) {
  return checked_int64_mul(picoseconds, 1000000, result);
}

/*
 * Converts a datetime from a datetimestruct to a datetime based
 * on a metadata unit. Returns -1 on and sets PyErr on error.
 */
npy_datetime npy_datetimestruct_to_datetime(NPY_DATETIMEUNIT base,
                                            const npy_datetimestruct *dts) {
  if ((base == NPY_FR_Y) || (base == NPY_FR_M)) {
    int64_t years;
    PD_CHECK_OVERFLOW(scaleYearToEpoch(dts->year, &years));

    if (base == NPY_FR_Y) {
      return years;
    }

    int64_t months;
    PD_CHECK_OVERFLOW(scaleYearsToMonths(years, &months));

    int64_t months_adder;
    PD_CHECK_OVERFLOW(checked_int64_sub(dts->month, 1, &months_adder));
    PD_CHECK_OVERFLOW(checked_int64_add(months, months_adder, &months));

    if (base == NPY_FR_M) {
      return months;
    }
  }

  const int64_t days = get_datetimestruct_days(dts);
  if (base == NPY_FR_D) {
    return days;
  }

  if (base == NPY_FR_W) {
    int64_t weeks;
    PD_CHECK_OVERFLOW(scaleDaysToWeeks(days, &weeks));
    return weeks;
  }

  int64_t hours;
  PD_CHECK_OVERFLOW(scaleDaysToHours(days, &hours));
  PD_CHECK_OVERFLOW(checked_int64_add(hours, dts->hour, &hours));

  if (base == NPY_FR_h) {
    return hours;
  }

  int64_t minutes;
  PD_CHECK_OVERFLOW(scaleHoursToMinutes(hours, &minutes));
  PD_CHECK_OVERFLOW(checked_int64_add(minutes, dts->min, &minutes));

  if (base == NPY_FR_m) {
    return minutes;
  }

  int64_t seconds;
  PD_CHECK_OVERFLOW(scaleMinutesToSeconds(minutes, &seconds));
  PD_CHECK_OVERFLOW(checked_int64_add(seconds, dts->sec, &seconds));

  if (base == NPY_FR_s) {
    return seconds;
  }

  if (base == NPY_FR_ms) {
    int64_t milliseconds;
    PD_CHECK_OVERFLOW(scaleSecondsToMilliseconds(seconds, &milliseconds));
    PD_CHECK_OVERFLOW(
        checked_int64_add(milliseconds, dts->us / 1000, &milliseconds));

    return milliseconds;
  }

  int64_t microseconds;
  PD_CHECK_OVERFLOW(scaleSecondsToMicroseconds(seconds, &microseconds));
  PD_CHECK_OVERFLOW(checked_int64_add(microseconds, dts->us, &microseconds));

  if (base == NPY_FR_us) {
    return microseconds;
  }

  if (base == NPY_FR_ns) {
    int64_t nanoseconds;
    PD_CHECK_OVERFLOW(
        scaleMicrosecondsToNanoseconds(microseconds, &nanoseconds));
    PD_CHECK_OVERFLOW(
        checked_int64_add(nanoseconds, dts->ps / 1000, &nanoseconds));

    return nanoseconds;
  }

  int64_t picoseconds;
  PD_CHECK_OVERFLOW(scaleMicrosecondsToPicoseconds(microseconds, &picoseconds));
  PD_CHECK_OVERFLOW(checked_int64_add(picoseconds, dts->ps, &picoseconds));

  if (base == NPY_FR_ps) {
    return picoseconds;
  }

  if (base == NPY_FR_fs) {
    int64_t femtoseconds;
    PD_CHECK_OVERFLOW(
        scalePicosecondsToFemtoseconds(picoseconds, &femtoseconds));
    PD_CHECK_OVERFLOW(
        checked_int64_add(femtoseconds, dts->as / 1000, &femtoseconds));
    return femtoseconds;
  }

  if (base == NPY_FR_as) {
    int64_t attoseconds;
    PD_CHECK_OVERFLOW(scalePicosecondsToAttoseconds(picoseconds, &attoseconds));
    PD_CHECK_OVERFLOW(checked_int64_add(attoseconds, dts->as, &attoseconds));
    return attoseconds;
  }

  /* Something got corrupted */
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyErr_SetString(PyExc_ValueError,
                  "NumPy datetime metadata with corrupt unit value");
  PyGILState_Release(gstate);

  return -1;
}

/*
 * Port numpy#13188 https://github.com/numpy/numpy/pull/13188/
 *
 * Computes the python `ret, d = divmod(d, unit)`.
 *
 * Note that GCC is smart enough at -O2 to eliminate the `if(*d < 0)` branch
 * for subsequent calls to this command - it is able to deduce that `*d >= 0`.
 */
npy_int64 extract_unit(npy_datetime *d, npy_datetime unit) {
  assert(unit > 0);
  npy_int64 div = *d / unit;
  npy_int64 mod = *d % unit;
  if (mod < 0) {
    mod += unit;
    div -= 1;
  }
  assert(mod >= 0);
  *d = mod;
  return div;
}

/*
 * Converts a datetime based on the given metadata into a datetimestruct
 */
void pandas_datetime_to_datetimestruct(npy_datetime dt, NPY_DATETIMEUNIT base,
                                       npy_datetimestruct *out) {
  npy_int64 perday;

  /* Initialize the output to all zeros */
  memset(out, 0, sizeof(npy_datetimestruct));
  out->year = 1970;
  out->month = 1;
  out->day = 1;

  /*
   * Note that care must be taken with the / and % operators
   * for negative values.
   */
  switch (base) {
  case NPY_FR_Y:
    out->year = 1970 + dt;
    break;

  case NPY_FR_M:
    out->year = 1970 + extract_unit(&dt, 12);
    out->month = (npy_int32)dt + 1;
    break;

  case NPY_FR_W:
    /* A week is 7 days */
    set_datetimestruct_days(dt * 7, out);
    break;

  case NPY_FR_D:
    set_datetimestruct_days(dt, out);
    break;

  case NPY_FR_h:
    perday = 24LL;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)dt;
    break;

  case NPY_FR_m:
    perday = 24LL * 60;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)extract_unit(&dt, 60);
    out->min = (npy_int32)dt;
    break;

  case NPY_FR_s:
    perday = 24LL * 60 * 60;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)extract_unit(&dt, 60 * 60);
    out->min = (npy_int32)extract_unit(&dt, 60);
    out->sec = (npy_int32)dt;
    break;

  case NPY_FR_ms:
    perday = 24LL * 60 * 60 * 1000;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)extract_unit(&dt, 1000LL * 60 * 60);
    out->min = (npy_int32)extract_unit(&dt, 1000LL * 60);
    out->sec = (npy_int32)extract_unit(&dt, 1000LL);
    out->us = (npy_int32)(dt * 1000);
    break;

  case NPY_FR_us:
    perday = 24LL * 60LL * 60LL * 1000LL * 1000LL;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 60 * 60);
    out->min = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 60);
    out->sec = (npy_int32)extract_unit(&dt, 1000LL * 1000);
    out->us = (npy_int32)dt;
    break;

  case NPY_FR_ns:
    perday = 24LL * 60LL * 60LL * 1000LL * 1000LL * 1000LL;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 60 * 60);
    out->min = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 60);
    out->sec = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000);
    out->us = (npy_int32)extract_unit(&dt, 1000LL);
    out->ps = (npy_int32)(dt * 1000);
    break;

  case NPY_FR_ps:
    perday = 24LL * 60 * 60 * 1000 * 1000 * 1000 * 1000;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 60 * 60);
    out->min = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 60);
    out->sec = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000);
    out->us = (npy_int32)extract_unit(&dt, 1000LL);
    out->ps = (npy_int32)(dt * 1000);
    break;

  case NPY_FR_fs:
    /* entire range is only +- 2.6 hours */
    out->hour = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 *
                                                 1000 * 60 * 60);
    if (out->hour < 0) {
      out->year = 1969;
      out->month = 12;
      out->day = 31;
      out->hour += 24;
      assert(out->hour >= 0);
    }
    out->min =
        (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000 * 60);
    out->sec = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000);
    out->us = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000);
    out->ps = (npy_int32)extract_unit(&dt, 1000LL);
    out->as = (npy_int32)(dt * 1000);
    break;

  case NPY_FR_as:
    /* entire range is only +- 9.2 seconds */
    out->sec =
        (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000 * 1000);
    if (out->sec < 0) {
      out->year = 1969;
      out->month = 12;
      out->day = 31;
      out->hour = 23;
      out->min = 59;
      out->sec += 60;
      assert(out->sec >= 0);
    }
    out->us = (npy_int32)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000);
    out->ps = (npy_int32)extract_unit(&dt, 1000LL * 1000);
    out->as = (npy_int32)dt;
    break;

  default:
    PyErr_SetString(PyExc_RuntimeError,
                    "NumPy datetime metadata is corrupted with invalid "
                    "base unit");
  }
}

/*
 * Converts a timedelta from a timedeltastruct to a timedelta based
 * on a metadata unit. The timedelta is assumed to be valid.
 *
 * Returns 0 on success, -1 on failure.
 */
void pandas_timedelta_to_timedeltastruct(npy_timedelta td,
                                         NPY_DATETIMEUNIT base,
                                         pandas_timedeltastruct *out) {
  npy_int64 frac;
  npy_int64 sfrac;
  npy_int64 ifrac;
  int sign;
  npy_int64 per_day;
  npy_int64 per_sec;

  /* Initialize the output to all zeros */
  memset(out, 0, sizeof(pandas_timedeltastruct));

  switch (base) {
  case NPY_FR_ns:

    per_day = 86400000000000LL;
    per_sec = 1000LL * 1000LL * 1000LL;

    // put frac in seconds
    if (td < 0 && td % per_sec != 0)
      frac = td / per_sec - 1;
    else
      frac = td / per_sec;

    if (frac < 0) {
      sign = -1;

      // even fraction
      if ((-frac % 86400LL) != 0) {
        out->days = -frac / 86400LL + 1;
        frac += 86400LL * out->days;
      } else {
        frac = -frac;
      }
    } else {
      sign = 1;
      out->days = 0;
    }

    if (frac >= 86400) {
      out->days += frac / 86400LL;
      frac -= out->days * 86400LL;
    }

    if (frac >= 3600) {
      out->hrs = (npy_int32)(frac / 3600LL);
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = (npy_int32)(frac / 60LL);
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = (npy_int32)frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = (npy_int32)(ifrac / (1000LL * 1000LL));
      ifrac -= out->ms * 1000LL * 1000LL;
      out->us = (npy_int32)(ifrac / 1000LL);
      ifrac -= out->us * 1000LL;
      out->ns = (npy_int32)ifrac;
    } else {
      out->ms = 0;
      out->us = 0;
      out->ns = 0;
    }
    break;

  case NPY_FR_us:

    per_day = 86400000000LL;
    per_sec = 1000LL * 1000LL;

    // put frac in seconds
    if (td < 0 && td % per_sec != 0)
      frac = td / per_sec - 1;
    else
      frac = td / per_sec;

    if (frac < 0) {
      sign = -1;

      // even fraction
      if ((-frac % 86400LL) != 0) {
        out->days = -frac / 86400LL + 1;
        frac += 86400LL * out->days;
      } else {
        frac = -frac;
      }
    } else {
      sign = 1;
      out->days = 0;
    }

    if (frac >= 86400) {
      out->days += frac / 86400LL;
      frac -= out->days * 86400LL;
    }

    if (frac >= 3600) {
      out->hrs = (npy_int32)(frac / 3600LL);
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = (npy_int32)(frac / 60LL);
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = (npy_int32)frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = (npy_int32)(ifrac / 1000LL);
      ifrac -= out->ms * 1000LL;
      out->us = (npy_int32)(ifrac / 1L);
      ifrac -= out->us * 1L;
      out->ns = (npy_int32)ifrac;
    } else {
      out->ms = 0;
      out->us = 0;
      out->ns = 0;
    }
    break;

  case NPY_FR_ms:

    per_day = 86400000LL;
    per_sec = 1000LL;

    // put frac in seconds
    if (td < 0 && td % per_sec != 0)
      frac = td / per_sec - 1;
    else
      frac = td / per_sec;

    if (frac < 0) {
      sign = -1;

      // even fraction
      if ((-frac % 86400LL) != 0) {
        out->days = -frac / 86400LL + 1;
        frac += 86400LL * out->days;
      } else {
        frac = -frac;
      }
    } else {
      sign = 1;
      out->days = 0;
    }

    if (frac >= 86400) {
      out->days += frac / 86400LL;
      frac -= out->days * 86400LL;
    }

    if (frac >= 3600) {
      out->hrs = (npy_int32)(frac / 3600LL);
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = (npy_int32)(frac / 60LL);
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = (npy_int32)frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = (npy_int32)ifrac;
      out->us = 0;
      out->ns = 0;
    } else {
      out->ms = 0;
      out->us = 0;
      out->ns = 0;
    }
    break;

  case NPY_FR_s:
    // special case where we can simplify many expressions bc per_sec=1

    per_day = 86400LL;
    per_sec = 1L;

    // put frac in seconds
    if (td < 0 && td % per_sec != 0)
      frac = td / per_sec - 1;
    else
      frac = td / per_sec;

    if (frac < 0) {
      sign = -1;

      // even fraction
      if ((-frac % 86400LL) != 0) {
        out->days = -frac / 86400LL + 1;
        frac += 86400LL * out->days;
      } else {
        frac = -frac;
      }
    } else {
      sign = 1;
      out->days = 0;
    }

    if (frac >= 86400) {
      out->days += frac / 86400LL;
      frac -= out->days * 86400LL;
    }

    if (frac >= 3600) {
      out->hrs = (npy_int32)(frac / 3600LL);
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = (npy_int32)(frac / 60LL);
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = (npy_int32)frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = 0;
      out->us = 0;
      out->ns = 0;
    } else {
      out->ms = 0;
      out->us = 0;
      out->ns = 0;
    }
    break;

  case NPY_FR_m:

    out->days = td / 1440LL;
    td -= out->days * 1440LL;
    out->hrs = (npy_int32)(td / 60LL);
    td -= out->hrs * 60LL;
    out->min = (npy_int32)td;

    out->sec = 0;
    out->ms = 0;
    out->us = 0;
    out->ns = 0;
    break;

  case NPY_FR_h:
    out->days = td / 24LL;
    td -= out->days * 24LL;
    out->hrs = (npy_int32)td;

    out->min = 0;
    out->sec = 0;
    out->ms = 0;
    out->us = 0;
    out->ns = 0;
    break;

  case NPY_FR_D:
    out->days = td;
    out->hrs = 0;
    out->min = 0;
    out->sec = 0;
    out->ms = 0;
    out->us = 0;
    out->ns = 0;
    break;

  case NPY_FR_W:
    out->days = 7 * td;
    out->hrs = 0;
    out->min = 0;
    out->sec = 0;
    out->ms = 0;
    out->us = 0;
    out->ns = 0;
    break;

  default:
    PyErr_SetString(PyExc_RuntimeError,
                    "NumPy timedelta metadata is corrupted with "
                    "invalid base unit");
  }

  out->seconds = out->hrs * 3600 + out->min * 60 + out->sec;
  out->microseconds = out->ms * 1000 + out->us;
  out->nanoseconds = out->ns;
}

/*
 * This function returns a pointer to the DateTimeMetaData
 * contained within the provided datetime dtype.
 *
 * Copied near-verbatim from numpy/core/src/multiarray/datetime.c
 */
PyArray_DatetimeMetaData
get_datetime_metadata_from_dtype(PyArray_Descr *dtype) {
  return (((PyArray_DatetimeDTypeMetaData *)dtype->c_metadata)->meta);
}
