/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

This file is derived from NumPy 1.7. See NUMPY_LICENSE.txt

*/

#define _PANDAS_DATETIME_IMPL

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "datetime.h"
#include "pd_datetime.h"

static const int days_per_month_table[2][12] = {
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

/*
 * Returns 1 if the given year is a leap year, 0 otherwise.
 */
static int is_leapyear(npy_int64 year) {
  return (year & 0x3) == 0 && /* year % 4 == 0 */
         ((year % 100) != 0 || (year % 400) == 0);
}

/*
 * Calculates the days offset from the 1970 epoch.
 */
static npy_int64 get_datetimestruct_days(const npy_datetimestruct *dts) {
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
 * Converts a datetime from a datetimestruct to a datetime based
 * on a metadata unit. The date is assumed to be valid.
 */
static npy_datetime
npy_datetimestruct_to_datetime(NPY_DATETIMEUNIT base,
                               const npy_datetimestruct *dts) {
  npy_datetime ret;

  if (base == NPY_FR_Y) {
    /* Truncate to the year */
    ret = dts->year - 1970;
  } else if (base == NPY_FR_M) {
    /* Truncate to the month */
    ret = 12 * (dts->year - 1970) + (dts->month - 1);
  } else {
    /* Otherwise calculate the number of days to start */
    npy_int64 days = get_datetimestruct_days(dts);

    switch (base) {
    case NPY_FR_W:
      /* Truncate to weeks */
      if (days >= 0) {
        ret = days / 7;
      } else {
        ret = (days - 6) / 7;
      }
      break;
    case NPY_FR_D:
      ret = days;
      break;
    case NPY_FR_h:
      ret = days * 24 + dts->hour;
      break;
    case NPY_FR_m:
      ret = (days * 24 + dts->hour) * 60 + dts->min;
      break;
    case NPY_FR_s:
      ret = ((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec;
      break;
    case NPY_FR_ms:
      ret = (((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec) * 1000 +
            dts->us / 1000;
      break;
    case NPY_FR_us:
      ret = (((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec) *
                1000000 +
            dts->us;
      break;
    case NPY_FR_ns:
      ret = ((((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec) *
                 1000000 +
             dts->us) *
                1000 +
            dts->ps / 1000;
      break;
    case NPY_FR_ps:
      ret = ((((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec) *
                 1000000 +
             dts->us) *
                1000000 +
            dts->ps;
      break;
    case NPY_FR_fs:
      /* only 2.6 hours */
      ret = (((((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec) *
                  1000000 +
              dts->us) *
                 1000000 +
             dts->ps) *
                1000 +
            dts->as / 1000;
      break;
    case NPY_FR_as:
      /* only 9.2 secs */
      ret = (((((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec) *
                  1000000 +
              dts->us) *
                 1000000 +
             dts->ps) *
                1000000 +
            dts->as;
      break;
    default:
      /* Something got corrupted */
      PyErr_SetString(PyExc_ValueError,
                      "NumPy datetime metadata with corrupt unit value");
      return -1;
    }
  }
  return ret;
}

/*
 * Function: scaleNanosecToUnit
 * -----------------------------
 *
 * Scales an integer value representing time in nanoseconds to provided unit.
 *
 * Mutates the provided value directly. Returns 0 on success, non-zero on
 * error.
 */
static int scaleNanosecToUnit(npy_int64 *value, NPY_DATETIMEUNIT unit) {
  switch (unit) {
  case NPY_FR_ns:
    break;
  case NPY_FR_us:
    *value /= 1000LL;
    break;
  case NPY_FR_ms:
    *value /= 1000000LL;
    break;
  case NPY_FR_s:
    *value /= 1000000000LL;
    break;
  default:
    return -1;
  }

  return 0;
}

/*
 * Compares two npy_datetimestruct objects chronologically
 */
static int cmp_npy_datetimestruct(const npy_datetimestruct *a,
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
static PyObject *extract_utc_offset(PyObject *obj) {
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
 *
 * Converts a Python datetime.datetime or datetime.date
 * object into a NumPy npy_datetimestruct.  Uses tzinfo (if present)
 * to convert to UTC time.
 *
 * The following implementation just asks for attributes, and thus
 * supports datetime duck typing. The tzinfo time zone conversion
 * requires this style of access as well.
 *
 * Returns -1 on error, 0 on success, and 1 (with no error set)
 * if obj doesn't have the needed date or datetime attributes.
 */
static int convert_pydatetime_to_datetimestruct(PyObject *dtobj,
                                                npy_datetimestruct *out) {
  // Assumes that obj is a valid datetime object
  PyObject *tmp;
  PyObject *obj = (PyObject *)dtobj;

  /* Initialize the output to all zeros */
  memset(out, 0, sizeof(npy_datetimestruct));
  out->month = 1;
  out->day = 1;

  out->year = PyLong_AsLong(PyObject_GetAttrString(obj, "year"));
  out->month = PyLong_AsLong(PyObject_GetAttrString(obj, "month"));
  out->day = PyLong_AsLong(PyObject_GetAttrString(obj, "day"));

  // TODO(anyone): If we can get PyDateTime_IMPORT to work, we could use
  // PyDateTime_Check here, and less verbose attribute lookups.

  /* Check for time attributes (if not there, return success as a date) */
  if (!PyObject_HasAttrString(obj, "hour") ||
      !PyObject_HasAttrString(obj, "minute") ||
      !PyObject_HasAttrString(obj, "second") ||
      !PyObject_HasAttrString(obj, "microsecond")) {
    return 0;
  }

  out->hour = PyLong_AsLong(PyObject_GetAttrString(obj, "hour"));
  out->min = PyLong_AsLong(PyObject_GetAttrString(obj, "minute"));
  out->sec = PyLong_AsLong(PyObject_GetAttrString(obj, "second"));
  out->us = PyLong_AsLong(PyObject_GetAttrString(obj, "microsecond"));

  if (PyObject_HasAttrString(obj, "tzinfo")) {
    PyObject *offset = extract_utc_offset(obj);
    /* Apply the time zone offset if datetime obj is tz-aware */
    if (offset != NULL) {
      if (offset == Py_None) {
        Py_DECREF(offset);
        return 0;
      }
      PyObject *tmp_int;
      int seconds_offset, minutes_offset;
      /*
       * The timedelta should have a function "total_seconds"
       * which contains the value we want.
       */
      tmp = PyObject_CallMethod(offset, "total_seconds", "");
      Py_DECREF(offset);
      if (tmp == NULL) {
        return -1;
      }
      tmp_int = PyNumber_Long(tmp);
      if (tmp_int == NULL) {
        Py_DECREF(tmp);
        return -1;
      }
      seconds_offset = PyLong_AsLong(tmp_int);
      if (seconds_offset == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp_int);
        Py_DECREF(tmp);
        return -1;
      }
      Py_DECREF(tmp_int);
      Py_DECREF(tmp);

      /* Convert to a minutes offset and apply it */
      minutes_offset = seconds_offset / 60;

      add_minutes_to_datetimestruct(out, -minutes_offset);
    }
  }

  return 0;
}

/*
 * Port numpy#13188 https://github.com/numpy/numpy/pull/13188/
 *
 * Computes the python `ret, d = divmod(d, unit)`.
 *
 * Note that GCC is smart enough at -O2 to eliminate the `if(*d < 0)` branch
 * for subsequent calls to this command - it is able to deduce that `*d >= 0`.
 */
static npy_int64 extract_unit(npy_datetime *d, npy_datetime unit) {
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
      dts->day = days + 1;
      return;
    } else {
      days -= month_lengths[i];
    }
  }
}

/*
 * Converts a datetime based on the given metadata into a datetimestruct
 */
static void pandas_datetime_to_datetimestruct(npy_datetime dt,
                                              NPY_DATETIMEUNIT base,
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
    out->month = dt + 1;
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
    out->hour = dt;
    break;

  case NPY_FR_m:
    perday = 24LL * 60;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (int)extract_unit(&dt, 60);
    out->min = (int)dt;
    break;

  case NPY_FR_s:
    perday = 24LL * 60 * 60;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (int)extract_unit(&dt, 60 * 60);
    out->min = (int)extract_unit(&dt, 60);
    out->sec = (int)dt;
    break;

  case NPY_FR_ms:
    perday = 24LL * 60 * 60 * 1000;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (int)extract_unit(&dt, 1000LL * 60 * 60);
    out->min = (int)extract_unit(&dt, 1000LL * 60);
    out->sec = (int)extract_unit(&dt, 1000LL);
    out->us = (int)(dt * 1000);
    break;

  case NPY_FR_us:
    perday = 24LL * 60LL * 60LL * 1000LL * 1000LL;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (int)extract_unit(&dt, 1000LL * 1000 * 60 * 60);
    out->min = (int)extract_unit(&dt, 1000LL * 1000 * 60);
    out->sec = (int)extract_unit(&dt, 1000LL * 1000);
    out->us = (int)dt;
    break;

  case NPY_FR_ns:
    perday = 24LL * 60LL * 60LL * 1000LL * 1000LL * 1000LL;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 60 * 60);
    out->min = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 60);
    out->sec = (int)extract_unit(&dt, 1000LL * 1000 * 1000);
    out->us = (int)extract_unit(&dt, 1000LL);
    out->ps = (int)(dt * 1000);
    break;

  case NPY_FR_ps:
    perday = 24LL * 60 * 60 * 1000 * 1000 * 1000 * 1000;

    set_datetimestruct_days(extract_unit(&dt, perday), out);
    out->hour = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 60 * 60);
    out->min = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 60);
    out->sec = (int)extract_unit(&dt, 1000LL * 1000 * 1000);
    out->us = (int)extract_unit(&dt, 1000LL);
    out->ps = (int)(dt * 1000);
    break;

  case NPY_FR_fs:
    /* entire range is only +- 2.6 hours */
    out->hour =
        (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000 * 60 * 60);
    if (out->hour < 0) {
      out->year = 1969;
      out->month = 12;
      out->day = 31;
      out->hour += 24;
      assert(out->hour >= 0);
    }
    out->min = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000 * 60);
    out->sec = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000);
    out->us = (int)extract_unit(&dt, 1000LL * 1000 * 1000);
    out->ps = (int)extract_unit(&dt, 1000LL);
    out->as = (int)(dt * 1000);
    break;

  case NPY_FR_as:
    /* entire range is only +- 9.2 seconds */
    out->sec =
        (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000 * 1000 * 1000);
    if (out->sec < 0) {
      out->year = 1969;
      out->month = 12;
      out->day = 31;
      out->hour = 23;
      out->min = 59;
      out->sec += 60;
      assert(out->sec >= 0);
    }
    out->us = (int)extract_unit(&dt, 1000LL * 1000 * 1000 * 1000);
    out->ps = (int)extract_unit(&dt, 1000LL * 1000);
    out->as = (int)dt;
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
static void pandas_timedelta_to_timedeltastruct(npy_timedelta td,
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
      out->hrs = frac / 3600LL;
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = frac / 60LL;
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = ifrac / (1000LL * 1000LL);
      ifrac -= out->ms * 1000LL * 1000LL;
      out->us = ifrac / 1000LL;
      ifrac -= out->us * 1000LL;
      out->ns = ifrac;
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
      out->hrs = frac / 3600LL;
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = frac / 60LL;
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = ifrac / 1000LL;
      ifrac -= out->ms * 1000LL;
      out->us = ifrac / 1L;
      ifrac -= out->us * 1L;
      out->ns = ifrac;
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
      out->hrs = frac / 3600LL;
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = frac / 60LL;
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = frac;
      frac -= out->sec;
    } else {
      out->sec = 0;
    }

    sfrac = (out->hrs * 3600LL + out->min * 60LL + out->sec) * per_sec;

    if (sign < 0)
      out->days = -out->days;

    ifrac = td - (out->days * per_day + sfrac);

    if (ifrac != 0) {
      out->ms = ifrac;
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
      out->hrs = frac / 3600LL;
      frac -= out->hrs * 3600LL;
    } else {
      out->hrs = 0;
    }

    if (frac >= 60) {
      out->min = frac / 60LL;
      frac -= out->min * 60LL;
    } else {
      out->min = 0;
    }

    if (frac >= 0) {
      out->sec = frac;
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
    out->hrs = td / 60LL;
    td -= out->hrs * 60LL;
    out->min = td;

    out->sec = 0;
    out->ms = 0;
    out->us = 0;
    out->ns = 0;
    break;

  case NPY_FR_h:
    out->days = td / 24LL;
    td -= out->days * 24LL;
    out->hrs = td;

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
 * Provides a string length to use for converting datetime
 * objects with the given local and unit settings.
 */
static int get_datetime_iso_8601_strlen(int local, NPY_DATETIMEUNIT base) {
  int len = 0;

  switch (base) {
  /* Generic units can only be used to represent NaT */
  /*    return 4;*/
  case NPY_FR_as:
    len += 3; /* "###" */
  case NPY_FR_fs:
    len += 3; /* "###" */
  case NPY_FR_ps:
    len += 3; /* "###" */
  case NPY_FR_ns:
    len += 3; /* "###" */
  case NPY_FR_us:
    len += 3; /* "###" */
  case NPY_FR_ms:
    len += 4; /* ".###" */
  case NPY_FR_s:
    len += 3; /* ":##" */
  case NPY_FR_m:
    len += 3; /* ":##" */
  case NPY_FR_h:
    len += 3; /* "T##" */
  case NPY_FR_D:
  case NPY_FR_W:
    len += 3; /* "-##" */
  case NPY_FR_M:
    len += 3; /* "-##" */
  case NPY_FR_Y:
    len += 21; /* 64-bit year */
    break;
  default:
    len += 3; /* handle the now defunct NPY_FR_B */
    break;
  }

  if (base >= NPY_FR_h) {
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
 * Converts an npy_datetimestruct to an (almost) ISO 8601
 * NULL-terminated string using timezone Z (UTC). If the string fits in
 * the space exactly, it leaves out the NULL terminator and returns success.
 *
 * The differences from ISO 8601 are the 'NaT' string, and
 * the number of year digits is >= 4 instead of strictly 4.
 *
 * 'base' restricts the output to that unit. Set 'base' to
 * -1 to auto-detect a base after which all the values are zero.
 *
 *  Returns 0 on success, -1 on failure (for example if the output
 *  string was too short).
 */
static int make_iso_8601_datetime(npy_datetimestruct *dts, char *outstr,
                                  int outlen, int utc, NPY_DATETIMEUNIT base) {
  char *substr = outstr;
  int sublen = outlen;
  int tmplen;

  /*
   * Print weeks with the same precision as days.
   *
   * TODO: Could print weeks with YYYY-Www format if the week
   *       epoch is a Monday.
   */
  if (base == NPY_FR_W) {
    base = NPY_FR_D;
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
#endif  // _WIN32
  /* If it ran out of space or there isn't space for the NULL terminator */
  if (tmplen < 0 || tmplen > sublen) {
    goto string_too_short;
  }
  substr += tmplen;
  sublen -= tmplen;

  /* Stop if the unit is years */
  if (base == NPY_FR_Y) {
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
  if (base == NPY_FR_M) {
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
  if (base == NPY_FR_D) {
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
  if (base == NPY_FR_h) {
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
  if (base == NPY_FR_m) {
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
  if (base == NPY_FR_s) {
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
  if (base == NPY_FR_ms) {
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
  if (base == NPY_FR_us) {
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
  if (base == NPY_FR_ns) {
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
  if (base == NPY_FR_ps) {
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
  if (base == NPY_FR_fs) {
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
  /* UTC "Zulu" time */
  if (utc) {
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

static int make_iso_8601_timedelta(pandas_timedeltastruct *tds, char *outstr,
                                   size_t *outlen) {
  *outlen = 0;
  *outlen += snprintf(outstr, 60, // NOLINT
                      "P%" NPY_INT64_FMT "DT%" NPY_INT32_FMT "H%" NPY_INT32_FMT
                      "M%" NPY_INT32_FMT,
                      tds->days, tds->hrs, tds->min, tds->sec);
  outstr += *outlen;

  if (tds->ns != 0) {
    *outlen += snprintf(outstr, 12, // NOLINT
                        ".%03" NPY_INT32_FMT "%03" NPY_INT32_FMT
                        "%03" NPY_INT32_FMT "S",
                        tds->ms, tds->us, tds->ns);
  } else if (tds->us != 0) {
    *outlen += snprintf(outstr, 9, // NOLINT
                        ".%03" NPY_INT32_FMT "%03" NPY_INT32_FMT "S", tds->ms,
                        tds->us);
  } else if (tds->ms != 0) {
    *outlen += snprintf(outstr, 6, // NOLINT
                        ".%03" NPY_INT32_FMT "S", tds->ms);
  } else {
    *outlen += snprintf(outstr, 2, // NOLINT
                        "%s", "S");
  }

  return 0;
}

/* Converts the int64_t representation of a datetime to ISO; mutates len */
static char *int64ToIso(int64_t value, NPY_DATETIMEUNIT base, size_t *len) {
  npy_datetimestruct dts;
  int ret_code;

  pandas_datetime_to_datetimestruct(value, NPY_FR_ns, &dts);

  *len = (size_t)get_datetime_iso_8601_strlen(0, base);
  char *result = PyObject_Malloc(*len);

  if (result == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  // datetime64 is always naive
  ret_code = make_iso_8601_datetime(&dts, result, *len, 0, base);
  if (ret_code != 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Could not convert datetime value to string");
    PyObject_Free(result);
  }

  // Note that get_datetime_iso_8601_strlen just gives a generic size
  // for ISO string conversion, not the actual size used
  *len = strlen(result);
  return result;
}

static npy_datetime NpyDateTimeToEpoch(npy_datetime dt, NPY_DATETIMEUNIT base) {
  scaleNanosecToUnit(&dt, base);
  return dt;
}

/* Convert PyDatetime To ISO C-string. mutates len */
static char *PyDateTimeToIso(PyObject *obj, NPY_DATETIMEUNIT base,
                             size_t *len) {
  npy_datetimestruct dts;
  int ret;

  ret = convert_pydatetime_to_datetimestruct(obj, &dts);
  if (ret != 0) {
    if (!PyErr_Occurred()) {
      PyErr_SetString(PyExc_ValueError,
                      "Could not convert PyDateTime to numpy datetime");
    }
    return NULL;
  }

  *len = (size_t)get_datetime_iso_8601_strlen(0, base);
  char *result = PyObject_Malloc(*len);
  // Check to see if PyDateTime has a timezone.
  // Don't convert to UTC if it doesn't.
  int is_tz_aware = 0;
  if (PyObject_HasAttrString(obj, "tzinfo")) {
    PyObject *offset = extract_utc_offset(obj);
    if (offset == NULL) {
      PyObject_Free(result);
      return NULL;
    }
    is_tz_aware = offset != Py_None;
    Py_DECREF(offset);
  }
  ret = make_iso_8601_datetime(&dts, result, *len, is_tz_aware, base);

  if (ret != 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Could not convert datetime value to string");
    PyObject_Free(result);
    return NULL;
  }

  // Note that get_datetime_iso_8601_strlen just gives a generic size
  // for ISO string conversion, not the actual size used
  *len = strlen(result);
  return result;
}

static npy_datetime PyDateTimeToEpoch(PyObject *dt, NPY_DATETIMEUNIT base) {
  npy_datetimestruct dts;
  int ret;

  ret = convert_pydatetime_to_datetimestruct(dt, &dts);
  if (ret != 0) {
    if (!PyErr_Occurred()) {
      PyErr_SetString(PyExc_ValueError,
                      "Could not convert PyDateTime to numpy datetime");
    }
    // TODO(username): is setting errMsg required?
    // ((JSONObjectEncoder *)tc->encoder)->errorMsg = "";
    // return NULL;
  }

  npy_datetime npy_dt = npy_datetimestruct_to_datetime(NPY_FR_ns, &dts);
  return NpyDateTimeToEpoch(npy_dt, base);
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
 *
 * 'out' gets filled with the parsed date-time.
 * 'out_local' gets set to 1 if the parsed time contains timezone,
 *      to 0 otherwise.
 * 'out_tzoffset' gets set to timezone offset by minutes
 *      if the parsed time was in local time,
 *      to 0 otherwise. The values 'now' and 'today' don't get counted
 *      as local, and neither do UTC +/-#### timezone offsets, because
 *      they aren't using the computer's local timezone offset.
 *
 * Returns 0 on success, -1 on failure.
 */

typedef enum {
  COMPARISON_SUCCESS,
  COMPLETED_PARTIAL_MATCH,
  COMPARISON_ERROR
} DatetimePartParseResult;
// This function will advance the pointer on format
// and decrement characters_remaining by n on success
// On failure will return COMPARISON_ERROR without incrementing
// If `format_requirement` is PARTIAL_MATCH, and the `format` string has
// been exhausted, then return COMPLETED_PARTIAL_MATCH.
static DatetimePartParseResult
compare_format(const char **format, int *characters_remaining,
               const char *compare_to, int n,
               const FormatRequirement format_requirement) {
  if (format_requirement == INFER_FORMAT) {
    return COMPARISON_SUCCESS;
  }
  if (*characters_remaining < 0) {
    return COMPARISON_ERROR;
  }
  if (format_requirement == PARTIAL_MATCH && *characters_remaining == 0) {
    return COMPLETED_PARTIAL_MATCH;
  }
  if (*characters_remaining < n) {
    // TODO(pandas-dev): PyErr to differentiate what went wrong
    return COMPARISON_ERROR;
  } else {
    if (strncmp(*format, compare_to, n)) {
      // TODO(pandas-dev): PyErr to differentiate what went wrong
      return COMPARISON_ERROR;
    } else {
      *format += n;
      *characters_remaining -= n;
      return COMPARISON_SUCCESS;
    }
  }
  return COMPARISON_SUCCESS;
}

static int parse_iso_8601_datetime(const char *str, int len, int want_exc,
                                   npy_datetimestruct *out,
                                   NPY_DATETIMEUNIT *out_bestunit,
                                   int *out_local, int *out_tzoffset,
                                   const char *format, int format_len,
                                   FormatRequirement format_requirement) {
  if (len < 0 || format_len < 0)
    goto parse_error;
  int year_leap = 0;
  int i, numdigits;
  const char *substr;
  int sublen;
  NPY_DATETIMEUNIT bestunit = NPY_FR_GENERIC;
  DatetimePartParseResult comparison;

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
  memset(out, 0, sizeof(npy_datetimestruct));
  out->month = 1;
  out->day = 1;

  substr = str;
  sublen = len;

  /* Skip leading whitespace */
  while (sublen > 0 && isspace(*substr)) {
    ++substr;
    --sublen;
    comparison =
        compare_format(&format, &format_len, " ", 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
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
  comparison =
      compare_format(&format, &format_len, "%Y", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }

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
    if (format_len) {
      goto parse_error;
    }
    bestunit = NPY_FR_Y;
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

    comparison =
        compare_format(&format, &format_len, &ymd_sep, 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
    /* Cannot have trailing separator */
    if (sublen == 0 || !isdigit(*substr)) {
      goto parse_error;
    }
  }

  /* PARSE THE MONTH */
  comparison =
      compare_format(&format, &format_len, "%m", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
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
    if (want_exc) {
      PyErr_Format(PyExc_ValueError,
                   "Month out of range in datetime string \"%s\"", str);
    }
    goto error;
  }

  /* Next character must be the separator, start of day, or end of string */
  if (sublen == 0) {
    bestunit = NPY_FR_M;
    /* Forbid YYYYMM. Parsed instead as YYMMDD by someone else. */
    if (!has_ymd_sep) {
      goto parse_error;
    }
    if (format_len) {
      goto parse_error;
    }
    if (out_local != NULL) {
      *out_local = 0;
    }
    goto finish;
  }

  if (has_ymd_sep) {
    /* Must have separator, but cannot be trailing */
    if (*substr != ymd_sep || sublen == 1) {
      goto parse_error;
    }
    ++substr;
    --sublen;
    comparison =
        compare_format(&format, &format_len, &ymd_sep, 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
  }

  /* PARSE THE DAY */
  comparison =
      compare_format(&format, &format_len, "%d", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
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
    if (want_exc) {
      PyErr_Format(PyExc_ValueError,
                   "Day out of range in datetime string \"%s\"", str);
    }
    goto error;
  }

  /* Next character must be a 'T', ' ', or end of string */
  if (sublen == 0) {
    if (out_local != NULL) {
      *out_local = 0;
    }
    if (format_len) {
      goto parse_error;
    }
    bestunit = NPY_FR_D;
    goto finish;
  }

  if ((*substr != 'T' && *substr != ' ') || sublen == 1) {
    goto parse_error;
  }
  comparison =
      compare_format(&format, &format_len, substr, 1, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
  ++substr;
  --sublen;

  /* PARSE THE HOURS */
  comparison =
      compare_format(&format, &format_len, "%H", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
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
      if (want_exc) {
        PyErr_Format(PyExc_ValueError,
                     "Hours out of range in datetime string \"%s\"", str);
      }
      goto error;
    }
  }

  /* Next character must be a ':' or the end of the string */
  if (sublen == 0) {
    if (!hour_was_2_digits) {
      goto parse_error;
    }
    if (format_len) {
      goto parse_error;
    }
    bestunit = NPY_FR_h;
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
    comparison =
        compare_format(&format, &format_len, ":", 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
  } else if (!isdigit(*substr)) {
    if (!hour_was_2_digits) {
      goto parse_error;
    }
    goto parse_timezone;
  }

  /* PARSE THE MINUTES */
  comparison =
      compare_format(&format, &format_len, "%M", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
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
      if (want_exc) {
        PyErr_Format(PyExc_ValueError,
                     "Minutes out of range in datetime string \"%s\"", str);
      }
      goto error;
    }
  } else if (!has_hms_sep) {
    goto parse_error;
  }

  if (sublen == 0) {
    bestunit = NPY_FR_m;
    if (format_len) {
      goto parse_error;
    }
    goto finish;
  }

  /* If we make it through this condition block, then the next
   * character is a digit. */
  if (has_hms_sep && *substr == ':') {
    comparison =
        compare_format(&format, &format_len, ":", 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
    ++substr;
    --sublen;
    /* Cannot have a trailing ':' */
    if (sublen == 0 || !isdigit(*substr)) {
      goto parse_error;
    }
  } else if (!has_hms_sep && isdigit(*substr)) {
  } else {
    goto parse_timezone;
  }

  /* PARSE THE SECONDS */
  comparison =
      compare_format(&format, &format_len, "%S", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
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
      if (want_exc) {
        PyErr_Format(PyExc_ValueError,
                     "Seconds out of range in datetime string \"%s\"", str);
      }
      goto error;
    }
  } else if (!has_hms_sep) {
    goto parse_error;
  }

  /* Next character may be a '.' indicating fractional seconds */
  if (sublen > 0 && *substr == '.') {
    ++substr;
    --sublen;
    comparison =
        compare_format(&format, &format_len, ".", 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
  } else {
    bestunit = NPY_FR_s;
    goto parse_timezone;
  }

  /* PARSE THE MICROSECONDS (0 to 6 digits) */
  comparison =
      compare_format(&format, &format_len, "%f", 2, format_requirement);
  if (comparison == COMPARISON_ERROR) {
    goto parse_error;
  } else if (comparison == COMPLETED_PARTIAL_MATCH) {
    goto finish;
  }
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
      bestunit = NPY_FR_us;
    } else {
      bestunit = NPY_FR_ms;
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
      bestunit = NPY_FR_ps;
    } else {
      bestunit = NPY_FR_ns;
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
    bestunit = NPY_FR_as;
  } else {
    bestunit = NPY_FR_fs;
  }

parse_timezone:
  /* trim any whitespace between time/timezone */
  while (sublen > 0 && isspace(*substr)) {
    ++substr;
    --sublen;
    comparison =
        compare_format(&format, &format_len, " ", 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
  }

  if (sublen == 0) {
    // Unlike NumPy, treating no time zone as naive
    if (format_len > 0) {
      goto parse_error;
    }
    goto finish;
  }

  /* UTC specifier */
  if (*substr == 'Z') {
    comparison =
        compare_format(&format, &format_len, "%z", 2, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
    /* "Z" should be equivalent to tz offset "+00:00" */
    if (out_local != NULL) {
      *out_local = 1;
    }

    if (out_tzoffset != NULL) {
      *out_tzoffset = 0;
    }

    if (sublen == 1) {
      if (format_len > 0) {
        goto parse_error;
      }
      goto finish;
    } else {
      ++substr;
      --sublen;
    }
  } else if (*substr == '-' || *substr == '+') {
    comparison =
        compare_format(&format, &format_len, "%z", 2, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
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
        if (want_exc) {
          PyErr_Format(PyExc_ValueError,
                       "Timezone hours offset out of range "
                       "in datetime string \"%s\"",
                       str);
        }
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
          if (want_exc) {
            PyErr_Format(PyExc_ValueError,
                         "Timezone minutes offset out of range "
                         "in datetime string \"%s\"",
                         str);
          }
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
    comparison =
        compare_format(&format, &format_len, " ", 1, format_requirement);
    if (comparison == COMPARISON_ERROR) {
      goto parse_error;
    } else if (comparison == COMPLETED_PARTIAL_MATCH) {
      goto finish;
    }
  }

  if ((sublen != 0) || (format_len != 0)) {
    goto parse_error;
  }

finish:
  if (out_bestunit != NULL) {
    *out_bestunit = bestunit;
  }
  return 0;

parse_error:
  if (want_exc) {
    PyErr_Format(PyExc_ValueError,
                 "Error parsing datetime string \"%s\" at position %d", str,
                 (int)(substr - str));
  }
  return -1;

error:
  return -1;
}

/* Converts the int64_t representation of a duration to ISO; mutates len */
static char *int64ToIsoDuration(int64_t value, size_t *len) {
  pandas_timedeltastruct tds;
  int ret_code;

  pandas_timedelta_to_timedeltastruct(value, NPY_FR_ns, &tds);

  // Max theoretical length of ISO Duration with 64 bit day
  // as the largest unit is 70 characters + 1 for a null terminator
  char *result = PyObject_Malloc(71);
  if (result == NULL) {
    PyErr_NoMemory();
    return NULL;
  }

  ret_code = make_iso_8601_timedelta(&tds, result, len);
  if (ret_code == -1) {
    PyErr_SetString(PyExc_ValueError,
                    "Could not convert timedelta value to string");
    PyObject_Free(result);
    return NULL;
  }

  return result;
}

/*
 * This function returns a pointer to the DateTimeMetaData
 * contained within the provided datetime dtype.
 *
 * Copied near-verbatim from numpy/core/src/multiarray/datetime.c
 */
static PyArray_DatetimeMetaData
get_datetime_metadata_from_dtype(PyArray_Descr *dtype) {
  return (((PyArray_DatetimeDTypeMetaData *)dtype->c_metadata)->meta);
}

static void pandas_datetime_destructor(PyObject *op) {
  void *ptr = PyCapsule_GetPointer(op, PandasDateTime_CAPSULE_NAME);
  PyMem_Free(ptr);
}

static int pandas_datetime_exec(PyObject *module) {
  PyDateTime_IMPORT;
  PandasDateTime_CAPI *capi = PyMem_Malloc(sizeof(PandasDateTime_CAPI));
  if (capi == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  capi->npy_datetimestruct_to_datetime = npy_datetimestruct_to_datetime;
  capi->scaleNanosecToUnit = scaleNanosecToUnit;
  capi->int64ToIso = int64ToIso;
  capi->NpyDateTimeToEpoch = NpyDateTimeToEpoch;
  capi->PyDateTimeToIso = PyDateTimeToIso;
  capi->PyDateTimeToEpoch = PyDateTimeToEpoch;
  capi->int64ToIsoDuration = int64ToIsoDuration;
  capi->pandas_datetime_to_datetimestruct = pandas_datetime_to_datetimestruct;
  capi->pandas_timedelta_to_timedeltastruct =
      pandas_timedelta_to_timedeltastruct;
  capi->convert_pydatetime_to_datetimestruct =
      convert_pydatetime_to_datetimestruct;
  capi->cmp_npy_datetimestruct = cmp_npy_datetimestruct;
  capi->get_datetime_metadata_from_dtype = get_datetime_metadata_from_dtype;
  capi->parse_iso_8601_datetime = parse_iso_8601_datetime;
  capi->get_datetime_iso_8601_strlen = get_datetime_iso_8601_strlen;
  capi->make_iso_8601_datetime = make_iso_8601_datetime;
  capi->make_iso_8601_timedelta = make_iso_8601_timedelta;

  PyObject *capsule = PyCapsule_New(capi, PandasDateTime_CAPSULE_NAME,
                                    pandas_datetime_destructor);
  if (capsule == NULL) {
    PyMem_Free(capi);
    return -1;
  }

  // Monkeypatch the top level pandas module to have an attribute for the
  // C-API. This is required because Python capsules do not support setting
  // this attribute on anything but the top level package. Ideally not
  // done when cpython gh-6898 gets implemented
  PyObject *pandas = PyImport_ImportModule("pandas");
  if (!pandas) {
    PyErr_SetString(PyExc_ImportError,
                    "pd_datetime.c could not import module pandas");
    Py_DECREF(capsule);
    // cpython doesn't PyMem_Free(capi) here - mistake or intentional?
    return -1;
  }

  if (PyModule_AddObject(pandas, "pandas_datetime_CAPI", capsule) < 0) {
    Py_DECREF(capsule);
    // cpython doesn't PyMem_Free(capi) here - mistake or intentional?
    return -1;
  }

  return 0;
}

static PyModuleDef_Slot pandas_datetime_slots[] = {
    {Py_mod_exec, pandas_datetime_exec}, {0, NULL}};

static struct PyModuleDef pandas_datetimemodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pandas._libs.pandas_datetime",

    .m_doc = "Internal module with datetime support for other extensions",
    .m_size = 0,
    .m_methods = NULL,
    .m_slots = pandas_datetime_slots};

PyMODINIT_FUNC PyInit_pandas_datetime(void) {
  return PyModuleDef_Init(&pandas_datetimemodule);
}
