#include "c_dates.h"
#include "c_freqs.h"
#include "c_convert.h"
#include "c_types.h"

#include <datetime.h>
#include <time.h>



#define INIT_ERR(errortype, errmsg) PyErr_SetString(errortype,errmsg);return -1


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



/* Returns absolute seconds from an hour, minute, and second
 */
#define secs_from_hms(hour, min, sec, multiplier) (\
  ((hour)*3600 + (min)*60 + (sec)) * (npy_int64)(multiplier)\
)

/* Returns the quarter */
#define month_to_quarter(month) (((month)-1)/3 + 1)
#define quarter_to_month(quarter) (((quarter)-1)*3 + 1)



/*
====================================================
== End of section borrowed from mx.DateTime       ==
====================================================
*/

void
init_metadata_from_unit(ts_metadata *meta, int unit) {
    meta->unit = (unit/1000)*1000;
    if (unit > FR_DAY) {
        meta->periods_per_day = highunits_per_day(unit);
        meta->secs_per_period = secs_per_highunits(unit, 1);
        meta->period_end_at = -1;
    }
    else {
        meta->periods_per_day = -1;
        meta->secs_per_period = -1;
        meta->period_end_at = unit - meta->unit;
    }
    meta->convert_to_start = 0;
    meta->timestep = 1;
}



// NPY_NO_EXPORT
ts_datetime
datetimestruct_to_tsdatetime(ts_metadata *meta, ts_datetimestruct *d)
{
    ts_datetime val;
    npy_int64 days=0;

    int unit = meta->unit;
    int period_end_at = meta->period_end_at;

    if ((unit > FR_MTH) || (unit == FR_UND)) {
        days = days_from_ymd(d->year, d->month, d->day);
        if (unit >= FR_HR){
            days -= HIGHFREQ_ORIG;
        };
    };
    if (unit == FR_ANN) {
        val = d->year;
    }
    else if (unit == FR_QTR) {
        npy_int64 quarter = ((d->month - 1)/3) + 1;
        if (period_end_at > 12) {
            // quarterly frequency with year determined by ending period
            val = d->year*4 + quarter;
        }
        else {
            /* quarterly frequency with year determined by ending period
                           or has December year end*/
            val = (d->year - 1)*4 + quarter;
        };
    }
    else if (unit == FR_MTH) {
        val = (d->year-1)*12 + d->month;
    }
    else if (unit == FR_WK) {
        ts_datetime end_week_day, adj_day;
        end_week_day = (7 - period_end_at) % 7;
        adj_day = days + ((7 - end_week_day) - days % 7) % 7;
        val = adj_day / 7;
    }
    else if (unit == FR_BUS) {
        npy_int64 weeks = days/7;
        val = days - weeks * 2;
        /*
        int dotw = day_of_week(days);
        if (dotw > 4){
            // Invalid business day
            val = 0;
        }
        else {
            npy_int64 x = days -2;
            val = 2 + (x/7)*5 + x%7;
        }
        */
    }
    else if ((unit == FR_DAY) || (unit==FR_UND)){
        val = days;
    }
    else if (unit == FR_HR) {
        val = days * 24 + d->hour;
    }
    else if (unit == FR_MIN){
        val = days * 1440 + d->hour * 60 + d->min;
    }
    else if (unit == FR_SEC){
        val = days *  (npy_int64)(86400)
            + secs_from_hms(d->hour, d->min, d->sec, 1);
    }
    else {
        /* Shouldn't get here */
        PyErr_SetString(PyExc_ValueError, "invalid internal frequency");
        val = -1;
    }
    return val;
}



// NPY_NO_EXPORT
ts_datetime
PyDatetime_to_tsdatetime(ts_metadata *meta, PyObject *datetime)
{
    ts_datetimestruct dinfo;
    ts_datetime val;

    if (!PyDateTime_Check(datetime) && !PyDate_Check(datetime)){
        PyObject *err_msg, *_type;
        _type = PyObject_Type(datetime);
        err_msg = PyString_FromString("Expected a datetime.date(time) object, received: ");
        PyString_ConcatAndDel(&err_msg, PyObject_Str(_type));
        PyErr_SetString(PyExc_TypeError, PyString_AsString(err_msg));
        Py_DECREF(_type);
        Py_DECREF(err_msg);
        val = -1;
    }
    else {
        dinfo.year = (npy_int64)PyDateTime_GET_YEAR(datetime);
        dinfo.month = PyDateTime_GET_MONTH(datetime);
        dinfo.day = (int)PyDateTime_GET_DAY(datetime);
        dinfo.hour = (int)PyDateTime_DATE_GET_HOUR(datetime);
        dinfo.min = (int)PyDateTime_DATE_GET_MINUTE(datetime);
        dinfo.sec = (int)PyDateTime_DATE_GET_SECOND(datetime);
        //
        val = datetimestruct_to_tsdatetime(meta, &dinfo);
    }
    return val;
}



// NPY_NO_EXPORT
ts_datetime
PyString_to_tsdatetime(ts_metadata *meta, PyObject *string)
{
    ts_datetime val;
    PyObject *pydatetime;
    PyObject *string_arg = PyTuple_New(1);

    PyTuple_SET_ITEM(string_arg, 0, string);
    Py_INCREF(string);
    if (meta->unit > FR_DAY) {
        pydatetime = PyEval_CallObject(DateTimeFromString, string_arg);
    } else {
        pydatetime = PyEval_CallObject(DateFromString, string_arg);
    };
    Py_DECREF(string_arg);
    val = PyDatetime_to_tsdatetime(meta, pydatetime);
    Py_DECREF(pydatetime);
    return val;
}



PyObject *
PyDateTime_FromFreqAndValue(ts_datetime val, ts_metadata *meta)
{
    PyObject *pydatetime;
    npy_int64 absdate, abstime=0;
    int hour=0, minute=0, second=0;
    ts_datetimestruct dinfo;

    int unit = meta->unit;

    conversion_function todays = get_converter_to_days(unit, 0);
    meta->convert_to_start = 0;

    absdate = todays(val, meta);
    if (unit > FR_DAY)
        abstime = _secs_from_midnight(val, unit);
    set_datetimestruct_from_days_and_secs(&dinfo, absdate, abstime);

    switch(unit) {
        case FR_HR:
            hour = dinfo.hour;
            break;
        case FR_MIN:
            hour = dinfo.hour;
            minute = dinfo.min;
            break;
        case FR_SEC:
            hour = dinfo.hour;
            minute = dinfo.min;
            second = (int)dinfo.sec;
            break;
    }
    pydatetime = PyDateTime_FromDateAndTime(dinfo.year, dinfo.month, dinfo.day,
                                            hour, minute, second, 0);
    return pydatetime;
}





/************************************************************
** Date type definition
************************************************************/

static void
DatetimeObject_dealloc(DatetimeObject* self) {
    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
DatetimeObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    DatetimeObject *self;
    self = (DatetimeObject*)type->tp_alloc(type, 0);

    // Initialize the metadata
    ts_metadata meta = {FR_UND, 1, 0, -1, -1, 0};

    if (self != NULL) {
        // initialize attributes that need initializing in here
        self->obval = -1;
        self->obmeta = meta;
    }
    return (PyObject *)self;
}


/* for use in C code */
static DatetimeObject *
DatetimeObject_New(void) {
    PyObject *dummy=NULL;
    return (DatetimeObject*)DatetimeObject_new(&DatetimeObject_Type, dummy, dummy);
}



static int
DatetimeObject_init(DatetimeObject *self, PyObject *args, PyObject *kwds) {

    PyObject *freq=NULL, *value=NULL, *datetime=NULL, *string=NULL;
    char *INSUFFICIENT_MSG = "insufficient parameters to initialize Date";

    int def_info=INT_ERR_CODE;

    int year=def_info, month=def_info, day=def_info, quarter=def_info,
        hour=def_info, minute=def_info, second=def_info;
    int unit;

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

    if (PyObject_HasAttrString(freq, "unit")) {
        PyObject *freq_attr = PyObject_GetAttrString(freq, "unit");
        unit = PyInt_AS_LONG(freq_attr);
        Py_DECREF(freq_attr);
    } else {
        if((unit = check_freq(freq)) == INT_ERR_CODE) return -1;
    }
    init_metadata_from_unit(&self->obmeta, unit);
    unit = self->obmeta.unit;

    // The input value is a date string...
    if ((value && PyString_Check(value)) || string) {
        if (!string)
            string = value;
        self->obval = PyString_to_tsdatetime(&(self->obmeta), string);
    }
    // The input is a datetime.date(time) object:
    else if ((value && (PyDateTime_Check(value) || PyDate_Check(value))) ||
             datetime) {
        if (!datetime)
            datetime = value;
        self->obval = PyDatetime_to_tsdatetime(&(self->obmeta), datetime);
    }
    // The input is a value
    else if (value) {
        self->obval = PyInt_AsLong(value);
    }
    // The input is a list of arguments
    else {
        // We always must have a year defined
        if (year == def_info) {
            INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
        };
        // For weekly to daily (undefined) units, we need month and day
        if (unit == FR_BUS || unit == FR_DAY || unit == FR_WK || unit == FR_UND) {
            if (month == def_info || day == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            };
                // if FR_BUS, check for week day
        }
        // Need a month for monthly
        else if (unit == FR_MTH) {
            if (month == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
        }
        else if (unit == FR_QTR) {
            if (quarter == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
            month = (quarter-1) * 3 + 1;
        }
        else if (unit == FR_SEC) {
            if (month == def_info || day == def_info || second == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
            if (hour == def_info) {
                hour = second/3600;
                minute = (second % 3600)/60;
                second = second % 60;
            }
            else if (minute == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
        }
        else if (unit == FR_MIN) {
            if (month == def_info || day == def_info || minute == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
            if (hour == def_info) {
                hour = minute/60;
                minute = minute % 60;
            }
        }
        else if (unit == FR_HR) {
            if (month == def_info || day == def_info || hour == def_info) {
                INIT_ERR(PyExc_ValueError, INSUFFICIENT_MSG);
            }
        };
        /*Fill the structure from the input, and get the corresponding datetime*/
        ts_datetimestruct dinfo;
        dinfo.year = year;
        dinfo.month = month;
        dinfo.day = day;
        dinfo.hour = hour;
        dinfo.min = minute;
        dinfo.sec = second;
        self->obval = datetimestruct_to_tsdatetime(&(self->obmeta), &dinfo);
    }
    return 0;
}


static PyMemberDef DatetimeObject_members[] = {
    {"value", T_INT, offsetof(DatetimeObject, obval), 0,
     "integer representation of the Date"},
    {NULL}  /* Sentinel */
};



typedef struct {
    PyObject_HEAD;
    int unit; /* frequency of date */
    int freq;
    long years;
    long months;
    long days;
    long seconds;
//    PyObject* cached_vals;
} TimeDeltaObject;

/* Forward declarations */
static PyTypeObject TimeDeltaType;
#define TimeDelta_Check(op) PyObject_TypeCheck(op, &TimeDeltaType)

static void
TimeDeltaObject_dealloc(TimeDeltaObject* self) {
//    Py_XDECREF(self->cached_vals);
    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
TimeDeltaObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    TimeDeltaObject *self;
    self = (TimeDeltaObject*)type->tp_alloc(type, 0);
    if (self != NULL) {
        // initialize attributes that need initializing in here
        self->unit = FR_UND;
        self->freq = self->unit;
        self->months = 0;
        self->days = 0;
        self->seconds = 0;
    }
    return (PyObject *)self;
}


static TimeDeltaObject *
TimeDeltaObject_New(void) {
    PyObject *dummy=NULL;
    return (TimeDeltaObject*)TimeDeltaObject_new(&TimeDeltaType, dummy, dummy);
}


#define INIT_ERR(errortype, errmsg) PyErr_SetString(errortype,errmsg);return -1

#define get_timedelta_months(o)  (((TimeDeltaObject *)(o))->months)
#define get_timedelta_days(o)    (((TimeDeltaObject *)(o))->days)
#define get_timedelta_seconds(o) (((TimeDeltaObject *)(o))->seconds)



static int
TimeDeltaObject_init(TimeDeltaObject *self, PyObject *args, PyObject *kwds) {

    PyObject *unit=NULL, *freq=NULL, *value=NULL, *delta=NULL;
    PyObject *py_years=NULL, *py_months=NULL, *py_days=NULL, *py_quarters=NULL;
    PyObject *py_hours=NULL, *py_minutes=NULL, *py_seconds=NULL;
    npy_int64 months=0, days=0, seconds=0;

    int fr_group;
    int free_dt=0;

    static char *kwlist[] = {"unit", "value",
                             "years", "months", "days", "quarters",
                             "hours", "minutes", "seconds",
                             "timedelta", "freq", NULL};

    // Check the parameters
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|OOOOOOOOO:__new__", kwlist,
                                      &unit, &value,
                                      &py_years, &py_months, &py_days, &py_quarters,
                                      &py_hours, &py_minutes, &py_seconds,
                                      &delta, &freq)) {
        return -1;
    }

    // If `freq` is defined, then overwrite `unit`
    if (freq){
        value = unit;
        unit = freq;
        freq = NULL;
    }

    // Process `unit`
    if (PyObject_HasAttrString(unit, "unit")) {
        PyObject *unit_attr = PyObject_GetAttrString(unit, "unit");
        self->unit = PyInt_AS_LONG(unit_attr);
        Py_DECREF(unit_attr);
    }
    else if (PyObject_HasAttrString(unit, "freq")) {
        PyObject *unit_attr = PyObject_GetAttrString(unit, "freq");
        self->unit = PyInt_AS_LONG(unit_attr);
        Py_DECREF(unit_attr);
    }
    else {
        if ((self->unit = check_freq(unit)) == INT_ERR_CODE)
            return -1;
    };
    self->freq = self->unit;
    fr_group = get_base_unit(self->unit);

    if (value && PyDelta_Check(value)) {
        if (!delta) delta = value;
        value = NULL;
    }; // datetime = (datetime||value), value = NULL

    if (value) {
        if (fr_group == FR_ANN)
            months = PyInt_AsLong(value) * 12;
        else if (fr_group == FR_QTR)
            months = PyInt_AsLong(value) * 3;
        else if (fr_group == FR_MTH)
            months = PyInt_AsLong(value);
        else if (fr_group == FR_WK)
            days = PyInt_AsLong(value) * 7;
        else if ((fr_group == FR_BUS) || (fr_group == FR_DAY) || (fr_group == FR_UND))
            days = PyInt_AsLong(value);
        else if (fr_group == FR_HR)
            seconds = PyInt_AsLong(value) * 3600;
        else if (fr_group == FR_MIN)
            seconds = PyInt_AsLong(value) * 60;
        else if (fr_group == FR_SEC)
            seconds = PyInt_AsLong(value);
        else
            days = PyInt_AsLong(value);
        Py_DECREF(value);
    };
    if (delta) {
        if (PyDelta_Check(delta)){
            self->days = ((PyDateTime_Delta *)(delta))->days;
            self->seconds = ((PyDateTime_Delta *)(delta))->seconds + \
                            ((PyDateTime_Delta *)(delta))->microseconds/1000000;
//            free_dt = 1;
        }
        else {
            PyObject *err_msg, *_type;
            _type = PyObject_Type(delta);
            err_msg = PyString_FromString("Expected timedelta object, received: ");
            PyString_ConcatAndDel(&err_msg, PyObject_Str(_type));
            PyErr_SetString(PyExc_TypeError, PyString_AsString(err_msg));
            Py_DECREF(_type);
            Py_DECREF(err_msg);
            return -1;
        }
    } else {
        if (py_years){
            months += PyInt_AsLong(py_years)*12;
//            Py_DECREF(py_years);
        };
        if (py_quarters){
            months += PyInt_AsLong(py_quarters)*3;
//            Py_DECREF(py_quarters);
        };
        if (py_months){
            months += PyInt_AsLong(py_months);
//            Py_DECREF(py_months);
        };
        if (py_days){
            days += PyInt_AsLong(py_days);
//            Py_DECREF(py_days);
        };
        if (py_hours){
            seconds += PyInt_AsLong(py_hours) * 3600;
//            Py_DECREF(py_hours);
        };
        if (py_minutes){
            seconds += PyInt_AsLong(py_minutes) * 60;
//            Py_DECREF(py_minutes);
        };
        if (py_seconds){
            seconds += PyInt_AsLong(py_seconds);
//            Py_DECREF(py_seconds);
        };
        normalize_days_secs(&days, &seconds);

        self->months = months;
        self->days = days;
        self->seconds = seconds;
        };

    if (free_dt) { Py_DECREF(delta); }
    return 0;
}



static PyMemberDef TimeDeltaObject_members[] = {
    {"unit", T_INT, offsetof(TimeDeltaObject, unit), 0,
     "frequency"},
    {"freq", T_INT, offsetof(TimeDeltaObject, freq), 0,
     "frequency"},
    {"months", T_INT, offsetof(TimeDeltaObject, months), 0,
     "months"},
    {"days", T_INT, offsetof(TimeDeltaObject, days), 0,
     "days"},
    {"seconds", T_INT, offsetof(TimeDeltaObject, seconds), 0,
     "seconds"},
    {NULL}  /* Sentinel */
};



static char DatetimeObject_toordinal_doc[] =
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
DatetimeObject_toordinal(DatetimeObject *self)
{
    conversion_function todays = get_converter_to_days(self->obmeta.unit, 0);
    return PyInt_FromLong(todays(self->obval, &(self->obmeta)));
};



static char DatetimeObject_asfreq_doc[] =
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


static PyObject*
DatetimeObject_convert(DatetimeObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *freq=NULL;
    char *relation_raw=NULL;
    char *relation_uc;
    conversion_function tomediator, frommediator;
    char relation;
    int invalid_relation=0;
    int tofreq;

    static char *kwlist[] = {"freq", "relation", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|s", kwlist,
                                      &freq, &relation_raw)) return NULL;

    if(relation_raw) {
        if (strlen(relation_raw) > 0) {
            if ((relation_uc = str_uppercase(relation_raw)) == NULL)
                return PyErr_NoMemory();
            // 'BEFORE' and 'AFTER' values for this parameter are deprecated
            if (strcmp(relation_uc, "END") == 0 ||
                strcmp(relation_uc, "E") == 0 ||
                strcmp(relation_uc, "START") == 0 ||
                strcmp(relation_uc, "S") == 0 ||
                strcmp(relation_uc, "BEFORE") == 0 ||
                strcmp(relation_uc, "B") == 0 ||
                strcmp(relation_uc, "AFTER") == 0 ||
                strcmp(relation_uc, "A") == 0) {
                 if (relation_uc[0] == 'E' || relation_uc[0] == 'A')
                     relation = 'E';
                 else
                     relation = 'S';
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

    if (relation == 'S')
        self->obmeta.convert_to_start = 1;
    else
        self->obmeta.convert_to_start = 0;

    if ((tofreq = check_freq(freq)) == INT_ERR_CODE)
        return NULL;

    DatetimeObject *result = DatetimeObject_New();

    init_metadata_from_unit(&result->obmeta, tofreq);
    result->obmeta.convert_to_start = self->obmeta.convert_to_start;

    int fromfreq = self->obmeta.unit;
    if (tofreq == fromfreq) {
        result->obval = self->obval;
        return (PyObject*)result;
    }

    /* Correction for business days */
    if ((tofreq == FR_BUS) && (fromfreq < FR_DAY))
        result->obmeta.convert_to_start = 1;

    tomediator = convert_to_mediator(fromfreq, tofreq, 1);
    frommediator = convert_from_mediator(fromfreq, tofreq, 1);
    result->obval = frommediator(tomediator(self->obval, &(self->obmeta)), &(result->obmeta));
    return (PyObject *)result;
}



static char DatetimeObject_strfmt_doc[] =
"Deprecated alias for strftime method";

static char DatetimeObject_strftime_doc[] =
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
DatetimeObject_strftime(DatetimeObject *self, PyObject *args)
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
    ymdstruct ymd;
    hmsstruct hms;
    npy_int64 absdate;
    double abstime;
    int i, result_len;
    PyObject *py_result;

    ts_metadata meta = self->obmeta;

    if (!PyArg_ParseTuple(args, "s:strftime(fmt)", &orig_fmt_str))
        return NULL;

    conversion_function convert_to_days = get_converter_to_days(meta.unit, 0);
    meta.convert_to_start = 0;

    absdate = convert_to_days(self->obval, &meta);
    ymd = days_to_ymdstruct(absdate, GREGORIAN_CALENDAR);
    abstime = _secs_from_midnight(self->obval, meta.unit);
    hms = seconds_to_hmsstruct(abstime);

    /* Populate standard C date struct with info from our date_info struct */
    c_date.tm_sec = hms.sec;
    c_date.tm_min = hms.min;
    c_date.tm_hour = hms.hour;
    c_date.tm_mday = ymd.day;
    c_date.tm_mon = ymd.month - 1;
    c_date.tm_year = ymd.year - 1900;
    // c_date.tm_wday = (tempDate.day_of_week + 1) % 7;
    // c_date.tm_yday = tempDate.day_of_year - 1;
    c_date.tm_isdst = -1;

    result_len = strlen(orig_fmt_str) + 50;
    if ((result = PyArray_malloc(result_len * sizeof(char))) == NULL)
        return PyErr_NoMemory();

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
            if (extra_fmts_found_one)
                free(tmp_str);

            if (fmt_str == NULL)
                return NULL;

            extra_fmts_found[i] = 1;
            extra_fmts_found_one = 1;
        }
    }

    strftime(result, result_len, fmt_str, &c_date);
    if (extra_fmts_found_one)
        free(fmt_str);

    // replace any place holders with the appropriate value
    for(i=0; i < num_extra_fmts; i++) {
        if (extra_fmts_found[i]) {
            char *tmp_str = result;
            char *extra_str;

            if (strcmp(extra_fmts[i][0], "%q") == 0 ||
                strcmp(extra_fmts[i][0], "%f") == 0 ||
                strcmp(extra_fmts[i][0], "%F") == 0) {

                ts_metadata meta_qtr;
                if (self->obmeta.unit == FR_QTR)
                    meta_qtr = self->obmeta;
                else
                    init_metadata_from_unit(&meta_qtr, FR_QTR);

                int year, quarter, qvals, year_len;
                meta_qtr.convert_to_start = 0;
//                conversion_info qtr_info;
//                metadata_to_conversioninfo(&meta_qtr, &qtr_info);
//                set_conversion_info(qtr_freq, 'E', &qtr_info);

                qvals = (*get_converter_from_days(meta_qtr.unit, 0))(absdate, &meta_qtr);
                quarter = qvals % 4;
                quarter = (quarter == 0 ? 4 : quarter);
                year = (qvals - quarter)/4 + 1;

                if(strcmp(extra_fmts[i][0], "%q") == 0) {
                    if ((extra_str = PyArray_malloc(2 * sizeof(char))) == NULL) {
                        free(tmp_str);
                        return PyErr_NoMemory();
                    }
                    sprintf(extra_str, "%i", quarter);
                } else {
                    if (meta_qtr.period_end_at > 12)
                        year -= 1;

                    if (strcmp(extra_fmts[i][0], "%f") == 0) {
                        year_len = 2;
                        year = year % 100;
                    }
                    else
                        year_len = 4;

                    if ((extra_str = PyArray_malloc((year_len+1) * sizeof(char))) == NULL) {
                        free(tmp_str);
                        return PyErr_NoMemory();
                    }

                    if (year_len == 2 && year < 10)
                        sprintf(extra_str, "0%i", year);
                    else
                        sprintf(extra_str, "%i", year);
                }
            }
            else {
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
DatetimeObject___str__(DatetimeObject* self)
{
    int unit = self->obmeta.unit;
    PyObject *string_arg, *retval;

    string_arg = NULL;
    if (unit == FR_UND) {
        retval = PyString_FromFormat("%ld", (long)(self->obval));
        return retval;
        }
    else if (unit == FR_ANN) { string_arg = Py_BuildValue("(s)", "%Y"); }
    else if (unit == FR_QTR) { string_arg = Py_BuildValue("(s)", "%FQ%q"); }
    else if (unit == FR_MTH) { string_arg = Py_BuildValue("(s)", "%b-%Y"); }
    else if (unit == FR_DAY ||
             unit == FR_BUS ||
             unit == FR_WK) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y"); }
    else if (unit == FR_HR) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:00"); }
    else if (unit == FR_MIN) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M"); }
    else if (unit == FR_SEC) { string_arg = Py_BuildValue("(s)", "%d-%b-%Y %H:%M:%S"); }

    if (string_arg == NULL) { return NULL; }
    retval = DatetimeObject_strftime(self, string_arg);
    Py_DECREF(string_arg);
    return retval;
}
static PyObject *
timedelta___str__(TimeDeltaObject* self)
{
    PyObject  *py_str;
    long months=get_timedelta_months(self);
    long days=get_timedelta_days(self);
    long seconds=get_timedelta_seconds(self);

    py_str = PyString_FromFormat("%ldm, %ldd, %lds",
                                  months, days, seconds);
    return py_str;
}



static PyObject *
DatetimeObject_freqstr(DatetimeObject *self, void *closure) {
    PyObject *key = PyInt_FromLong(self->obmeta.unit);
    PyObject *freq_aliases = PyDict_GetItem(freq_dict, key);
    PyObject *main_alias = PyTuple_GET_ITEM(freq_aliases, 0);
    Py_DECREF(key);
    Py_INCREF(main_alias);
    return main_alias;
}
static PyObject *
timedelta_freqstr(TimeDeltaObject *self, void *closure) {
    PyObject *key = PyInt_FromLong(self->unit);
    PyObject *freq_aliases = PyDict_GetItem(freq_dict, key);
    PyObject *main_alias = PyTuple_GET_ITEM(freq_aliases, 0);
    Py_DECREF(key);
    Py_INCREF(main_alias);
    return main_alias;
}



static PyObject *
DatetimeObject___repr__(DatetimeObject* self)
{
    PyObject *py_str_rep, *py_freqstr, *py_repr;
    char *str_rep, *freqstr, *repr;
    int repr_len;

    py_str_rep = DatetimeObject___str__(self);
    if (py_str_rep == NULL) return NULL;
    str_rep = PyString_AsString(py_str_rep);

    py_freqstr = DatetimeObject_freqstr(self, NULL);
    freqstr = PyString_AsString(py_freqstr);

    repr_len = strlen(str_rep) + strlen(freqstr) + 6;
    MEM_CHECK((repr = PyArray_malloc((repr_len + 1) * sizeof(char))));

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

static PyObject *
timedelta___repr__(TimeDeltaObject* self)
{
    PyObject *py_freqstr, *py_repr;
    char *freqstr;
    py_freqstr = timedelta_freqstr(self, NULL);
    freqstr = PyString_AsString(py_freqstr);

    if (get_timedelta_seconds(self) != 0){
        py_repr = PyString_FromFormat("%s(%ld, %ld, %ld)[%s]",
                                      self->ob_type->tp_name,
                                      self->months,
                                      self->days,
                                      self->seconds,
                                      freqstr);
        }
    else if (self->days != 0){
        py_repr = PyString_FromFormat("%s(%ld, %ld)[%s]",
                                      self->ob_type->tp_name,
                                      self->months,
                                      self->days,
                                      freqstr);
        }
    else {
        py_repr = PyString_FromFormat("%s(%ld)[%s]",
                                      self->ob_type->tp_name,
                                      self->months,
                                      freqstr);
    }
    Py_DECREF(py_freqstr);
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


DatetimeObject *
DatetimeObject_FromFreqAndValue(int freq, npy_int64 value) {
    DatetimeObject *result = DatetimeObject_New();
    init_metadata_from_unit(&result->obmeta, freq);
    result->obval = value;
    return result;
}


static TimeDeltaObject *
timedelta_fromMDS(int unit,
                  npy_int64 months, npy_int64 days, npy_int64 seconds)
{
    TimeDeltaObject *result = TimeDeltaObject_New();
    normalize_days_secs(&days, &seconds);
    result->unit = unit;
    result->months = months;
    result->days = days;
    result->seconds = seconds;
    return result;
}



static PyObject *
DatetimeObject_date_plus_int(PyObject *date, PyObject *pyint) {
    DatetimeObject *dateobj = (DatetimeObject*)date;
    if (!PyInt_Check(pyint) && !PyObject_HasAttrString(pyint, "__int__")) {
        return Py_NotImplemented;
    }
    DatetimeObject *result = DatetimeObject_New();
    result->obmeta = dateobj->obmeta;
    result->obval = dateobj->obval + PyInt_AsLong(pyint);
    return (PyObject*)result;
}


static PyObject *
date_plus_timedelta(PyObject *datearg, PyObject *deltaarg){
    DatetimeObject *date = (DatetimeObject*)datearg;
    TimeDeltaObject *delta = (TimeDeltaObject*)deltaarg;
    DatetimeObject *result;

    ts_datetimestruct dtinfo;
    npy_int64 seconds, days, months, years=0, absdate, abstime=0;
    // Get the info from the delta
    seconds = delta->seconds;
    days = delta->days;
    months = delta->months;

    if (date->obmeta.unit > FR_DAY) {
        abstime = _secs_from_midnight(date->obval, date->obmeta.unit);
        // Update the days and secs from the seconds since midnight
        seconds += abstime;
    }
    // Make sure the seconds are between -86400 and +86400
    normalize_days_secs(&days, &seconds);
    set_datetimestruct_from_secs(&dtinfo, seconds);

    // Convert the current date to days and set the datetimestruct
    conversion_function todays = get_converter_to_days(date->obmeta.unit, 0);
    date->obmeta.convert_to_start = 1;
    absdate = todays(date->obval, &(date->obmeta));
    set_datetimestruct_from_days(&dtinfo, absdate + days);

    // Update the datetime info with the months and normalize the months
    years = dtinfo.year;
    months += dtinfo.month;
    normalize_years_months(&years, &months);

    // Update the structure
    dtinfo.year = years;
    dtinfo.month = months;

    // Convert to datetime
    result = DatetimeObject_New();
    result->obmeta = date->obmeta;
    result->obval = datetimestruct_to_tsdatetime(&(date->obmeta), &dtinfo);
    return (PyObject*)result;

};

static PyObject *
date_plus_delta(PyObject *datearg, PyObject *deltaarg){
    DatetimeObject *date = (DatetimeObject*)datearg;
    PyDateTime_Delta *delta = (PyDateTime_Delta*)deltaarg;
    PyObject *timedelta=NULL;
    timedelta = (PyObject*)timedelta_fromMDS(date->obmeta.unit,
                                             0, delta->days, delta->seconds);
    PyObject *result=NULL;
    result = (PyObject *)date_plus_timedelta(datearg, timedelta);
    Py_DECREF(timedelta);
    return result;
}



static PyObject *
DatetimeObject___add__(PyObject *left, PyObject *right)
{
    if (DatetimeObject_Check(left) && DatetimeObject_Check(right)) {
        PyErr_SetString(PyExc_TypeError, "Cannot add Date to Date");
        return NULL;
    }
    PyObject *result = Py_NotImplemented;

    if (DatetimeObject_Check(left)) {
        if (TimeDelta_Check(right)){
            result = date_plus_timedelta(left, right);
        } else if (PyDelta_Check(right)) {
            result = date_plus_delta(left,right);
        } else {
            result = DatetimeObject_date_plus_int(left, right);
        };
    } else {
        if (TimeDelta_Check(left)){
            result = date_plus_timedelta(right, left);
        } else {
            result = DatetimeObject_date_plus_int(right, left);
        };
    };
    if (result == Py_NotImplemented)
        Py_INCREF(result);
    return result;
};


static PyObject *
DatetimeObject___subtract__(PyObject *left, PyObject *right)
{
    PyObject *result = Py_NotImplemented;

    DatetimeObject *dleft;
    if (!DatetimeObject_Check(left)) {
        PyErr_SetString(PyExc_ValueError, "Cannot subtract a Date from a non-Date object.");
        return NULL;
    }
    dleft = (DatetimeObject*)left;
    if (DatetimeObject_Check(right)) {
        DatetimeObject *dright = (DatetimeObject*)right;
        if (dleft->obmeta.unit != dright->obmeta.unit) {
            PyErr_SetString(PyExc_ValueError, "Cannot subtract Date objects with different frequencies.");
            result = NULL;
        }
        result = PyInt_FromLong(dleft->obval - dright->obval);
    }
    else if (PyDate_Check(right) || PyDateTime_Check(right)) {
        DatetimeObject *dright=DatetimeObject_New();
        PyObject *args = PyTuple_New(0);
        PyObject *kw = PyDict_New();
        PyObject *py_unit = PyInt_FromLong(dleft->obmeta.unit);
        PyDict_SetItemString(kw, "freq", py_unit);
        PyDict_SetItemString(kw, "datetime", right);
        Py_DECREF(py_unit);
        ///
        DatetimeObject_init(dright, args, kw);
        Py_DECREF(args);
        Py_DECREF(kw);
        result = PyInt_FromLong(dleft->obval - dright->obval);
        Py_DECREF(dright);
    }
    else if (TimeDelta_Check(right)){
        PyObject *minus_right = PyNumber_Negative(right);
        if (minus_right){
            result = date_plus_timedelta(left, minus_right);
            Py_DECREF(minus_right);
        } else {
            result = NULL;
        }
    }
    else {
        DatetimeObject *dtresult = DatetimeObject_New();
        dtresult->obmeta = dleft->obmeta;
        dtresult->obval = dleft->obval - PyLong_AsLong(right);
        result = (PyObject *)dtresult;
    }

    if (result == Py_NotImplemented)
        Py_INCREF(result);
    return result;
}

static PyObject *
timedelta_negative(TimeDeltaObject *self){
    TimeDeltaObject *result = TimeDeltaObject_New();
    result->unit = self->unit;
    result->months = -get_timedelta_months(self);
    result->days = -get_timedelta_days(self);
    result->seconds = -get_timedelta_seconds(self);
    return (PyObject *)result;
}


static PyObject *
timedelta_plus_timedelta(PyObject *tdaobj, PyObject *tdbobj) {
    TimeDeltaObject *tda = (TimeDeltaObject*)tdaobj;
    TimeDeltaObject *tdb = (TimeDeltaObject*)tdbobj;
    npy_int64 months = get_timedelta_months(tda) + get_timedelta_months(tdb);
    npy_int64 days = get_timedelta_days(tda) + get_timedelta_days(tdb);
    npy_int64 seconds = get_timedelta_seconds(tda) + get_timedelta_seconds(tdb);
    normalize_days_secs(&days, &seconds);
    //
    return (PyObject*)timedelta_fromMDS(tda->unit, months, days, seconds);

};


static PyObject *
timedelta_plus_int(PyObject *timedelta, PyObject *pyint) {
    TimeDeltaObject *deltaobj = (TimeDeltaObject*)timedelta;
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
        sprintf(err_str, "Cannot add TimeDelta and %s", type_str);
        Py_DECREF(obj_type);
        Py_DECREF(type_repr);
        PyErr_SetString(PyExc_TypeError, err_str);
        free(err_str);
        return NULL;
    }
    int freq_group = get_base_unit(deltaobj->unit);
    long years=0, months=0, days=0, seconds=0;
    switch(freq_group){
        case FR_ANN:
            years = PyInt_AsLong(pyint);
            break;
        case FR_QTR:
            months = PyInt_AsLong(pyint);
            months *= 3;
            break;
        case FR_MTH:
            months = PyInt_AsLong(pyint);
            break;
        case FR_WK:
            days = PyInt_AsLong(pyint);
            days *= 7;
            break;
        case FR_BUS:
            days = PyInt_AsLong(pyint);
            break;
        case FR_DAY:
            days = PyInt_AsLong(pyint);
            break;
        case FR_HR:
            seconds = PyInt_AsLong(pyint);
            seconds *= 3600;
            break;
        case FR_MIN:
            seconds = PyInt_AsLong(pyint);
            seconds *= 60;
            break;
        case FR_SEC:
            seconds = PyInt_AsLong(pyint);
            break;
        default:
            days = PyInt_AsLong(pyint);
            break;
    }
    return (PyObject*)timedelta_fromMDS(deltaobj->unit,
                                        deltaobj->months+months+12*years,
                                        deltaobj->days+days,
                                        deltaobj->seconds+seconds);
};

static PyObject *
timedelta_plus_delta(PyObject *left, PyObject *right){
    TimeDeltaObject *oleft = (TimeDeltaObject*)left;
    PyDateTime_Delta *oright = (PyDateTime_Delta*)right;
    PyObject *result=NULL;
    npy_int64 days = oleft->days + oright->days;
    npy_int64 seconds = oleft->seconds + oright->seconds;
    normalize_days_secs(&days, &seconds);
    result = (PyObject*)timedelta_fromMDS(oleft->unit,
                                          oleft->months, days, seconds);
    return result;
}


static PyObject *
timedelta_add(PyObject *left, PyObject *right)
{
    PyObject *result = Py_NotImplemented;

    if (TimeDelta_Check(left)) {
        if (TimeDelta_Check(right))
            result = timedelta_plus_timedelta(left, right);
        else if (PyDelta_Check(right))
            result = timedelta_plus_delta(left, right);
        else if (PyInt_Check(right) || PyLong_Check(right))
            result = timedelta_plus_int(left, right);
        }
    else if (PyDelta_Check(left))
        result = timedelta_plus_delta(right, left);
    else if (PyInt_Check(left) || PyLong_Check(left))
        result = timedelta_plus_int(right,left);

    if (result == Py_NotImplemented)
        Py_INCREF(result);
    return result;
};


static PyObject *
timedelta_subtract(PyObject *left, PyObject *right)
{
    PyObject *result = Py_NotImplemented;

    if (!TimeDelta_Check(left))
        PyErr_SetString(PyExc_ValueError, "Cannot subtract a TimeDelta from a non-TimeDelta object.");

    if (TimeDelta_Check(right)) {
        PyObject *minus_right = PyNumber_Negative(right);
        result = timedelta_plus_timedelta(left, minus_right);
    }
    else {
        PyObject *minus_right = PyNumber_Negative(right);
        if (minus_right) {
            if (PyDelta_Check(right))
                result = timedelta_plus_delta(left, minus_right);
            else if (TimeDelta_Check(right))
                result = timedelta_plus_timedelta(left, minus_right);
            else
                result = timedelta_plus_int(left, minus_right);
            Py_DECREF(minus_right);
        }
        else {
            result = NULL;
        };
    };
    if (result == Py_NotImplemented)
        Py_INCREF(result);
    return result;
}


static PyObject *
timedelta_times_int(PyObject *delta, PyObject *py_int)
{
    TimeDeltaObject *deltaobj = (TimeDeltaObject*)delta;
    TimeDeltaObject *result = TimeDeltaObject_New();
    long factor = PyInt_AsLong(py_int);
    result->unit = deltaobj->unit;
    result->months = get_timedelta_months(deltaobj) * factor;
    result->days = get_timedelta_days(deltaobj) * factor;
    result->seconds = get_timedelta_seconds(deltaobj) * factor;
    return (PyObject *)result;
};


#define NUM_CHECK(o) (PyInt_Check(o) || PyLong_Check(o) || PyFloat_Check(o))


static PyObject *
timedelta_multiply(PyObject *left, PyObject *right) {
    PyObject *result = Py_NotImplemented;

    if (TimeDelta_Check(left)) {
        if (NUM_CHECK(right))
            result = timedelta_times_int(left, right);
    }
    else if (NUM_CHECK(left)) {
        result = timedelta_times_int(right, left);
    };
    if (result == Py_NotImplemented)
        Py_INCREF(result);
    return result;
}




static int
DatetimeObject___compare__(DatetimeObject * obj1, DatetimeObject * obj2)
{
    if (obj1->obmeta.unit != obj2->obmeta.unit) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot compare Date objects with different frequencies.");
        return -1;
    }
    if (obj1->obval < obj2->obval) return -1;
    if (obj1->obval > obj2->obval) return 1;
    if (obj1->obval == obj2->obval) return 0;
    return -1;
}

static int
TimeDeltaObject___compare__(TimeDeltaObject * obj1, TimeDeltaObject * obj2)
{
//    if (obj1->unit != obj2->unit) {
//        PyErr_SetString(PyExc_ValueError,
//                        "Cannot compare Date objects with different frequencies.");
//        return -1;
//    }
//    if (obj1->value < obj2->value) return -1;
//    if (obj1->value > obj2->value) return 1;
//    if (obj1->value == obj2->value) return 0;
//    return -1;
    return -1;
}



static long
DatetimeObject___hash__(DatetimeObject *self)
{
    /* within a given frequency, hash values are guaranteed to be unique
       for different dates. For different frequencies, we make a reasonable
       effort to ensure hash values will be unique, but it is not guaranteed */
    if (self->obmeta.unit == FR_BUS) {
        return self->obval + 10000000;
    } else if (self->obmeta.unit == FR_WK) {
        return self->obval + 100000000;
    } else {
        return self->obval;
    };
}
static long
TimeDeltaObject___hash__(TimeDeltaObject *self)
{
//    register int freq_group = get_base_unit(self->freq);
//    /* within a given frequency, hash values are guaranteed to be unique
//       for different dates. For different frequencies, we make a reasonable
//       effort to ensure hash values will be unique, but it is not guaranteed */
//    if (freq_group == FR_BUS) {
//        return self->value + 10000000;
//    } else if (freq_group == FR_WK) {
//        return self->value + 100000000;
//    } else {
//        return self->value;
//    };
    return -1;
}


static PyObject *
DatetimeObject_toint(DatetimeObject *self) {
    return PyInt_FromLong(self->obval);
}

static PyObject *
DatetimeObject_tofloat(DatetimeObject *self) {
    return PyFloat_FromDouble((double)(self->obval));
}

static PyObject *
DatetimeObject_tolong(DatetimeObject *self) {
    return PyLong_FromLong(self->obval);
}



/***************************************************
           ====== Date Properties ======
****************************************************/

static PyObject *
DatetimeObject_unit(DatetimeObject *self) {
    int unit = self->obmeta.unit;
    int period_end_at = self->obmeta.period_end_at;
    if (period_end_at > 0)
        return PyInt_FromLong(unit + period_end_at);
    return PyInt_FromLong(unit);
}
static PyObject *
DatetimeObject_timestep(DatetimeObject *self) {
    return PyInt_FromLong(self->obmeta.timestep);
}


// helper function for date property funcs
static int
DatetimeObject_set_datestruct(DatetimeObject *self, ts_datetimestruct *dinfo) {
    conversion_function todays = get_converter_to_days(self->obmeta.unit, 0);
    self->obmeta.convert_to_start = 0;
    npy_int64 days = todays(self->obval, &(self->obmeta));
    set_datetimestruct_from_days(dinfo, days);
    return 0;
}

// helper function for date property funcs
static int
DatetimeObject_set_datetimestruct(DatetimeObject *self, ts_datetimestruct *dinfo) {
    conversion_function todays = get_converter_to_days(self->obmeta.unit, 0);
    self->obmeta.convert_to_start = 0;
    npy_int64 absdate = todays(self->obval, &(self->obmeta));
    npy_int64 abstime = _secs_from_midnight(self->obval, self->obmeta.unit);
    set_datetimestruct_from_days_and_secs(dinfo, absdate, abstime);
    return 0;
}


static PyObject *
DatetimeObject_year(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    DatetimeObject_set_datestruct(self, &dinfo);
    return PyInt_FromLong(dinfo.year);
}
PyObject *
_loop_get_year(npy_int64 value, int unit,
               conversion_function todays, ts_metadata *meta,
               ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(dinfo->year);
}


static PyObject *
DatetimeObject_qyear(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    DatetimeObject_set_datestruct(self, &dinfo);
    int year = dinfo.year;

    if (self->obmeta.unit == FR_QTR) {
        int end_month = self->obmeta.period_end_at;
        if (end_month == 0)
            end_month = 12;
        else if (end_month > 12) {
            end_month -= 12;
            year -= 1;
        }
        if (dinfo.month > end_month)
            year += 1;
    }
    return PyInt_FromLong(year);
}
PyObject *
_loop_get_qyear(npy_int64 value, int unit,
                conversion_function todays, ts_metadata *meta,
                ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(dinfo->year);
}
PyObject *
_loop_get_qyear_from_qtr(npy_int64 value, int unit,
                         conversion_function todays, ts_metadata *meta,
                         ts_datetimestruct *dinfo)
{
    int ending_month = meta->period_end_at;
    if (ending_month == 0)
        ending_month = 12;
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    if (dinfo->month > ending_month)
        return PyInt_FromLong(dinfo->year+1);
    return PyInt_FromLong(dinfo->year);
}


static PyObject *
DatetimeObject_quarter(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    DatetimeObject_set_datestruct(self, &dinfo);
    int month = dinfo.month;
    int end_month = ending_month(&(self->obmeta));

    if (self->obmeta.unit == FR_QTR)
        month = dinfo.month - end_month;
        if (month <= 0)
            month += 12;
    return PyInt_FromLong(month_to_quarter(month));
}
PyObject *
_loop_get_quarter(npy_int64 value, int unit,
                  conversion_function todays, ts_metadata *meta,
                  ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(month_to_quarter(dinfo->month));
}
PyObject *
_loop_get_quarter_from_qtr(npy_int64 value, int unit,
                           conversion_function todays, ts_metadata *meta,
                           ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    int month = dinfo->month - ending_month(meta);
    if (month <= 0)
        month += 12;
    return PyInt_FromLong(month_to_quarter(month));
}


static PyObject *
DatetimeObject_month(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datestruct(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.month);
}
PyObject *
_loop_get_month(npy_int64 value, int unit,
                conversion_function todays, ts_metadata *meta,
                ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(dinfo->month);
}

static PyObject *
DatetimeObject_day(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datestruct(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.day);
}
PyObject *
_loop_get_day(npy_int64 value, int unit,
              conversion_function todays, ts_metadata *meta,
              ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(dinfo->day);
}


static PyObject *
DatetimeObject_day_of_week(DatetimeObject *self, void *closure) {
    PyObject *daily_obj;
    daily_obj = DatetimeObject_toordinal(self);
    npy_int64 absdate = PyInt_AsLong(daily_obj);
    Py_DECREF(daily_obj);
    return PyInt_FromLong(day_of_week(absdate));
}
static PyObject *
DatetimeObject_weekday(DatetimeObject *self, void *closure) {
    return DatetimeObject_day_of_week(self, closure);
}
PyObject *
_loop_get_day_of_week(npy_int64 value, int unit,
                      conversion_function todays, ts_metadata *meta,
                      ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    return PyInt_FromLong(day_of_week(absdate));
}


static PyObject *
DatetimeObject_day_of_year(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datestruct(self, &dinfo) == -1)
        return NULL;
    return PyInt_FromLong(dinfo.day_of_year);
}
PyObject *
_loop_get_day_of_year(npy_int64 value, int unit,
                      conversion_function todays, ts_metadata *meta,
                      ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(dinfo->day_of_year);
}


static PyObject *
DatetimeObject_week(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datestruct(self, &dinfo) == -1)
        return NULL;
    return PyInt_FromLong(isoweek_from_datetimestruct(&dinfo));
}
PyObject *
_loop_get_week(npy_int64 value, int unit,
               conversion_function todays, ts_metadata *meta,
               ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    set_datetimestruct_from_days(dinfo, absdate);
    return PyInt_FromLong(isoweek_from_datetimestruct(dinfo));
}


static PyObject *
DatetimeObject_hour(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datetimestruct(self, &dinfo) == -1)
        return NULL;
    return PyInt_FromLong(dinfo.hour);
}
PyObject *
_loop_get_hour(npy_int64 value, int unit,
               conversion_function todays, ts_metadata *meta,
               ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    npy_int64 abstime = _secs_from_highfreq(value, meta) % 86400;
    set_datetimestruct_from_days_and_secs(dinfo, absdate, abstime);
    return PyInt_FromLong(dinfo->hour);
}


static PyObject *
DatetimeObject_minute(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datetimestruct(self, &dinfo) == -1)
        return NULL;
    return PyInt_FromLong(dinfo.min);
}
PyObject *
_loop_get_minute(npy_int64 value, int unit,
                 conversion_function todays, ts_metadata *meta,
                 ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    npy_int64 abstime = _secs_from_highfreq(value, meta) % 86400;
    set_datetimestruct_from_days_and_secs(dinfo, absdate, abstime);
    return PyInt_FromLong(dinfo->min);
}


static PyObject *
DatetimeObject_second(DatetimeObject *self, void *closure) {
    ts_datetimestruct dinfo;
    if(DatetimeObject_set_datetimestruct(self, &dinfo) == -1) return NULL;
    return PyInt_FromLong(dinfo.sec);
}
PyObject *
_loop_get_second(npy_int64 value, int unit,
                 conversion_function todays, ts_metadata *meta,
                 ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    npy_int64 abstime = _secs_from_highfreq(value, meta) % 86400;
    set_datetimestruct_from_days_and_secs(dinfo, absdate, abstime);
    return PyInt_FromLong(dinfo->sec);
}


//static PyObject *
//DatetimeObject_ordinal(DatetimeObject *self, void *closure){
//    return PyInt_FromLong((long)DatetimeObject_toordinal(self));
//}
PyObject *
_loop_get_ordinal(npy_int64 value, int unit,
                  conversion_function todays, ts_metadata *meta,
                  ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    return PyInt_FromLong(absdate);
}


static PyObject *
DatetimeObject_datetime(DatetimeObject *self, void *closure) {
    PyObject *datetime;
    int hour=0, minute=0, second=0;
    int unit;
    ts_datetimestruct dinfo;

    if(DatetimeObject_set_datetimestruct(self, &dinfo) == -1) return NULL;
    unit = self->obmeta.unit;

    switch(unit) {
        case FR_HR:
            hour = dinfo.hour;
            break;
        case FR_MIN:
            hour = dinfo.hour;
            minute = dinfo.min;
            break;
        case FR_SEC:
            hour = dinfo.hour;
            minute = dinfo.min;
            second = (int)dinfo.sec;
            break;
    }
    datetime = PyDateTime_FromDateAndTime(
                dinfo.year, dinfo.month, dinfo.day, hour, minute, second, 0);
    return datetime;
}
PyObject *
_loop_get_datetime(npy_int64 value, int unit,
                   conversion_function todays, ts_metadata *meta,
                   ts_datetimestruct *dinfo)
{
    npy_int64 absdate = (*todays)(value, meta);
    npy_int64 abstime = 0;
    if (unit > FR_DAY)
        abstime = _secs_from_midnight(value, unit);
    set_datetimestruct_from_days_and_secs(dinfo, absdate, abstime);
    return PyDateTime_FromDateAndTime(dinfo->year, dinfo->month, dinfo->day,
                                      dinfo->hour, dinfo->min, dinfo->sec, 0);
}



static PyObject *
TimeDeltaObject_timedelta(TimeDeltaObject *self, void *closure) {
    PyObject *timedelta;
    npy_int64 days = 30.4375 * get_timedelta_months(self);
    days += get_timedelta_days(self);
    npy_int64 seconds=self->seconds;
    timedelta = PyDelta_FromDSU(days, seconds, 0);
    return timedelta;
}



static int
DatetimeObject_ReadOnlyErr(DatetimeObject *self, PyObject *value, void *closure) {
   PyErr_SetString(PyExc_AttributeError, "Cannot set read-only property");
   return -1;
}
static int
TimeDeltaObject_ReadOnlyErr(TimeDeltaObject *self, PyObject *value, void *closure) {
   PyErr_SetString(PyExc_AttributeError, "Cannot set read-only property");
   return -1;
}

static PyGetSetDef DatetimeObject_getseters[] = {
    {"unit", (getter)DatetimeObject_unit, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the unit.", NULL},
    {"freq", (getter)DatetimeObject_unit, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the unit.", NULL},
    {"timestep", (getter)DatetimeObject_timestep, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the timestep.", NULL},
    {"year", (getter)DatetimeObject_year, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the year.", NULL},
    {"qyear", (getter)DatetimeObject_qyear, (setter)DatetimeObject_ReadOnlyErr,
            "For quarterly frequency dates, returns the year corresponding to the\n"
            "year end (start) month. When using QTR or QTR-E based quarterly\n"
            "frequencies, this is the fiscal year in a financial context.\n\n"
            "For non-quarterly dates, this simply returns the year of the date.",
            NULL},
    {"quarter", (getter)DatetimeObject_quarter, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the quarter.", NULL},
    {"month", (getter)DatetimeObject_month, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the month.", NULL},
    {"week", (getter)DatetimeObject_week, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the week.", NULL},
    {"day", (getter)DatetimeObject_day, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the day of month.", NULL},
    {"weekday", (getter)DatetimeObject_weekday, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the day of week.", NULL},
    // deprecated alias for weekday property
    {"day_of_week", (getter)DatetimeObject_weekday, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the day of week.", NULL},
    {"day_of_year", (getter)DatetimeObject_day_of_year, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the day of year.", NULL},
    {"second", (getter)DatetimeObject_second, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the second.", NULL},
    {"minute", (getter)DatetimeObject_minute, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the minute.", NULL},
    {"hour", (getter)DatetimeObject_hour, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the hour.", NULL},
    {"freqstr", (getter)DatetimeObject_freqstr, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the string representation of frequency.", NULL},
    {"datetime", (getter)DatetimeObject_datetime, (setter)DatetimeObject_ReadOnlyErr,
            "Returns the Date object converted to standard python datetime object",
            NULL},
    {NULL}  /* Sentinel */
};
static PyGetSetDef TimeDeltaObject_getseters[] = {
//    {"year", (getter)DeltaObject_year, (setter)DeltaObject_ReadOnlyErr,
//            "Returns the year.", NULL},
//    {"quarter", (getter)DeltaObject_quarter, (setter)DeltaObject_ReadOnlyErr,
//            "Returns the quarter.", NULL},
//    {"month", (getter)DeltaObject_month, (setter)DeltaObject_ReadOnlyErr,
//            "Returns the month.", NULL},
//    {"week", (getter)DeltaObject_week, (setter)DeltaObject_ReadOnlyErr,
//            "Returns the week.", NULL},
//    // deprecated alias for weekday property
//    {"minute", (getter)DeltaObject_minute, (setter)DeltaObject_ReadOnlyErr,
//            "Returns the minute.", NULL},
//    {"hour", (getter)DeltaObject_hour, (setter)DeltaObject_ReadOnlyErr,
//            "Returns the hour.", NULL},
    {"freqstr", (getter)timedelta_freqstr, (setter)TimeDeltaObject_ReadOnlyErr,
            "Returns the string representation of frequency.", NULL},
    {"timedelta", (getter)TimeDeltaObject_timedelta, (setter)TimeDeltaObject_ReadOnlyErr,
            "Returns the Delta object converted to standard python timedelta object",
            NULL},
    {NULL}  /* Sentinel */
};




static PyNumberMethods DatetimeObject_as_number = {
    (binaryfunc)DatetimeObject___add__,      /* nb_add */
    (binaryfunc)DatetimeObject___subtract__, /* nb_subtract */
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
    (unaryfunc)DatetimeObject_toint,       /* nb_int */
    (unaryfunc)DatetimeObject_tolong,      /* nb_long */
    (unaryfunc)DatetimeObject_tofloat,     /* nb_float */
    (unaryfunc)0,                        /* nb_oct */
    (unaryfunc)0,                        /* nb_hex */
};
static PyNumberMethods TimeDeltaObject_as_number = {
    (binaryfunc)timedelta_add,       /* nb_add */
    (binaryfunc)timedelta_subtract,  /* nb_subtract */
    (binaryfunc)timedelta_multiply,  /* nb_multiply */
    0,                               /* nb_divide */
    0,                               /* nb_remainder */
    0,                               /* nb_divmod */
    0,                               /* nb_power */
    (unaryfunc)timedelta_negative,   /* nb_negative */
    0,                               /* nb_positive */
    0,                               /* nb_absolute */
    0,                               /* nb_nonzero */
    0,                               /* nb_invert */
    0,                               /* nb_lshift */
    0,                               /* nb_rshift */
    0,                               /* nb_and */
    0,                               /* nb_xor */
    0,                               /* nb_or */
    0,                               /* nb_coerce */
    0,                               /* nb_int */
    0,                               /* nb_long */
    0,                               /* nb_float */
    (unaryfunc)0,                    /* nb_oct */
    (unaryfunc)0,                    /* nb_hex */
};



static PyMethodDef DatetimeObject_methods[] = {
    {"toordinal", (PyCFunction)DatetimeObject_toordinal, METH_NOARGS,
     DatetimeObject_toordinal_doc},
    {"strftime", (PyCFunction)DatetimeObject_strftime, METH_VARARGS,
     DatetimeObject_strftime_doc},
    // deprecated alias for strftime
    {"strfmt", (PyCFunction)DatetimeObject_strftime, METH_VARARGS,
     DatetimeObject_strfmt_doc},
    {"asfreq", (PyCFunction)DatetimeObject_convert, METH_VARARGS | METH_KEYWORDS,
     DatetimeObject_asfreq_doc},
    {"convert", (PyCFunction)DatetimeObject_convert, METH_VARARGS | METH_KEYWORDS,
     DatetimeObject_asfreq_doc},
    {NULL}  /* Sentinel */
};
//static PyMethodDef TimeDeltaObject_methods[] = {
////    {"toordinal", (PyCFunction)DeltaObject_toordinal, METH_NOARGS,
////     DeltaObject_toordinal_doc},
//    {NULL}  /* Sentinel */
//};



PyTypeObject DatetimeObject_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                               /* ob_size */
    "timeseries.Date",               /* tp_name */
    sizeof(DatetimeObject),              /* tp_basicsize */
    0,                               /* tp_itemsize */
    (destructor)DatetimeObject_dealloc,  /* tp_dealloc */
    0,                               /* tp_print */
    0,                               /* tp_getattr */
    0,                               /* tp_setattr */
    (cmpfunc)DatetimeObject___compare__, /* tp_compare */
    (reprfunc)DatetimeObject___repr__,   /* tp_repr */
    &DatetimeObject_as_number,           /* tp_as_number */
    0,                               /* tp_as_sequence */
    0,                               /* tp_as_mapping */
    (hashfunc)DatetimeObject___hash__,   /* tp_hash */
    0,                               /* tp_call*/
    (reprfunc)DatetimeObject___str__,    /* tp_str */
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
    DatetimeObject_methods,              /* tp_methods */
    DatetimeObject_members,              /* tp_members */
    DatetimeObject_getseters,            /* tp_getset */
    0,                               /* tp_base */
    0,                               /* tp_dict */
    0,                               /* tp_descr_get */
    0,                               /* tp_descr_set */
    0,                               /* tp_dictoffset */
    (initproc)DatetimeObject_init,       /* tp_init */
    0,                               /* tp_alloc */
    DatetimeObject_new,                  /* tp_new */
};

static PyTypeObject TimeDeltaType = {
    PyObject_HEAD_INIT(NULL)
    0,                                    /* ob_size */
    "timeseries.TimeDelta",               /* tp_name */
    sizeof(TimeDeltaObject),              /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)TimeDeltaObject_dealloc,  /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    (cmpfunc)TimeDeltaObject___compare__, /* tp_compare */
    (reprfunc)timedelta___repr__,   /* tp_repr */
    &TimeDeltaObject_as_number,           /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    (hashfunc)TimeDeltaObject___hash__,   /* tp_hash */
    0,                                    /* tp_call*/
    (reprfunc)timedelta___str__,    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |                  /* tp_flags */
    Py_TPFLAGS_CHECKTYPES |
    Py_TPFLAGS_BASETYPE,
    "Defines a Delta object, as the combination of a date and a frequency.\n",  /* tp_doc */
    0,                                /* tp_traverse */
    0,                                /* tp_clear */
    0,                                /* tp_richcompare */
    0,                                /* tp_weaklistoffset */
    0,                                /* tp_iter */
    0,                                /* tp_iternext */
    0,              /* tp_methods */
//    DeltaObject_methods,              /* tp_methods */
    TimeDeltaObject_members,              /* tp_members */
    TimeDeltaObject_getseters,            /* tp_getset */
    0,                                /* tp_base */
    0,                                /* tp_dict */
    0,                                /* tp_descr_get */
    0,                                /* tp_descr_set */
    0,                                /* tp_dictoffset */
    (initproc)TimeDeltaObject_init,       /* tp_init */
    0,                                    /* tp_alloc */
    TimeDeltaObject_new,                  /* tp_new */
};



///////////////////////////////////////////////////////////////////////
#define DICT_SETINT_STRKEY(dict, key, val) {\
        PyObject *pyval = PyInt_FromLong(val); \
        PyDict_SetItemString(dict, key, pyval); \
        Py_DECREF(pyval); }


PyObject *
c_dates_now(PyObject *self, PyObject *args) {

    PyObject *freq; //, *init_args, *init_kwargs;
    struct tm *timeinfo;
    ts_datetimestruct dinfo;
    int freq_val;
    DatetimeObject *date;


#ifdef WIN32
    __time64_t rawtime;
#else
    time_t rawtime;
#endif

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

    //
    dinfo.year = timeinfo->tm_year + 1900;
    dinfo.month = timeinfo->tm_mon + 1;
    dinfo.day = timeinfo->tm_mday;
    dinfo.hour = timeinfo->tm_hour;
    dinfo.min = timeinfo->tm_min;
    dinfo.sec = timeinfo->tm_sec;

    date = DatetimeObject_New();
    init_metadata_from_unit(&date->obmeta, freq_val);
    date->obval = datetimestruct_to_tsdatetime(&(date->obmeta), &dinfo);
    return (PyObject*)date;
}



PyObject *
DateArray_asfreq(PyObject *self, PyObject *args)
{
    PyArrayObject *fromDates, *toDates;
    PyArrayIterObject *iterFrom, *iterTo;
    PyObject *fromDateObj, *toDateObj;
    char *relation;
    char relation_from, relation_to;
    conversion_function converterfrom, converterto;
//    conversion_info infofrom, infoto;
    int fromfreq, tofreq;
    npy_int64 fromDate, toDate, tmpDate;

    ts_metadata input_meta, output_meta;

    if (!PyArg_ParseTuple(args,
                "Oiis:asfreq(fromDates, fromfreq, tofreq, relation)",
                &fromDates, &fromfreq, &tofreq, &relation)) return NULL;

    relation_from = relation[0];
    if ((tofreq == FR_BUS) && (fromfreq < FR_DAY))
        relation_to = 'S';
    else
        relation_to = relation_from;

    toDates = (PyArrayObject *)PyArray_Copy(fromDates);

    iterFrom = (PyArrayIterObject *)PyArray_IterNew((PyObject *)fromDates);
    if (iterFrom == NULL)
        return NULL;

    iterTo = (PyArrayIterObject *)PyArray_IterNew((PyObject *)toDates);
    if (iterTo == NULL)
        return NULL;

    init_metadata_from_unit(&input_meta, fromfreq);
    if (relation_from == 'S')
        input_meta.convert_to_start = 1;
    fromfreq = input_meta.unit;
    converterfrom = convert_to_mediator(fromfreq, tofreq, 0);

    init_metadata_from_unit(&output_meta, tofreq);
    if (relation_to == 'S')
        output_meta.convert_to_start = 1;
    tofreq = output_meta.unit;
    converterto = convert_from_mediator(fromfreq, tofreq, 0);


    while (iterFrom->index < iterFrom->size) {

        fromDateObj = PyArray_GETITEM(fromDates, iterFrom->dataptr);
        fromDate = PyInt_AsLong(fromDateObj);
        tmpDate = converterfrom(fromDate, &input_meta);
        toDate = converterto(tmpDate, &output_meta);
        // ERR_CHECK(toDate = asfreq_main(fromDate, relation[0], &af_info));
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
int __skip_periods_year(int freq) {
    int freq_group = get_base_unit(freq);
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
int __skip_periods_quarter(int freq) {
    int freq_group = get_base_unit(freq);
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
int __skip_periods_month(int freq) {
    int freq_group = get_base_unit(freq);
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
int __skip_periods_day(int freq) {
    switch(freq)
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
int __skip_periods_week(int freq) {
    switch(freq)
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
int __skip_periods_hour(int freq) {
    switch(freq)
    {
        case FR_MIN:
            return 60;
        case FR_SEC:
            return 60*60;
        default:
            return 1;
    }
}
int __skip_periods_minute(int freq) {
    switch(freq)
    {
        case FR_SEC:
            return 60;
        default:
            return 1;
    }
}


PyObject *
DateArray_getdateinfo(PyObject *self, PyObject *args)
{
    int freq, is_full, skip_periods, counter=1, val_changed=0;
    char *infochar;

    ts_metadata meta;

    PyObject *prev_val=NULL;
    PyArrayObject *input, *output;
    PyArrayIterObject *iterin, *iterout;

//    PyObject* (*getdateparam)(DatetimeObject*, void*) = NULL;
    PyObject* (*getdateparam)(npy_int64, int,
                              conversion_function, ts_metadata*,
                              ts_datetimestruct*) = NULL;

    if (!PyArg_ParseTuple(args,
                          "Oisi:getDateInfo(array, freq, info, is_full)",
                          &input, &freq, &infochar, &is_full))
        return NULL;
    if (*infochar == 'P') {
        output = (PyArrayObject *)PyArray_SimpleNew(input->nd, input->dimensions, NPY_OBJECT);
    }
    else {
        output = (PyArrayObject *)PyArray_Copy(input);
    };

    conversion_function todays = get_converter_to_days(freq, 1);
    init_metadata_from_unit(&meta, freq);
    meta.convert_to_start = 0;
    ts_datetimestruct dinfo;


    iterin = (PyArrayIterObject *)PyArray_IterNew((PyObject *)input);
    iterout = (PyArrayIterObject *)PyArray_IterNew((PyObject *)output);

    switch(*infochar)
    {
        case 'Y': //year
            getdateparam = &_loop_get_year;
            skip_periods = __skip_periods_year(freq);
            break;
        case 'F': //"fiscal" year
            if (get_base_unit(freq) == FR_QTR)
                getdateparam = &_loop_get_qyear_from_qtr;
            else
                getdateparam = &_loop_get_qyear;
            skip_periods = __skip_periods_year(freq);
            break;
        case 'Q': //quarter
            if (get_base_unit(freq) == FR_QTR)
                getdateparam = &_loop_get_quarter_from_qtr;
            else
                getdateparam = &_loop_get_quarter;
            skip_periods = __skip_periods_quarter(freq);
            break;
        case 'M': //month
            getdateparam = &_loop_get_month;
            skip_periods = __skip_periods_month(freq);
            break;
        case 'D': //day
            getdateparam = &_loop_get_day;
            skip_periods = __skip_periods_day(freq);
            break;
        case 'R': //day of year
            getdateparam = &_loop_get_day_of_year;
            skip_periods = __skip_periods_day(freq);
            break;
        case 'W': //day of week
            getdateparam = &_loop_get_day_of_week;
            skip_periods = __skip_periods_day(freq);
            break;
        case 'I': //week of year
            getdateparam = &_loop_get_week;
            skip_periods = __skip_periods_week(freq);
            break;
        case 'H': //hour
            getdateparam = &_loop_get_hour;
            skip_periods = __skip_periods_hour(freq);
            break;
        case 'T': //minute
            getdateparam = &_loop_get_minute;
            skip_periods = __skip_periods_minute(freq);
            break;
        case 'S': //second
            getdateparam = &_loop_get_second;
            skip_periods = 1;
            break;
        case 'O': //toordinal
            getdateparam = &_loop_get_ordinal;
            skip_periods = __skip_periods_day(freq);
            break;
        default:
            return NULL;
    }

    {
        PyObject *val, *result=NULL;
        while (iterin->index < iterin->size) {

            if ((val_changed == 0) ||
                (is_full == 0) ||
                (prev_val == NULL) ||
                (counter >= skip_periods)) {

                   val = PyArray_GETITEM(input, iterin->dataptr);
                   result = getdateparam(PyInt_AsLong(val), freq,
                                         todays, &meta, &dinfo);

                   if ((prev_val != NULL) &&
                       (PyLong_AsLong(prev_val) != PyLong_AsLong(result))) {
                       val_changed = 1;
                       counter = 0;
                   }
                   Py_DECREF(val);
//                   Py_DECREF(curr_date);
                   if (prev_val != NULL) {
                       Py_DECREF(prev_val);
                   }
                   prev_val = result;
            }

            PyArray_SETITEM(output, iterout->dataptr, result);

            PyArray_ITER_NEXT(iterin);
            PyArray_ITER_NEXT(iterout);

            counter += 1;
        }
    }
    if (prev_val != NULL) {
        Py_DECREF(prev_val);
    }
    Py_DECREF(iterin);
    Py_DECREF(iterout);
    return (PyObject *) output;
}


PyObject *
DateArray_getdatetime(PyObject *self, PyObject *args)
{
    int freq;

    PyArrayObject *input, *output;
    PyArrayIterObject *iterin, *iterout;

    if (!PyArg_ParseTuple(args,
                          "Oi:getdatetime(array, freq)", &input, &freq))
        return NULL;
    output = (PyArrayObject *)PyArray_SimpleNew(input->nd, input->dimensions, NPY_OBJECT);

    conversion_function todays = get_converter_to_days(freq, 1);
    ts_metadata meta = ((DatetimeObject *)self)->obmeta;
    meta.convert_to_start = 1;
    ts_datetimestruct dinfo;

    iterin = (PyArrayIterObject *)PyArray_IterNew((PyObject *)input);
    iterout = (PyArrayIterObject *)PyArray_IterNew((PyObject *)output);

    {
    PyObject *val, *result=NULL;
    while (iterin->index < iterin->size) {
        val = PyArray_GETITEM(input, iterin->dataptr);
        result = _loop_get_datetime(PyInt_AsLong(val), freq,
                                    todays, &meta, &dinfo);
        PyArray_SETITEM(output, iterout->dataptr, result);
        PyArray_ITER_NEXT(iterin);
        PyArray_ITER_NEXT(iterout);
        }
    }
    Py_DECREF(iterin);
    Py_DECREF(iterout);
    return (PyObject *) output;
}




void import_c_dates(PyObject *m)
{

    if (PyType_Ready(&DatetimeObject_Type) < 0) return;
    if (PyType_Ready(&TimeDeltaType) < 0) return;

    DateCalc_Error = PyErr_NewException("c_dates.DateCalc_Error", NULL, NULL);
    DateCalc_RangeError = PyErr_NewException("c_dates.DateCalc_RangeError",
                                             NULL, NULL);

    import_array();
    PyDateTime_IMPORT;

    Py_INCREF(&DatetimeObject_Type);
    PyModule_AddObject(m, "Date", (PyObject *)(&DatetimeObject_Type));
    Py_INCREF(&TimeDeltaType);
    PyModule_AddObject(m, "TimeDelta", (PyObject *)(&TimeDeltaType));

    PyModule_AddObject(m, "DateCalc_Error", DateCalc_Error);
    PyModule_AddObject(m, "DateCalc_RangeError", DateCalc_RangeError);

}
