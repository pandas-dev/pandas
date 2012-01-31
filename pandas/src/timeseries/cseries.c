#include "c_lib.h"
#include "c_dates.h"
#include "c_tseries.h"

static PyMethodDef cseries_methods[] = {

    {"MA_mov_sum", (PyCFunction)MaskedArray_mov_sum,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_median", (PyCFunction)MaskedArray_mov_median,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_min", (PyCFunction)MaskedArray_mov_min,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_max", (PyCFunction)MaskedArray_mov_max,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_average_expw", (PyCFunction)MaskedArray_mov_average_expw,
     METH_VARARGS | METH_KEYWORDS, ""},

    {"TS_convert", (PyCFunction)TimeSeries_convert,
     METH_VARARGS, ""},

    {"DA_asfreq", (PyCFunction)DateArray_asfreq,
     METH_VARARGS, ""},
    {"DA_getDateInfo", (PyCFunction)DateArray_getDateInfo,
     METH_VARARGS, ""},


    {"now", (PyCFunction)c_dates_now,
     METH_VARARGS,
        "now(freq)\n"
        "\n"
        "Returns the current date/time, at the given frequency\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "freq : {freq_spec}\n"
        "   Frequency to convert the Date to. Accepts any valid frequency\n"
        "   specification (string or integer)\n"},

    {"check_freq", (PyCFunction)c_dates_check_freq,
     METH_VARARGS,
        "Translates a user specified frequency into the corresponding frequency constant"},

    {"check_freq_str", (PyCFunction)c_dates_check_freq_str,
     METH_VARARGS,
        "Translates a user specified frequency into standard string representation"},

    {"get_freq_group", (PyCFunction)c_dates_get_freq_group,
     METH_VARARGS,
        "translate user specified frequency into frequency group constant"},


    {"set_callback_DateFromString", (PyCFunction)set_callback_DateFromString,
     METH_VARARGS, ""},
    {"set_callback_DateTimeFromString", (PyCFunction)set_callback_DateTimeFromString,
     METH_VARARGS, ""},

    {NULL, NULL}
};

PyMODINIT_FUNC
initcseries(void)
{
    PyObject *m;

    m = Py_InitModule("cseries", cseries_methods);
    if (m == NULL)
      return;

    import_c_lib(m);
    import_c_dates(m);
    import_c_tseries(m);

}
