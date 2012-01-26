#include "c_freqs.h"



//static PyObject *freq_dict, *freq_dict_rev, *freq_constants;

#define DICT_SETINT_STRKEY(dict, key, val) {\
        PyObject *pyval = PyInt_FromLong(val); \
        PyDict_SetItemString(dict, key, pyval); \
        Py_DECREF(pyval); }

#define ADD_FREQ_CONSTANT(const_name, val) \
    DICT_SETINT_STRKEY(freq_constants, const_name, val)

#define INIT_FREQ(const_name, key, aliases) \
    {PyObject *pykey = PyInt_FromLong(key); \
     PyDict_SetItem(freq_dict, pykey, aliases); \
     PyDict_SetItemString(freq_constants, const_name, pykey); \
     Py_DECREF(pykey); \
     Py_DECREF(aliases); }




static int init_freq_group(int num_items, int num_roots, int base_const,
                           char item_abbrevs[][2][10], char group_prefixes[][15],
                           char item_const_names[][15]) {
    int i;

    for (i = 0; i < num_items; i++) {

        PyObject *aliases;
        int j, size, k;

        if (i == 0) { k = 3; } else { k = 2; }

        size = num_roots * k;

        aliases = PyTuple_New(size);

        for (j = 0; j < num_roots; j++) {
            PyObject *alias_v1, *alias_v2;
            char *root, *alt;

            if ((root = PyArray_malloc((30) * sizeof(char))) == NULL) return INT_ERR_CODE;
            if ((alt = PyArray_malloc((30) * sizeof(char))) == NULL) return INT_ERR_CODE;

            strcpy(root, group_prefixes[j]);
            strcpy(alt, group_prefixes[j]);

            if (i == 0) {
                PyObject *alias = PyString_FromString(root);
                PyTuple_SET_ITEM(aliases, j*k + 2, alias);
            }

            strcat(root, "-");
            strcat(root, item_abbrevs[i][0]);
            strcat(alt, "-");
            strcat(alt, item_abbrevs[i][1]);

            alias_v1 = PyString_FromString(root);
            alias_v2 = PyString_FromString(alt);

            free(root);
            free(alt);

            PyTuple_SET_ITEM(aliases, j*k, alias_v1);
            PyTuple_SET_ITEM(aliases, j*k + 1, alias_v2);
        }

        INIT_FREQ(item_const_names[i], base_const+i, aliases);
    }

    return 0;
}

/* take a dictionary with integer keys and tuples of strings for values,
   and populate a dictionary with all the strings as keys and integers
   for values */
static int reverse_dict(PyObject *source, PyObject *dest) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(source, &pos, &key, &value)) {
        PyObject *tuple_iter;
        PyObject *item;

        if((tuple_iter = PyObject_GetIter(value)) == NULL) return INT_ERR_CODE;

        while ((item = PyIter_Next(tuple_iter)) != NULL) {
            PyDict_SetItem(dest, item, key);
            Py_DECREF(item);
        }
        Py_DECREF(tuple_iter);
    }
    return 0;
}



static int build_freq_dict(void) {

    char ANN_prefixes[8][15] = { "A", "Y", "ANN", "ANNUAL", "ANNUALLY",
                                 "YR", "YEAR", "YEARLY" };

    char QTRE_prefixes[8][15] = { "Q", "QTR", "QUARTER", "QUARTERLY", "Q-E",
                                  "QTR-E", "QUARTER-E", "QUARTERLY-E"};
    char QTRS_prefixes[4][15] = { "Q-S", "QTR-S", "QUARTER-S", "QUARTERLY-S" };

    char WK_prefixes[4][15] =  { "W", "WK", "WEEK", "WEEKLY" };

    /* Note: order of this array must match up with how the Annual
       frequency constants are lined up */
    char month_names[12][2][10] = {
        { "DEC", "DECEMBER" },
        { "JAN", "JANUARY" },
        { "FEB", "FEBRUARY" },
        { "MAR", "MARCH" },
        { "APR", "APRIL" },
        { "MAY", "MAY" },
        { "JUN", "JUNE" },
        { "JUL", "JULY" },
        { "AUG", "AUGUST" },
        { "SEP", "SEPTEMBER" },
        { "OCT", "OCTOBER" },
        { "NOV", "NOVEMBER" }};

    char day_names[7][2][10] = {
        { "SUN", "SUNDAY" },
        { "MON", "MONDAY" },
        { "TUE", "TUESDAY" },
        { "WED", "WEDNESDAY" },
        { "THU", "THURSDAY" },
        { "FRI", "FRIDAY" },
        { "SAT", "SATURDAY" }};

    char ANN_const_names[12][15] = {
        "FR_ANNDEC",
        "FR_ANNJAN",
        "FR_ANNFEB",
        "FR_ANNMAR",
        "FR_ANNAPR",
        "FR_ANNMAY",
        "FR_ANNJUN",
        "FR_ANNJUL",
        "FR_ANNAUG",
        "FR_ANNSEP",
        "FR_ANNOCT",
        "FR_ANNNOV"};

    char QTRE_const_names[12][15] = {
        "FR_QTREDEC",
        "FR_QTREJAN",
        "FR_QTREFEB",
        "FR_QTREMAR",
        "FR_QTREAPR",
        "FR_QTREMAY",
        "FR_QTREJUN",
        "FR_QTREJUL",
        "FR_QTREAUG",
        "FR_QTRESEP",
        "FR_QTREOCT",
        "FR_QTRENOV"};

    char QTRS_const_names[12][15] = {
        "FR_QTRSDEC",
        "FR_QTRSJAN",
        "FR_QTRSFEB",
        "FR_QTRSMAR",
        "FR_QTRSAPR",
        "FR_QTRSMAY",
        "FR_QTRSJUN",
        "FR_QTRSJUL",
        "FR_QTRSAUG",
        "FR_QTRSSEP",
        "FR_QTRSOCT",
        "FR_QTRSNOV"};

    char WK_const_names[7][15] = {
        "FR_WKSUN",
        "FR_WKMON",
        "FR_WKTUE",
        "FR_WKWED",
        "FR_WKTHU",
        "FR_WKFRI",
        "FR_WKSAT"};

    PyObject *aliases;

    freq_dict = PyDict_New();
    freq_dict_rev = PyDict_New();
    freq_constants = PyDict_New();

    aliases = Py_BuildValue("(ssss)", "M", "MTH", "MONTH", "MONTHLY");
    INIT_FREQ("FR_MTH", FR_MTH, aliases);

    aliases = Py_BuildValue("(ssss)", "B", "BUS", "BUSINESS", "BUSINESSLY");
    INIT_FREQ("FR_BUS", FR_BUS, aliases);

    aliases = Py_BuildValue("(ssss)", "D", "DAY", "DLY", "DAILY");
    INIT_FREQ("FR_DAY", FR_DAY, aliases);

    aliases = Py_BuildValue("(sssss)", "H", "HR", "HOUR", "HRLY", "HOURLY");
    INIT_FREQ("FR_HR", FR_HR, aliases);

    aliases = Py_BuildValue("(ssss)", "T", "MIN", "MINUTE", "MINUTELY");
    INIT_FREQ("FR_MIN", FR_MIN, aliases);

    aliases = Py_BuildValue("(ssss)", "S", "SEC", "SECOND", "SECONDLY");
    INIT_FREQ("FR_SEC", FR_SEC, aliases);

    aliases = Py_BuildValue("(ssss)", "U", "UND", "UNDEF", "UNDEFINED");
    INIT_FREQ("FR_UND", FR_UND, aliases);

    ADD_FREQ_CONSTANT("FR_ANN", FR_ANN);

    if(init_freq_group(12, 8, FR_ANN,
        month_names, ANN_prefixes, ANN_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    ADD_FREQ_CONSTANT("FR_QTR", FR_QTR);

    if(init_freq_group(12, 8, FR_QTREDEC,
        month_names, QTRE_prefixes, QTRE_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    if(init_freq_group(12, 4, FR_QTRSDEC,
        month_names, QTRS_prefixes, QTRS_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    ADD_FREQ_CONSTANT("FR_WK", FR_WK);

    if(init_freq_group(7, 4, FR_WK,
                    day_names, WK_prefixes, WK_const_names) == INT_ERR_CODE) {
            return INT_ERR_CODE;
    }

    if(reverse_dict(freq_dict, freq_dict_rev) == INT_ERR_CODE) {
        return INT_ERR_CODE;
    }

    return 0;
}


/* take user specified frequency and convert to int representation
   of the frequency */
int check_freq(PyObject *freq_spec) {

    if (PyInt_Check(freq_spec)) {
        return (int)PyInt_AsLong(freq_spec);
    } else if (PyString_Check(freq_spec)) {
        char *freq_str, *freq_str_uc;
        PyObject *freq_val;

        freq_str = PyString_AsString(freq_spec);
        if((freq_str_uc = str_uppercase(freq_str)) == NULL) {return INT_ERR_CODE;}

        freq_val = PyDict_GetItemString(freq_dict_rev, freq_str_uc);

        free(freq_str_uc);

        if (freq_val == NULL) {
            PyErr_SetString(PyExc_ValueError, "invalid frequency specification");
            return INT_ERR_CODE;
        } else {
            int ret_val = (int)PyInt_AsLong(freq_val);
            return ret_val;
        }
    } else if (freq_spec == Py_None) {
        return FR_UND;
    } else {
        int retval = (int)PyInt_AsLong(freq_spec);
        if (PyErr_Occurred()) {
            PyErr_SetString(PyExc_ValueError, "invalid frequency specification");
            return INT_ERR_CODE;
        } else {
            return retval; }
    }

}




PyObject *
c_freqs_check_freq(PyObject *self, PyObject *args) {
    PyObject *freq;
    int freq_val;

    if (!PyArg_ParseTuple(args, "O:check_freq(freq)", &freq)) 
        return NULL;
    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) 
        return NULL;
    return PyInt_FromLong(freq_val);
}


PyObject *
c_freqs_check_freq_str(PyObject *self, PyObject *args) {
    PyObject *alias_tuple, *result, *freq_key;

    if ((freq_key = c_freqs_check_freq(self, args)) == NULL) 
        return NULL;

    alias_tuple = PyDict_GetItem(freq_dict, freq_key);
    result = PyTuple_GET_ITEM(alias_tuple, 0);

    Py_INCREF(result);
    Py_DECREF(freq_key);

    return result;
}

PyObject *
c_freqs_get_freq_group(PyObject *self, PyObject *args) {
    PyObject *freq;
    int freq_val;
    if (!PyArg_ParseTuple(args, "O:get_freq_group(freq)", &freq)) 
        return NULL;
    if ((freq_val = check_freq(freq)) == INT_ERR_CODE) 
        return NULL;
    return PyInt_FromLong(get_base_unit(freq_val));
}


void import_c_freqs(PyObject *m)
{
    if(build_freq_dict() == INT_ERR_CODE) {
        PyErr_SetString(PyExc_ImportError,              \
                        "initialization of module timeseries.c_dates failed");
        return;
    };

    PyModule_AddObject(m, "freq_dict", freq_dict);
    PyModule_AddObject(m, "freq_dict_rev", freq_dict_rev);
    PyModule_AddObject(m, "freq_constants", freq_constants);

}

