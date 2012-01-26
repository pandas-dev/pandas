#include "c_freqs.h"
#include "c_convert.h"
#include "c_dates.h"
#include "c_tseries.h"

/* Helper function for TimeSeries_convert:
    determine the size of the second dimension for the resulting
    converted array */
static long get_width(int fromFreq, int toFreq) {

    int maxBusDaysPerYear = 262;
    int maxBusDaysPerQuarter = 66;
    int maxBusDaysPerMonth = 23;

    int maxDaysPerYear = 366;
    int maxDaysPerQuarter = 92;
    int maxDaysPerMonth = 31;

    int fromGroup = get_base_unit(fromFreq);
    int toGroup = get_base_unit(toFreq);

    if (fromGroup == FR_UND) { fromGroup = FR_DAY; }

    switch(fromGroup)
    {
        case FR_ANN: return 1;
        case FR_QTR:
            switch(toGroup)
            {
                case FR_ANN: return 4;
                default: return 1;
            }
        case FR_MTH: //monthly
            switch(toGroup)
            {
                case FR_ANN: return 12;
                case FR_QTR: return 3;
                default: return 1;
            }
        case FR_WK: //weekly
            switch(toGroup)
            {
                case FR_ANN: return 53;
                case FR_QTR: return 13;
                case FR_MTH: return 4;
                default: return 1;
            }
        case FR_BUS: //business
            switch(toGroup)
            {
                case FR_ANN: return maxBusDaysPerYear;;
                case FR_QTR: return maxBusDaysPerQuarter;
                case FR_MTH: return maxBusDaysPerMonth;
                case FR_WK: return 5;
                default: return 1;
            }
        case FR_DAY: //daily
            switch(toGroup)
            {
                case FR_ANN: return maxDaysPerYear;;
                case FR_QTR: return maxDaysPerQuarter;
                case FR_MTH: return maxDaysPerMonth;
                case FR_WK: return 7;
                default: return 1;
            }
        case FR_HR: //hourly
            switch(toGroup)
            {
                case FR_ANN: return 24 * maxDaysPerYear;;
                case FR_QTR: return 24 * maxDaysPerQuarter;
                case FR_MTH: return 24 * maxDaysPerMonth;
                case FR_WK: return 24 * 7;
                case FR_DAY: return 24;
                case FR_BUS: return 24;
                default: return 1;
            }
        case FR_MIN: //minutely
            switch(toGroup)
            {
                case FR_ANN: return 24 * 60 * maxDaysPerYear;;
                case FR_QTR: return 24 * 60 * maxDaysPerQuarter;
                case FR_MTH: return 24 * 60 * maxDaysPerMonth;
                case FR_WK: return 24 * 60 * 7;
                case FR_DAY: return 24 * 60;
                case FR_BUS: return 24 * 60;
                case FR_HR: return 60;
                default: return 1;
            }
        case FR_SEC: //minutely
            switch(toGroup)
            {
                case FR_ANN: return 24 * 60 * 60 * maxDaysPerYear;
                case FR_QTR: return 24 * 60 * 60 * maxDaysPerQuarter;
                case FR_MTH: return 24 * 60 * 60 * maxDaysPerMonth;
                case FR_WK: return 24 * 60 * 60 * 7;
                case FR_DAY: return 24 * 60 * 60;
                case FR_BUS: return 24 * 60 * 60;
                case FR_HR: return 60 * 60;
                case FR_MIN: return 60;
                default: return 1;
            }
        default: return 1;
    }
}

PyObject *
TimeSeries_convert(PyObject *self, PyObject *args)
{
    PyObject *arrayTest;
    PyArrayObject *array, *newArray;
    PyArrayObject *mask, *newMask;

    PyObject *returnVal = NULL;
    PyObject *start_index_retval;

    long period;
    long startIndex, endIndex;
    npy_int64 newStart, newStartTemp;
    npy_int64 newEnd, newEndTemp;
    long newLen, newWidth;
    long currIndex, prevIndex;
    long nd;
    npy_intp *dim, *newIdx;
    long currPerLen=0;
    char *position;
    PyObject *fromFreq_arg, *toFreq_arg;
    int fromFreq, toFreq;
    char relation_from, relation_to;
    int i;
    conversion_function totmp, fromtmp;
    ts_metadata metato, metafrom;

    PyObject *val, *valMask;

    // long (*asfreq_main)(long, char, asfreq_info*) = NULL;
    // long (*asfreq_endpoints)(long, char, asfreq_info*) = NULL;
    // long (*asfreq_reverse)(long, char, asfreq_info*) = NULL;

    returnVal = PyDict_New();

    if (!PyArg_ParseTuple(args,
        "OOlOslO:convert(array, fromfreq, period, tofreq, position, startindex, mask)",
        &array, &fromFreq_arg, &period, &toFreq_arg,
        &position, &startIndex, &mask)) return NULL;

    if((fromFreq = check_freq(fromFreq_arg)) == INT_ERR_CODE)
        return NULL;
    if((toFreq = check_freq(toFreq_arg)) == INT_ERR_CODE)
        return NULL;

    if (toFreq == fromFreq) {
        PyObject *sidx;
        newArray = (PyArrayObject *)PyArray_Copy(array);
        newMask = (PyArrayObject *)PyArray_Copy(mask);
        sidx = PyInt_FromLong(startIndex);

        PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
        PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);
        PyDict_SetItemString(returnVal, "startindex", sidx);

        Py_DECREF(newArray);
        Py_DECREF(newMask);
        Py_DECREF(sidx);

        return returnVal;
    }

    switch(position[0])
    {
        case 'S':
            // start -> before
            relation_to = 'S';
            break;
        case 'E':
            // end -> after
            relation_to = 'E';
            break;
        default:
            return NULL;
            break;
    }
    if ((toFreq == FR_BUS) && (fromFreq < FR_DAY))
        relation_from = 'S';
    else
        relation_from = relation_to;

    totmp = convert_to_mediator(fromFreq, toFreq, 1);
    init_metadata_from_unit(&metato, fromFreq);
    metato.convert_to_start = 1;
//    set_conversion_info(fromFreq, 'S', &infoto);
    fromtmp = convert_from_mediator(fromFreq, toFreq, 1);
    init_metadata_from_unit(&metafrom, toFreq);
    metafrom.convert_to_start = 1;
//    set_conversion_info(toFreq, 'S', &infofrom);



    // get_asfreq_info(fromFreq, toFreq, &af_info);

    // asfreq_main = get_asfreq_func(fromFreq, toFreq, 1);
    // asfreq_endpoints = get_asfreq_func(fromFreq, toFreq, 0);

    //convert start index to new frequency
    ERR_CHECK(newStartTemp = fromtmp(totmp(startIndex, &metato), &metafrom));
// asfreq_main(startIndex, 'S', &af_info));
    newStart = newStartTemp;
//    if (newStartTemp < 1) {
//        ERR_CHECK(newStart = asfreq_endpoints(startIndex, 'E', &af_info));
//    } else {
//        newStart = newStartTemp;
//    };
//    if (newStart < 1) {
//        PyErr_SetString(PyExc_ValueError,
//                        "start_date outside allowable range for destination frequency");
//        return NULL;
//    };

    //convert end index to new frequency
    endIndex = startIndex + (array->dimensions[0] - 1)*period;

    metato.convert_to_start = (int)0;
    metafrom.convert_to_start = (int)0;
    ERR_CHECK(newEndTemp = fromtmp(totmp(endIndex, &metato), &metafrom));
    // ERR_CHECK(newEndTemp = asfreq_main(endIndex, 'E', &af_info));
//    if (newEndTemp < 1) {
//        ERR_CHECK(newEnd = asfreq_endpoints(endIndex, 'S', &af_info));
//    } else { newEnd = newEndTemp; }
    newEnd = newEndTemp;
    newLen = newEnd - newStart + 1;
    newWidth = get_width(fromFreq, toFreq);
    if (newWidth % period > 0){
        newWidth = newWidth / period + 1;
    } else {
        newWidth /= period;
    }

    if (newWidth > 1) {
        long tempval;
        conversion_function totmprev, fromtmprev;
        ts_metadata metatorev, metafromrev;


        // get_asfreq_info(toFreq, fromFreq, &af_info_rev);
        // asfreq_reverse = get_asfreq_func(toFreq, fromFreq, 0);
        totmprev = convert_to_mediator(toFreq, fromFreq, 0);
        init_metadata_from_unit(&metatorev, toFreq);
        metatorev.convert_to_start = 1;
//        set_conversion_info(toFreq, 'S', &metatorev);
        fromtmprev = convert_from_mediator(toFreq, fromFreq, 0);
//        set_conversion_info(fromFreq, 'S', &infofromrev);
        init_metadata_from_unit(&metafromrev, fromFreq);
        metafromrev.convert_to_start = 1;

        // ERR_CHECK(tempval = asfreq_reverse(newStart, 'S', &af_info_rev));
        ERR_CHECK(tempval = fromtmprev(totmprev(newStart, &metatorev), &metafromrev));
        currPerLen = startIndex - tempval;

        nd = 2;
        dim = PyDimMem_NEW(nd);
        dim[0] = (npy_intp)newLen;
        dim[1] = (npy_intp)newWidth;
    } else {
        nd = 1;
        dim = PyDimMem_NEW(nd);
        dim[0] = (npy_intp)newLen;
    }

    newIdx = PyDimMem_NEW(nd);
    arrayTest = PyArray_SimpleNew(nd, dim, array->descr->type_num);
    if (arrayTest == NULL) { return NULL; }
    newArray = (PyArrayObject*)arrayTest;
    newMask  = (PyArrayObject*)PyArray_SimpleNew(nd, dim, mask->descr->type_num);

    PyDimMem_FREE(dim);

    PyArray_FILLWBYTE(newArray,0);
    PyArray_FILLWBYTE(newMask,1);

    prevIndex = newStart;

    metafrom.convert_to_start = (relation_from == 'S');
    metato.convert_to_start = (relation_to == 'S');


    //set values in the new array

    for (i = 0; i < array->dimensions[0]; i++) {

        npy_intp idx = (npy_intp)i;

        val = PyArray_GETITEM(array, PyArray_GetPtr(array, &idx));
        valMask = PyArray_GETITEM(mask, PyArray_GetPtr(mask, &idx));

        // ERR_CHECK(currIndex = asfreq_main(startIndex + i*period, relation, &af_info));
        ERR_CHECK(currIndex = fromtmp(totmp(startIndex + i*period, &metato), &metafrom));

        newIdx[0] = (npy_intp)(currIndex-newStart);

        if (newWidth > 1) {
            if (currIndex != prevIndex) {
                //reset period length
                currPerLen = 0;
                prevIndex = currIndex;
            }
            newIdx[1] = (npy_intp)currPerLen;
            currPerLen++;
        }

        if (newIdx[0] > -1) {
            PyArray_SETITEM(newArray, PyArray_GetPtr(newArray, newIdx), val);
            PyArray_SETITEM(newMask, PyArray_GetPtr(newMask, newIdx), valMask);
        }

        Py_DECREF(val);
        Py_DECREF(valMask);

    }

    PyDimMem_FREE(newIdx);

    start_index_retval = (PyObject*)PyInt_FromLong(newStart);

    PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
    PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);
    PyDict_SetItemString(returnVal, "startindex", start_index_retval);

    Py_DECREF(newArray);
    Py_DECREF(newMask);
    Py_DECREF(start_index_retval);

    return returnVal;
}


/* This function is directly copied from the numpy source  */
/* Return typenumber from dtype2 unless it is NULL, then return
   NPY_DOUBLE if dtype1->type_num is integer or bool
   and dtype1->type_num otherwise.
*/
static int
_get_type_num_double(PyArray_Descr *dtype1, PyArray_Descr *dtype2)
{
    if (dtype2 != NULL) {
        return dtype2->type_num;
    }

    /* For integer or bool data-types */
    if (dtype1->type_num < NPY_FLOAT) {
        return NPY_DOUBLE;
    }
    else {
        return dtype1->type_num;
    }
}

static int
_get_type_num(PyArray_Descr *dtype1, PyArray_Descr *dtype2)
{
    if (dtype2 != NULL) {
        return dtype2->type_num;
    } else {
        return dtype1->type_num;
    }
}


/* validates the standard arguments to moving functions and set the original
   mask, original ndarray, and mask for the result */
static PyObject *
check_mov_args(
    PyObject *orig_arrayobj, int span, int min_win_size,
    PyObject **orig_ndarray, PyObject **orig_mask, PyObject **result_mask
) {

    PyArrayObject **orig_ndarray_tmp, **result_mask_tmp;
    int *raw_result_mask;

    if (!PyArray_Check(orig_arrayobj)) {
        PyErr_SetString(PyExc_ValueError, "array must be a valid subtype of ndarray");
        return NULL;
    }

    // check if array has a mask, and if that mask is an array
    if (PyObject_HasAttrString(orig_arrayobj, "_mask")) {
        PyObject *tempMask = PyObject_GetAttrString(orig_arrayobj, "_mask");
        if (PyArray_Check(tempMask)) {
            *orig_mask = PyArray_EnsureArray(tempMask);
        } else {
            Py_DECREF(tempMask);
        }
    }

    *orig_ndarray = PyArray_EnsureArray(orig_arrayobj);
    orig_ndarray_tmp = (PyArrayObject**)orig_ndarray;

    if ((*orig_ndarray_tmp)->nd != 1) {
        PyErr_SetString(PyExc_ValueError, "array must be 1 dimensional");
        return NULL;
    }

    if (span < min_win_size) {
        char *error_str;
        error_str = PyArray_malloc(60 * sizeof(char));
        MEM_CHECK(error_str);
        sprintf(error_str,
                "span must be greater than or equal to %i",
                min_win_size);
        PyErr_SetString(PyExc_ValueError, error_str);
        free(error_str);
        return NULL;
    }

    raw_result_mask = PyArray_malloc((*orig_ndarray_tmp)->dimensions[0] * sizeof(int));
    MEM_CHECK(raw_result_mask);

    {
        PyArrayObject *orig_mask_tmp;
        int i, valid_points=0, is_masked;

        orig_mask_tmp = (PyArrayObject*)(*orig_mask);

        for (i=0; i<((*orig_ndarray_tmp)->dimensions[0]); i++) {

            npy_intp idx = (npy_intp)i;
            is_masked=0;

            if (*orig_mask != NULL) {
                PyObject *valMask;
                valMask = PyArray_GETITEM(orig_mask_tmp,
                                          PyArray_GetPtr(orig_mask_tmp, &idx));
                is_masked = (int)PyInt_AsLong(valMask);
                Py_DECREF(valMask);
            }

            if (is_masked) {
                valid_points=0;
            } else {
                if (valid_points < span) { valid_points += 1; }
                if (valid_points < span) { is_masked = 1; }
            }

            raw_result_mask[i] = is_masked;
        }
    }

    *result_mask = PyArray_SimpleNewFromData(
                             1, (*orig_ndarray_tmp)->dimensions,
                             PyArray_INT32, raw_result_mask);
    MEM_CHECK(*result_mask);
    result_mask_tmp = (PyArrayObject**)result_mask;
    (*result_mask_tmp)->flags = ((*result_mask_tmp)->flags) | NPY_OWNDATA;
    return 0;
}

// check if value at specified index is masked
static int
_is_masked(PyArrayObject *mask, npy_intp idx) {

    if (mask != NULL) {
        PyObject *val_mask;
        int is_masked;

        val_mask = PyArray_GETITEM(mask, PyArray_GetPtr(mask, &idx));
        is_masked = (int)PyInt_AsLong(val_mask);
        Py_DECREF(val_mask);
        return is_masked;
    } else {
        return 0;
    }

}

/* computation portion of moving sum. Appropriate mask is overlayed on top
   afterwards */
static PyObject*
calc_mov_sum(
    PyArrayObject *orig_ndarray, PyArrayObject *orig_mask, int span, int rtype)
{
    PyArrayObject *result_ndarray=NULL;
    int i=0, non_masked=0;

    result_ndarray = (PyArrayObject*)PyArray_ZEROS(
                                       orig_ndarray->nd,
                                       orig_ndarray->dimensions,
                                       rtype, 0);
    NULL_CHECK(result_ndarray);

    for (i=0; i<orig_ndarray->dimensions[0]; i++) {

        PyObject *val=NULL, *mov_sum_val=NULL;
        npy_intp idx = (npy_intp)i;
        int curr_val_masked;

        curr_val_masked = _is_masked(orig_mask, idx);

        val = PyArray_GETITEM(
            orig_ndarray, PyArray_GetPtr(orig_ndarray, &idx));

        if (curr_val_masked == 0) {
            non_masked += 1;
        } else {
            non_masked = 0;
        }

        if (
            ((i == 0) || (curr_val_masked == 1)) ||
            ((i > 0) && (_is_masked(orig_mask, i-1) == 1))
        ) {
            // if current or previous value is masked, reset moving sum
            mov_sum_val = val;
        } else {
            PyObject *mov_sum_prevval;

            idx = (npy_intp)(i-1);
            mov_sum_prevval= PyArray_GETITEM(result_ndarray,
                                   PyArray_GetPtr(result_ndarray, &idx));
            mov_sum_val = np_add(val, mov_sum_prevval);
            Py_DECREF(mov_sum_prevval);
            NULL_CHECK(mov_sum_val);

            if (non_masked > span) {
                PyObject *temp_val, *rem_val;
                idx = (npy_intp)(i-span);

                if (_is_masked(orig_mask, idx) == 0) {
                    // don't subtract off old value if it was masked because it
                    // is not included in moving sum

                    temp_val = mov_sum_val;
                    rem_val = PyArray_GETITEM(orig_ndarray,
                                       PyArray_GetPtr(orig_ndarray, &idx));

                    mov_sum_val = np_subtract(temp_val, rem_val);
                    NULL_CHECK(mov_sum_val);

                    Py_DECREF(temp_val);
                    Py_DECREF(rem_val);
                }
            }
        }

        idx = (npy_intp)i;

        PyArray_SETITEM(result_ndarray,
                        PyArray_GetPtr(result_ndarray, &idx),
                        mov_sum_val);

        if (mov_sum_val != val) { Py_DECREF(val); }

        Py_DECREF(mov_sum_val);
    }

    return (PyObject*)result_ndarray;

}

PyObject *
MaskedArray_mov_sum(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL, *orig_mask=NULL,
             *result_ndarray=NULL, *result_mask=NULL,
             *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span, type_num_double;

    static char *kwlist[] = {"array", "span", "type_num_double", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oii|O&:mov_sum(array, span, type_num_double , dtype)", kwlist,
                &orig_arrayobj, &span, &type_num_double,
                PyArray_DescrConverter2, &dtype)) return NULL;

    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &orig_mask, &result_mask);

    if (type_num_double) {
        /* if the moving sum is being used as an intermediate step in something
        like a standard deviation calculation, etc... then _get_type_num_double
        should be used to determine the appropriate return type. */
        rtype = _get_type_num_double(((PyArrayObject*)orig_ndarray)->descr, dtype);
    } else {
        rtype = _get_type_num(((PyArrayObject*)orig_ndarray)->descr, dtype);
    }

    result_ndarray = calc_mov_sum(
        (PyArrayObject*)orig_ndarray, (PyArrayObject*)orig_mask,
        span, rtype
    );
    NULL_CHECK(result_ndarray);

    result_dict = PyDict_New();
    MEM_CHECK(result_dict);
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

PyObject*
calc_mov_ranked(PyArrayObject *orig_ndarray, int span, int rtype, char rank_type)
{
    PyArrayObject *result_ndarray=NULL;
    PyObject **result_array, **ref_array, **even_array=NULL;
    PyObject *new_val, *old_val;
    PyObject *temp_add, *one_half;
    int a, i, k, R, arr_size, z;
    int *r;
    npy_intp idx;

    arr_size = (int)(orig_ndarray->dimensions[0]);

    result_ndarray = (PyArrayObject*)PyArray_ZEROS(
                                       orig_ndarray->nd,
                                       orig_ndarray->dimensions,
                                       rtype, 0);
    NULL_CHECK(result_ndarray);

    if (arr_size >= span) {
        result_array = calloc(arr_size, sizeof(PyObject*));
        MEM_CHECK(result_array);

        /* this array will be used for quick access to the data in the original
           array (so PyArray_GETITEM doesn't have to be used over and over in the
           main loop) */
        ref_array = PyArray_malloc(arr_size * sizeof(PyObject*));
        MEM_CHECK(ref_array);

        for (i=0; i<arr_size; i++) {
            idx = (npy_intp)i;
            ref_array[i] = PyArray_GETITEM(orig_ndarray, PyArray_GetPtr(orig_ndarray, &idx));
        }

        /* this array wll be used for keeping track of the "ranks" of the values
           in the current window */
        r = PyArray_malloc(span * sizeof(int));
        MEM_CHECK(r);

        for (i=0; i < span; i++) {
            r[i] = 1;
        }

        if (rank_type == 'E' && ((span % 2) == 0)) {
            // array to store two median values when span is an even #
            even_array = calloc(2, sizeof(PyObject*));
            MEM_CHECK(even_array);
        }

        switch(rank_type) {
            case 'E': // median
                R = (span + 1)/2;
                break;
            case 'I': // min
                R = 1;
                break;
            case 'A': // max
                R = span;
                break;
            default:
            {
                PyErr_SetString(PyExc_RuntimeError, "unexpected rank type");
                return NULL;
            }
        }

        one_half = PyFloat_FromDouble(0.5);

        z = arr_size - span;

        /* Calculate initial ranks "r" */
        for (i=0; i < span; i++) {

            for (k=0;   k < i;  k++) {
                if (np_greater_equal(ref_array[z+i], ref_array[z+k])) {
                    r[i]++;
                }
            }
            for (k=i+1; k < span; k++) {
                if (np_greater(ref_array[z+i], ref_array[z+k])) {
                    r[i]++;
                }
            }

            /* If rank=R, this is the median */
            if (even_array != NULL) {
                if (r[i]==R) {
                    even_array[0] = ref_array[z+i];
                } else if (r[i] == (R+1)) {
                    even_array[1] = ref_array[z+i];
                }
            } else {
                if (r[i]==R) {
                    result_array[arr_size-1] = ref_array[z+i];
                }
            }
        }

        if (even_array != NULL) {
            temp_add = np_add(even_array[0], even_array[1]);
            result_array[arr_size-1] = np_multiply(temp_add, one_half);
            Py_DECREF(temp_add);
        }

        for (i=arr_size-2; i >= span-1; i--) {
            a = span;
            z = i - span + 1;
            old_val = ref_array[i+1];
            new_val = ref_array[i-span+1];

            for (k=span-1; k > 0; k--) {
                r[k] = r[k-1]; /* Shift previous iteration's ranks */
                if (np_greater_equal(ref_array[z+k], new_val)) {r[k]++; a--;}
                if (np_greater(ref_array[z+k], old_val)) {r[k]--;}

                if (r[k]==R) {
                    result_array[i] = ref_array[z+k];
                }

                if (even_array != NULL) {
                    if (r[k]==R) {
                        even_array[0] = ref_array[z+k];
                    } else if (r[k] == (R+1)) {
                        even_array[1] = ref_array[z+k];
                    }
                } else {
                    if (r[k]==R) {
                        result_array[i] = ref_array[z+k];
                    }
                }

            }

            r[0] = a;

            if (even_array != NULL) {
                if (a==R) {
                    even_array[0] = new_val;
                } else if (a == (R+1)) {
                    even_array[1] = new_val;
                }

                temp_add = np_add(even_array[0], even_array[1]);
                result_array[i] = np_multiply(temp_add, one_half);;
                Py_DECREF(temp_add);

            } else {
                if (a==R) {
                    result_array[i] = new_val;
                }
            }

        }

        Py_DECREF(one_half);

        for (i=span-1; i<arr_size; i++) {
            idx = (npy_intp)i;
            PyArray_SETITEM(result_ndarray,
                            PyArray_GetPtr(result_ndarray, &idx),
                            result_array[i]);
        }

        for (i=0; i<arr_size; i++) {
            Py_DECREF(ref_array[i]);
        }

        if (even_array != NULL) {
            for (i=span-1; i<arr_size; i++) {
                Py_DECREF(result_array[i]);
            }
        }

        free(ref_array);
        free(result_array);
    }

    return (PyObject*)result_ndarray;

}

PyObject *
MaskedArray_mov_median(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL, *orig_mask=NULL,
             *result_ndarray=NULL, *result_mask=NULL, *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_median(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;

    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &orig_mask, &result_mask);

    if ((span % 2) == 0) {
        rtype = _get_type_num_double(((PyArrayObject*)orig_ndarray)->descr, dtype);
    } else {
        rtype = _get_type_num(((PyArrayObject*)orig_ndarray)->descr, dtype);
    }

    result_ndarray = calc_mov_ranked((PyArrayObject*)orig_ndarray,
                                     span, rtype, 'E');
    NULL_CHECK(result_ndarray);

    result_dict = PyDict_New();
    MEM_CHECK(result_dict);
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

PyObject *
MaskedArray_mov_min(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL, *orig_mask=NULL,
             *result_ndarray=NULL, *result_mask=NULL, *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_min(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;

    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &orig_mask, &result_mask);

    rtype = _get_type_num(((PyArrayObject*)orig_ndarray)->descr, dtype);

    result_ndarray = calc_mov_ranked((PyArrayObject*)orig_ndarray,
                                     span, rtype, 'I');
    NULL_CHECK(result_ndarray);

    result_dict = PyDict_New();
    MEM_CHECK(result_dict);
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

PyObject *
MaskedArray_mov_max(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL, *orig_mask=NULL,
             *result_ndarray=NULL, *result_mask=NULL, *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_max(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;

    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &orig_mask, &result_mask);

    rtype = _get_type_num(((PyArrayObject*)orig_ndarray)->descr, dtype);

    result_ndarray = calc_mov_ranked((PyArrayObject*)orig_ndarray,
                                     span, rtype, 'A');
    NULL_CHECK(result_ndarray);

    result_dict = PyDict_New();
    MEM_CHECK(result_dict);
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

/* computation portion of exponentially weighted moving average. Appropriate
   mask is overlayed on top afterwards */
static PyObject*
calc_mov_average_expw(
    PyArrayObject *orig_ndarray, PyArrayObject *orig_mask, int span, int rtype)
{
    PyArrayObject *result_ndarray=NULL;
    PyObject *decay_factor=NULL;
    int i=0, initialized=0;

    result_ndarray = (PyArrayObject*)PyArray_ZEROS(
                                       orig_ndarray->nd,
                                       orig_ndarray->dimensions,
                                       rtype, 0);
    NULL_CHECK(result_ndarray);

    decay_factor = PyFloat_FromDouble(2.0/((double)(span + 1)));

    for (i=0; i<orig_ndarray->dimensions[0]; i++) {

        PyObject *val=NULL, *mov_avg_val=NULL;
        npy_intp idx = (npy_intp)i;
        int curr_val_masked;

        curr_val_masked = _is_masked(orig_mask, idx);

        val = PyArray_GETITEM(
            orig_ndarray, PyArray_GetPtr(orig_ndarray, &idx));

        if (initialized == 0) {
            mov_avg_val = val;
            if (curr_val_masked == 0) {
                initialized = 1;
            }
        } else {
            PyObject *mov_avg_prevval, *temp_val_a, *temp_val_b;
            idx = (npy_intp)(i-1);
            mov_avg_prevval = PyArray_GETITEM(result_ndarray,
                               PyArray_GetPtr(result_ndarray, &idx));

            if (curr_val_masked == 0) {
                temp_val_a = np_subtract(val, mov_avg_prevval);
                temp_val_b = np_multiply(decay_factor, temp_val_a);
                mov_avg_val = np_add(mov_avg_prevval, temp_val_b);

                Py_DECREF(mov_avg_prevval);
                Py_DECREF(temp_val_a);
                Py_DECREF(temp_val_b);
                NULL_CHECK(mov_avg_val);
            } else {
                mov_avg_val = mov_avg_prevval;
            }
        }

        idx = (npy_intp)i;

        PyArray_SETITEM(result_ndarray,
                        PyArray_GetPtr(result_ndarray, &idx),
                        mov_avg_val);

        if (mov_avg_val != val) { Py_DECREF(val); }

        Py_DECREF(mov_avg_val);
    }

    return (PyObject*)result_ndarray;

}

PyObject *
MaskedArray_mov_average_expw(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL, *orig_mask=NULL,
             *result_ndarray=NULL, *result_mask=NULL,
             *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_average_expw(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;

    // note: we do not actually use the "result_mask" in this case
    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &orig_mask, &result_mask);

    rtype = _get_type_num_double(((PyArrayObject*)orig_ndarray)->descr, dtype);

    result_ndarray = calc_mov_average_expw(
        (PyArrayObject*)orig_ndarray, (PyArrayObject*)orig_mask,
        span, rtype
    );
    NULL_CHECK(result_ndarray);

    result_dict = PyDict_New();
    MEM_CHECK(result_dict);
    PyDict_SetItemString(result_dict, "array", result_ndarray);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

void import_c_tseries(PyObject *m) { import_array(); }
