#include "c_freqs.h"
#include "c_convert.h"
#include "c_dates.h"
#include "c_datearray.h"
#include <datetime.h>
#include <time.h>

#include "c_lib.h"


static PyTypeObject DatetimeArray_Type;

int PyArray_TS_DATETIME;

#define TS_METADATA_DTSTR "timeunit"

typedef struct {
   PyObject_HEAD;
   ts_datetime obval;
   ts_metadata obmeta;
} DatetimeScalarObject;

//NPY_NO_EXPORT PyTypeObject DatetimeArrType_Type = {
//#if defined(NPY_PY3K)
//    PyVarObject_HEAD_INIT(NULL, 0)
//#else
//    PyObject_HEAD_INIT(NULL)
//    0,                                          /* ob_size */
//#endif
////    "timeseries.datetime" _THIS_SIZE,                  /* tp_name*/
//    "timeseries.datetime",                  /* tp_name*/
//    sizeof(DatetimeScalarObject),               /* tp_basicsize*/
//    0,                                          /* tp_itemsize */
//    0,                                          /* tp_dealloc */
//    0,                                          /* tp_print */
//    0,                                          /* tp_getattr */
//    0,                                          /* tp_setattr */
//#if defined(NPY_PY3K)
//    0,                                          /* tp_reserved */
//#else
//    0,                                          /* tp_compare */
//#endif
//    0,                                          /* tp_repr */
//    0,                                          /* tp_as_number */
//    0,                                          /* tp_as_sequence */
//    0,                                          /* tp_as_mapping */
//    0,                                          /* tp_hash */
//    0,                                          /* tp_call */
//    0,                                          /* tp_str */
//    0,                                          /* tp_getattro */
//    0,                                          /* tp_setattro */
//    0,                                          /* tp_as_buffer */
//    0,                                          /* tp_flags */
//    0,                                          /* tp_doc */
//    0,                                          /* tp_traverse */
//    0,                                          /* tp_clear */
//    0,                                          /* tp_richcompare */
//    0,                                          /* tp_weaklistoffset */
//    0,                                          /* tp_iter */
//    0,                                          /* tp_iternext */
//    0,                                          /* tp_methods */
//    0,                                          /* tp_members */
//    0,                                          /* tp_getset */
//    0,                                          /* tp_base */
//    0,                                          /* tp_dict */
//    0,                                          /* tp_descr_get */
//    0,                                          /* tp_descr_set */
//    0,                                          /* tp_dictoffset */
//    0,                                          /* tp_init */
//    0,                                          /* tp_alloc */
//    0,                                          /* tp_new */
//    0,                                          /* tp_free */
//    0,                                          /* tp_is_gc */
//    0,                                          /* tp_bases */
//    0,                                          /* tp_mro */
//    0,                                          /* tp_cache */
//    0,                                          /* tp_subclasses */
//    0,                                          /* tp_weaklist */
//    0,                                          /* tp_del */
//#if PY_VERSION_HEX >= 0x02060000
//    0,                                          /* tp_version_tag */
//#endif
//};
//
//#undef _THIS_SIZE
///**/




#if PY_VERSION_HEX >= 0x02070000
#define get_metadata_from_descr(descr)  \
    ((descr->metadata == NULL) ? \
     NULL :                                       \
    ((ts_metadata *)(PyCapsule_GetPointer(                   \
                     PyDict_GetItemString(descr->metadata, TS_METADATA_DTSTR), \
                     NULL))))
#else
#define get_metadata_from_descr(descr)  \
    ((descr->metadata == NULL) ? \
     NULL :                                       \
     ((ts_metadata *)(PyCObject_AsVoidPtr(                    \
                      PyDict_GetItemString(descr->metadata, TS_METADATA_DTSTR)))))
#endif

#define asarray(self) ( ((PyArrayObject *)self) )
#define get_base(self) ( ((PyArrayObject *)self)->base )
#define asndarray(self) ( asarray(get_base(self)) )
#define get_descr(self) ( ((PyArrayObject *)self)->descr )
#define get_metadata_from_array(self) (get_metadata_from_descr(get_descr(self)))
//#define get_timestep(self) (get_metadata_from_array(self)->timestep)


#define TS_METADATA_DTSTR "timeunit"

//----------------------------------------------------------------------------
/* from private/npy_3kcompat.h */
#if PY_VERSION_HEX >= 0x02070000

static NPY_INLINE PyObject *
NpyCapsule_FromVoidPtr(void *ptr, void (*dtor)(PyObject *))
{
    PyObject *ret = PyCapsule_New(ptr, NULL, dtor);
    if (ret == NULL) {
        PyErr_Clear();
    }
    return ret;
}

static void
simple_capsule_dtor(PyObject *cap)
{
    PyArray_free(PyCapsule_GetPointer(cap, NULL));
}

#else

static NPY_INLINE PyObject *
NpyCapsule_FromVoidPtr(void *ptr, void (*dtor)(void *))
{
    return PyCObject_FromVoidPtr(ptr, dtor);
}

static void
simple_capsule_dtor(void *ptr)
{
    PyArray_free(ptr);
}

#endif
/**/


#include "numpy/noprefix.h"

static void
init_descr_metadata(PyArray_Descr *descr)
{
    ts_metadata *dt_data;
    PyObject *cobj;

    dt_data = _pya_malloc(sizeof(ts_metadata));
    dt_data->unit = FR_UND;
    dt_data->timestep = 1;
    dt_data->period_end_at = 0;
    dt_data->periods_per_day = -1;
    dt_data->secs_per_period = -1;
    dt_data->convert_to_start = 0;

/* FIXME
 * There is no error check here and no way to indicate an error
 * until the metadata turns up NULL.
 */
    cobj = NpyCapsule_FromVoidPtr((void *)dt_data, simple_capsule_dtor);
    descr->metadata = PyDict_New();
    PyDict_SetItemString(descr->metadata, TS_METADATA_DTSTR, cobj);
    Py_DECREF(cobj);

}

static void
update_descr_metadata(PyArray_Descr *descr, ts_metadata *meta) {
    PyObject *cobj;
    cobj = NpyCapsule_FromVoidPtr((void *)meta, simple_capsule_dtor);
    descr->metadata = PyDict_New();
    PyDict_SetItemString(descr->metadata, TS_METADATA_DTSTR, cobj);
    Py_DECREF(cobj);
}


//-----------------------------------------------------------------------------

static PyObject *
DatetimeArray_new(PyTypeObject *cls, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"object", "unit", "timestep", "freq", NULL};

    PyObject *obj;
    PyArrayObject *arr = NULL;
    PyObject *unit = NULL, *freq=NULL;
    int timestep = 1;
    DatetimeArrayObject *self;
    PyArray_Descr *descr;

    if(!PyArg_ParseTupleAndKeywords(args, kw,"O|OiO",kwlist,
                                    &obj,
                                    &unit, &timestep, &freq))
        return NULL;

    arr = (PyArrayObject *)PyArray_FROM_O(obj);
    if(arr == NULL)
        return NULL;
//    DEBUGPRINTF("We have an array...");

    descr = PyArray_DescrNewFromType(PyArray_INT64);
    if (descr == NULL)
        return NULL;
    Py_INCREF(descr);
    init_descr_metadata(descr);

    self = (DatetimeArrayObject *)PyArray_NewFromDescr(&DatetimeArray_Type,
                                                       descr,
                                                       arr->nd, arr->dimensions,
                                                       arr->strides,
                                                       arr->data,
                                                       arr->flags,
                                                       (PyObject *)arr);
    if(self == NULL)
        return NULL;
    Py_INCREF(arr);
    PyArray_BASE(self) = (PyObject *)arr;

    if (PyObject_SetAttrString((PyObject *)self, "dtype", (PyObject *)descr) < 0) {
        goto fail;
    }

    ts_metadata *obmeta = get_metadata_from_descr(descr);

    if (unit == NULL){
        if (freq == NULL)
            freq = PyInt_FromLong(FR_UND);
        unit = freq;
    }
    int u = check_freq(unit);
    if (u == -1)
        goto fail;
    init_metadata_from_unit(obmeta, u);

    obmeta->timestep = timestep;

    return (PyObject *)self;

 fail:
    DEBUGPRINTF("Dropping it..");
    Py_XDECREF(unit);
    Py_XDECREF((PyObject*)self);
    return NULL;
}

static DatetimeArrayObject *
DatetimeArray_new_from_array_and_unit(PyArrayObject *data, int unit)
{
    DatetimeArrayObject *self;
    PyArray_Descr *descr;
    ts_metadata *obmeta;

    descr = PyArray_DescrNewFromType(PyArray_INT64);
    if (descr == NULL)
        return NULL;
    Py_INCREF(descr);
    init_descr_metadata(descr);
    obmeta = get_metadata_from_descr(descr);
    init_metadata_from_unit(obmeta, unit);

    self = (DatetimeArrayObject *)PyArray_NewFromDescr(&DatetimeArray_Type,
                                                       descr,
                                                       data->nd,
                                                       data->dimensions,
                                                       data->strides,
                                                       data->data,
                                                       data->flags,
                                                       (PyObject *)data);
    Py_INCREF(data);
    PyArray_BASE(self) = (PyObject *)data;
    return self;
}

static void
DatetimeArray_dealloc(DatetimeArrayObject *self)
{
//    DEBUGPRINTF("Dropping cache");
//    Py_XDECREF(self->cached_vals);
    DEBUGPRINTF("Dropping object");
    self->base.ob_type->tp_free((PyObject*)self);
}

static PyObject*
DatetimeArray_finalize(DatetimeArrayObject *self, PyObject *args)
{
    DatetimeArrayObject *context;
    if(PyArg_ParseTuple(args, "O", &context))
    {
        if (DatetimeArray_Check(context)){
            DEBUGPRINTF("in context from DTA");
            PyArray_Descr *descr = get_descr(self);
            Py_INCREF(descr);
            ts_metadata *meta_context = get_metadata_from_array(context);
            update_descr_metadata(descr, meta_context);
        } else {
            DEBUGPRINTF("in context from scratch");
            init_descr_metadata(get_descr(self));
        };
        ts_timestatus default_status = {-1, -1, -1};
        self->status = default_status;
    }
    PyErr_Clear();
    DEBUGPRINTF("Returning w/ base unit %i...", get_metadata_from_array(self)->unit);
    Py_RETURN_NONE;
}



static int
_get_unit_from_descr(PyArray_Descr *descr) {
    ts_metadata *meta = get_metadata_from_descr(descr);
    return meta->unit + meta->period_end_at;
}
static int
_get_unit_from_array(DatetimeArrayObject *self) {
    ts_metadata *meta = get_metadata_from_descr(((PyArrayObject *)self)->descr);
    return meta->unit + meta->period_end_at;
}
static PyObject *
DatetimeArray_unit(DatetimeArrayObject *self){
    int unit = _get_unit_from_array(self);
    return PyInt_FromLong(unit);
}
static PyObject *
DatetimeArray_timestep(DatetimeArrayObject *self){
    ts_metadata *meta = get_metadata_from_descr(((PyArrayObject *)self)->descr);
    return PyInt_FromLong(meta->timestep);
}
static PyObject *
DatetimeArray_freqstr(DatetimeArrayObject *self) {
    PyObject *key = DatetimeArray_unit(self);
    PyObject *freq_aliases = PyDict_GetItem(freq_dict, key);
    PyObject *main_alias = PyTuple_GET_ITEM(freq_aliases, 0);
    Py_DECREF(key);
    return main_alias;
}

static PyObject *
DatetimeArray_steps(DatetimeArrayObject *self){
    PyArrayObject *steps=NULL;
    PyArrayIterObject *self_iter=NULL, *steps_iter=NULL;
    npy_intp size = PyArray_SIZE(self) - 1;

    steps = (PyArrayObject*)PyArray_ZEROS(1,
                                          &size,
                                          PyArray_INT64, 0);
    if (steps == NULL)
        goto fail;

    steps_iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)steps);
    if (steps_iter == NULL)
        goto fail;
    self_iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
    if (self_iter == NULL)
        goto fail;

    PyObject *val=NULL, *prev=NULL, *diff=NULL;
    prev = PyArray_GETITEM(self, self_iter->dataptr);
    PyArray_ITER_NEXT(self_iter);
    while (steps_iter->index < steps_iter->size) {
        val = PyArray_GETITEM(self, self_iter->dataptr);
        diff = PyNumber_Subtract(val, prev);
        PyArray_SETITEM(steps, steps_iter->dataptr, diff);
        PyArray_ITER_NEXT(self_iter);
        PyArray_ITER_NEXT(steps_iter);
        prev = val;
    };
    Py_DECREF(self_iter);
    Py_DECREF(steps_iter);
    Py_XDECREF(prev);
    Py_XDECREF(val);
    Py_XDECREF(diff);
    return (PyObject *)steps;

 fail:
    DEBUGPRINTF("DatetimeArray.steps: Oops...");
    Py_XDECREF(steps);
    Py_XDECREF(steps_iter);
    Py_XDECREF(self_iter);
    return NULL;
}



static int
DatetimeArray_check_status(DatetimeArrayObject *self)
{
    PyArrayIterObject *self_iter=NULL;
    npy_int64 timestep, diff;
    int is_chrono = 1, has_dups=0, has_missing=0;

    timestep = get_metadata_from_array(self)->timestep;
    self_iter = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
    if (self_iter == NULL) {
        Py_XDECREF(self_iter);
        return -1;
    }

    PyObject *val=NULL, *prev=NULL, *odiff=NULL;
    prev = PyArray_GETITEM(self, self_iter->dataptr);
    PyArray_ITER_NEXT(self_iter);
    while (self_iter->index < self_iter->size) {
        val = PyArray_GETITEM(self, self_iter->dataptr);
        odiff = PyNumber_Subtract(val, prev);
        diff = PyInt_AsLong(odiff);
        if (diff < 0)
            is_chrono = 0;
        else if (diff == 0)
            has_dups = 1;
        else if (diff > timestep)
            has_missing = 1;
        if (has_dups && has_missing)
            break;
        PyArray_ITER_NEXT(self_iter);
        prev = val;
    }
    Py_XDECREF(self_iter);
    Py_XDECREF(odiff);
    Py_XDECREF(prev);
    Py_XDECREF(val);
    // Set the status
//    self->status.has_dups = has_dups;
//    self->status.has_missing = has_missing;
//    self->status.is_chrono = is_chrono;
    ts_timestatus status = {has_dups, has_missing, is_chrono};
    self->status = status;
    return 0;
}
static PyObject *
DatetimeArray_has_dups(DatetimeArrayObject *self)
{
    if (self->status.has_dups == -1)
        if (DatetimeArray_check_status(self) < 0)
            return NULL;
    if (self->status.has_dups == 0)
        Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}
static PyObject *
DatetimeArray_has_missing(DatetimeArrayObject *self)
{
    if (self->status.has_missing == -1)
        if (DatetimeArray_check_status(self) < 0)
            return NULL;
    if (self->status.has_missing == 0)
        Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}
static PyObject *
DatetimeArray_is_chrono(DatetimeArrayObject *self)
{
    if (self->status.is_chrono == -1)
        if (DatetimeArray_check_status(self) < 0)
            return NULL;
    if (self->status.is_chrono == 0)
        Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}
static PyObject *
DatetimeArray_is_full(DatetimeArrayObject *self)
{
    if (self->status.has_dups == -1)
        if (DatetimeArray_check_status(self) < 0)
            return NULL;
    if (self->status.has_dups)
        Py_RETURN_FALSE;
    if (self->status.has_missing)
        Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}
static PyObject *
DatetimeArray_is_valid(DatetimeArrayObject *self)
{
    ts_timestatus status = self->status;
    if (status.has_dups == -1)
        if (DatetimeArray_check_status(self) < 0)
            return NULL;
    if (status.has_missing)
        Py_RETURN_FALSE;
    if (status.has_dups)
        Py_RETURN_FALSE;
    if (! status.is_chrono)
        Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}

static PyMemberDef DatetimeArray_members[] = {
//     {"cached_vals", T_OBJECT_EX, offsetof(DateTimeArray, cached_vals), 0,
//      "cached_values"},
    {NULL}  /* Sentinel */
};


static char *
DEBUGGETTYPE(PyObject *obj){
    char *type_str;
    PyObject *type_repr, *obj_type;
    obj_type = PyObject_Type(obj);
    type_repr = PyObject_Repr(obj_type);
    type_str = PyString_AsString(type_repr);
//    DEBUGPRINTF("get_tsdatetime_from_object got %s [%i]", type_str, meta->unit);
    Py_DECREF(obj_type);
    Py_DECREF(type_repr);
    return type_str;
}


static ts_datetime
get_tsdatetime_from_object(ts_metadata *meta, PyObject *date){
    ts_datetime value;
    //
    if (PyString_Check(date)) {
        value = PyString_to_tsdatetime(meta, date);
//        DEBUGPRINTF("get_tsdatetime_from_object.from string: %ld", value);
    }
    else if (PyDateTime_Check(date) || PyDate_Check(date)) {
        value = PyDatetime_to_tsdatetime(meta, date);
//        DEBUGPRINTF("get_tsdatetime_from_object.from datetime.datetime: %ld", value);
    }
    else if (DatetimeObject_Check(date)) {
        value = ((DatetimeObject *)date)->obval;
//        DEBUGPRINTF("get_tsdatetime_from_object.from tsdatetime: %ld", value);
    }
    else if (PyInt_Check(date) || PyLong_Check(date) || PyFloat_Check(date)) {
        value = (ts_datetime)PyInt_AsLong(date);
//        DEBUGPRINTF("get_tsdatetime_from_object.from number: %ld", value);
    }
    else {
        value = -1;
    }
    return value;
}



static PyObject *
DatetimeArray_single_date_to_index(DatetimeArrayObject *self, PyObject *date){
    intp count=0, i, size;
    int nd = ((PyArrayObject *)self)->nd, j, comparison;

    PyArrayIterObject *itr = NULL;
    PyObject *result = NULL, *item;
    intp *dptr[MAX_DIMS];

    itr = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
    if (itr == NULL)
        return NULL;

    ts_metadata *meta = get_metadata_from_array(self);
    ts_datetime value = get_tsdatetime_from_object(meta, date);
    if (value < 0) {
        goto fail;
    }

    PyArray_CompareFunc *cmprf = ((PyArrayObject *)self)->descr->f->compare;

    /*Count the valid elements*/
    size = itr->size;
    for (i = 0; i < size; i++) {
        comparison = cmprf(itr->dataptr, &value, self);
        if (comparison == 0)
            count++;
        PyArray_ITER_NEXT(itr);
    }

    PyArray_ITER_RESET(itr);
    result = PyTuple_New(nd);
    if (result == NULL)
        goto fail;
    for (j = 0; j < nd; j++) {
        item = PyArray_New(Py_TYPE(self), 1, &count,
                           PyArray_INTP, NULL, NULL, 0, 0,
                           (PyObject *)self);
        if (item == NULL)
            goto fail;
        PyTuple_SET_ITEM(result, j, item);
        dptr[j] = (intp *)PyArray_DATA(item);
    }
    if (nd == 1) {
        for (i = 0; i < size; i++){
            comparison = cmprf(itr->dataptr, &value, self);
            if (comparison == 0)
                *(dptr[0])++ = i;
            PyArray_ITER_NEXT(itr);
        }
    }
    else {
        itr->contiguous = 0;
        for (i = 0; i < size; i++){
            comparison = cmprf(itr->dataptr, &value, self);
            if (comparison == 0) {
                for (j = 0; j < nd; j++)
                    *(dptr[j])++ = itr->coordinates[j];
            }
            PyArray_ITER_NEXT(itr);
        }
    }
    Py_DECREF(itr);
    return result;
 fail:
    Py_XDECREF(result);
    Py_XDECREF(itr);
    return NULL;
}



static PyObject *
DatetimeArray_date_to_index(DatetimeArrayObject *self, PyObject *dateargs){
    PyObject *result=NULL, *date=NULL;
    ts_datetime value;
    Py_ssize_t i;

    /* Make sure we have at least 1 argument */
    Py_ssize_t nbargs = PyObject_Length(dateargs);
    if (nbargs < 1) {
        PyErr_SetString(PyExc_ValueError, "there should be at least one argument");
        goto fail;
    }

    ts_metadata *meta = get_metadata_from_array(self);
//    ts_timestatus status = self->status;
//    int is_valid = ((! status.has_missing) && (! status.has_dups) && (status.is_chrono));

    result = PyList_New(0);
    if (result == NULL)
        goto fail;

    PyArrayIterObject *itr = NULL;
    PyObject *indexlist = NULL;
    int comparison, empty;

    for (i=0; i < nbargs; i++){
        date = PyTuple_GetItem(dateargs, i);
        value = get_tsdatetime_from_object(meta, date);
        if (value < 0) {
            PyErr_SetString(PyExc_ValueError, "unable to retrieve date");
            Py_XDECREF(itr);
            Py_XDECREF(indexlist);
            goto fail;
        }

        indexlist = PyList_New(0);
        if (indexlist == NULL) {
            Py_XDECREF(itr);
            Py_XDECREF(indexlist);
            goto fail;
        }
        itr = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
        if (itr == NULL) {
            Py_XDECREF(itr);
            Py_XDECREF(indexlist);
            goto fail;
        }
        PyArray_CompareFunc *cmprf = ((PyArrayObject *)self)->descr->f->compare;
        while (itr->index < itr->size) {
            comparison = cmprf(itr->dataptr, &value, self);
            if (comparison == 0) {
                PyObject *coords = PyObject_GetAttrString((PyObject *)itr, "coords");
                PyList_Append(indexlist, coords);
                empty = 0;
            }
            PyArray_ITER_NEXT(itr);
        };
        if (empty) {
            indexlist = Py_None;
        };
        PyList_Append(result, indexlist);
    };
    Py_DECREF(itr);
    Py_DECREF(indexlist);

    if (nbargs == 1)
        return PyList_GetItem(result, 0);
    return result;

 fail:
    Py_XDECREF(result);
    Py_XDECREF(date);
    return NULL;
}





static PyObject *
DatetimeArray_getitem(DatetimeArrayObject *self, PyObject *op)
{
//    int reset_full=1, keep_chrono=0;
//    DEBUGPRINTF("in __getitem__ w %s", DEBUGGETTYPE(op));
    PyObject *idx;

    if (DatetimeObject_Check(op) || PyString_Check(op) || PyDateTime_Check(op)) {
        if (DatetimeObject_Check(op)) {
            DEBUGPRINTF("index is Date");
        }
        else if (PyString_Check(op)) {
            DEBUGPRINTF("index is string");
        }
        else if (PyDateTime_Check(op)) {
            DEBUGPRINTF("index is datetime");
        };
        idx = DatetimeArray_single_date_to_index(self, op);
        if (idx == NULL) {
            PyErr_SetString(PyExc_IndexError, "date out of bounds");
            return NULL;
        }
    }
    else {
        idx = op;
    }

    PyObject *r, *result;
    r = ((PyArrayObject *)self)->ob_type->tp_base->tp_as_mapping->mp_subscript((PyObject *)self, idx);
    if (r == NULL) {
        return NULL;
    }
//    DEBUGPRINTF("r is %s", DEBUGGETTYPE(r));
    ts_datetime obval;
    if (PyArray_IsScalar(r, Integer)) {
        int unit = _get_unit_from_descr(get_descr(self));

        obval = (ts_datetime)(PyInt_AsLong(r));
        result = (PyObject *)DatetimeObject_FromFreqAndValue(unit, PyInt_AsLong(r));
        Py_DECREF(r);
    }
    else {
        result = r;
        ((DatetimeArrayObject *)r)->status.is_chrono = self->status.is_chrono;
    }
    Py_DECREF(idx);
    return result;
}



NPY_NO_EXPORT PyMappingMethods DatetimeArray_as_mapping = {
    NULL,              /*mp_length*/
    (binaryfunc)&DatetimeArray_getitem,        /*mp_subscript*/
    NULL, /*mp_ass_subscript*/
};





/* Date & Time Information */
PyObject *
DatetimeArray_getdateinfo(DatetimeArrayObject *self, char *infochar)
{
    int skip_periods, counter=1, val_changed=0;

    PyObject *prev_val=NULL;
    PyArrayObject *output=NULL;
    PyArrayIterObject *iterin=NULL, *iterout=NULL;


    ts_metadata *meta = get_metadata_from_array(self);
    int unit = meta->unit;
    ts_timestatus status = self->status;
    int is_valid = ((! status.has_missing) && (! status.has_dups) && (status.is_chrono));

    output = (PyArrayObject *)PyArray_SimpleNew(((PyArrayObject *)self)->nd,
                                                 ((PyArrayObject *)self)->dimensions,
                                                 NPY_INT);


    conversion_function todays = get_converter_to_days(unit, 1);
    init_metadata_from_unit(meta, unit);
    meta->convert_to_start = 0;
    ts_datetimestruct dinfo;


    iterin = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
    iterout = (PyArrayIterObject *)PyArray_IterNew((PyObject *)output);

    PyObject* (*getdateparam)(npy_int64, int,
                              conversion_function, ts_metadata*,
                              ts_datetimestruct*) = NULL;
    switch(*infochar)
    {
        case 'Y': //year
            getdateparam = &_loop_get_year;
            skip_periods = __skip_periods_year(unit);
            break;
        case 'F': //"fiscal" year
            if (unit == FR_QTR)
                getdateparam = &_loop_get_qyear_from_qtr;
            else
                getdateparam = &_loop_get_qyear;
            skip_periods = __skip_periods_year(unit);
            break;
        case 'Q': //quarter
            if (unit == FR_QTR)
                getdateparam = &_loop_get_quarter_from_qtr;
            else
                getdateparam = &_loop_get_quarter;
            skip_periods = __skip_periods_quarter(unit);
            break;
        case 'M': //month
            getdateparam = &_loop_get_month;
            skip_periods = __skip_periods_month(unit);
            break;
        case 'D': //day
            getdateparam = &_loop_get_day;
            skip_periods = __skip_periods_day(unit);
            break;
        case 'R': //day of year
            getdateparam = &_loop_get_day_of_year;
            skip_periods = __skip_periods_day(unit);
            break;
        case 'W': //day of week
            getdateparam = &_loop_get_day_of_week;
            skip_periods = __skip_periods_day(unit);
            break;
        case 'I': //week of year
            getdateparam = &_loop_get_week;
            skip_periods = __skip_periods_week(unit);
            break;
        case 'H': //hour
            getdateparam = &_loop_get_hour;
            skip_periods = __skip_periods_hour(unit);
            break;
        case 'T': //minute
            getdateparam = &_loop_get_minute;
            skip_periods = __skip_periods_minute(unit);
            break;
        case 'S': //second
            getdateparam = &_loop_get_second;
            skip_periods = 1;
            break;
        case 'O': //toordinal
            getdateparam = &_loop_get_ordinal;
            skip_periods = __skip_periods_day(unit);
            break;
        default:
            return NULL;
    }

    {
    PyObject *val, *result=NULL;
    while (iterin->index < iterin->size) {

        if ((val_changed == 0) ||
            (is_valid == 0) ||
            (prev_val == NULL) ||
            (counter >= skip_periods)) {

               val = PyArray_GETITEM(self, iterin->dataptr);
               result = getdateparam(PyInt_AsLong(val), unit,
                                     todays, meta, &dinfo);

               if ((prev_val != NULL) &&
                   (PyLong_AsLong(prev_val) != PyLong_AsLong(result))) {
                   val_changed = 1;
                   counter = 0;
               }
               Py_DECREF(val);
               if (prev_val != NULL) {
                   Py_DECREF(prev_val);
               }
               prev_val = result;
        }
        PyArray_SETITEM(output, iterout->dataptr, result);

        PyArray_ITER_NEXT(iterin);
        PyArray_ITER_NEXT(iterout);
        counter++;
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
DatetimeArray_year(DatetimeArrayObject *self){
    char infochar = 'Y';
    return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_qyear(DatetimeArrayObject *self){
     char infochar = 'F';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_quarter(DatetimeArrayObject *self){
     char infochar = 'Q';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_month(DatetimeArrayObject *self){
     char infochar = 'M';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_week(DatetimeArrayObject *self){
     char infochar = 'I';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_day(DatetimeArrayObject *self){
     char infochar = 'D';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_day_of_week(DatetimeArrayObject *self){
     char infochar = 'W';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_day_of_year(DatetimeArrayObject *self){
     char infochar = 'R';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_hour(DatetimeArrayObject *self){
     char infochar = 'H';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_minute(DatetimeArrayObject *self){
     char infochar = 'T';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_second(DatetimeArrayObject *self){
     char infochar = 'S';
     return DatetimeArray_getdateinfo(self, &infochar);
}
PyObject *
DatetimeArray_ordinal(DatetimeArrayObject *self){
    char infochar = 'O';
    return DatetimeArray_getdateinfo(self, &infochar);
}


PyObject *
DatetimeArray_datetime(DatetimeArrayObject *self)
{
    PyArrayObject *output=NULL;
    PyArrayIterObject *iterin=NULL, *iterout=NULL;

    ts_metadata *meta = get_metadata_from_array(self);
    int unit = meta->unit;

    output = (PyArrayObject *)PyArray_SimpleNew(((PyArrayObject *)self)->nd,
                                                ((PyArrayObject *)self)->dimensions,
                                                NPY_OBJECT);

    conversion_function todays = get_converter_to_days(unit, 1);
    meta->convert_to_start = 1;
    ts_datetimestruct dinfo;

    iterin = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
    iterout = (PyArrayIterObject *)PyArray_IterNew((PyObject *)output);

    {
    PyObject *val, *result=NULL;
    while (iterin->index < iterin->size) {
        val = PyArray_GETITEM(self, iterin->dataptr);
        result = _loop_get_datetime(PyInt_AsLong(val), unit,
                                    todays, meta, &dinfo);
        PyArray_SETITEM(output, iterout->dataptr, result);
        PyArray_ITER_NEXT(iterin);
        PyArray_ITER_NEXT(iterout);
        }
    }
    Py_DECREF(iterin);
    Py_DECREF(iterout);
    return (PyObject *) output;
}

static PyObject *
DatetimeArray_start_date(DatetimeArrayObject *self) {
    PyObject *minobj, *result=NULL;
    minobj = PyArray_Min(asndarray(self), MAX_DIMS, NULL);
    ts_datetime val = PyInt_AsLong(minobj);
    int unit = _get_unit_from_descr(get_descr(self));
    result = (PyObject *)DatetimeObject_FromFreqAndValue(unit, val);
    Py_DECREF(minobj);
    return result;
}
static PyObject *
DatetimeArray_end_date(DatetimeArrayObject *self) {
    PyObject *maxobj, *result=NULL;
    maxobj = PyArray_Max(asndarray(self), MAX_DIMS, NULL);
    ts_datetime val = PyInt_AsLong(maxobj);
    int unit = _get_unit_from_descr(get_descr(self));
    result = (PyObject *)DatetimeObject_FromFreqAndValue(unit, val);
    Py_DECREF(maxobj);
    return result;
}


static PyObject *
DatetimeArray_tovalues(DatetimeArrayObject *self) {
    return get_base(self);
}
static PyObject *
DatetimeArray_toordinals(DatetimeArrayObject *self) {
    char infochar = 'O';
    return DatetimeArray_getdateinfo(self, &infochar);
}
static PyObject *
DatetimeArray_tolist(DatetimeArrayObject *self)
{
    PyObject *dtarray=NULL, *output=NULL;
    dtarray = DatetimeArray_datetime(self);
    output = PyArray_ToList((PyArrayObject *)self);
    Py_DECREF(dtarray);
    return output;
}


/*
 * PROPERTIES
 */

static int
DatetimeArray_ReadOnlyErr(DatetimeArrayObject *self, PyObject *value, void *closure) {
   PyErr_SetString(PyExc_AttributeError, "Cannot set read-only property");
   return -1;
};

static PyGetSetDef DatetimeArray_getseters[] = {
    {"unit", (getter)DatetimeArray_unit, (setter)DatetimeArray_ReadOnlyErr,
     "Returns the frequency.", NULL},
    {"timestep", (getter)DatetimeArray_timestep, (setter)DatetimeArray_ReadOnlyErr,
     "", NULL},
    {"freqstr", (getter)DatetimeArray_freqstr, (setter)DatetimeArray_ReadOnlyErr,
     "Returns the string representation of frequency.", NULL},
    {"steps", (getter)DatetimeArray_steps, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"year", (getter)DatetimeArray_year, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"qyear", (getter)DatetimeArray_qyear, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"quarter", (getter)DatetimeArray_quarter, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"month", (getter)DatetimeArray_month, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"week", (getter)DatetimeArray_week, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"day", (getter)DatetimeArray_day, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"day_of_week", (getter)DatetimeArray_day_of_week, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"day_of_year", (getter)DatetimeArray_day_of_year, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"hour", (getter)DatetimeArray_hour, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"minute", (getter)DatetimeArray_minute, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"second", (getter)DatetimeArray_second, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
    {"datetime", (getter)DatetimeArray_datetime, (setter)DatetimeArray_ReadOnlyErr,
     "time steps", NULL},
     {"start_date", (getter)DatetimeArray_start_date, (setter)DatetimeArray_ReadOnlyErr,
      "time steps", NULL},
     {"end_date", (getter)DatetimeArray_end_date, (setter)DatetimeArray_ReadOnlyErr,
      "time steps", NULL},
    {NULL, NULL, NULL, NULL, NULL}  /* Sentinel */
};

/*
 * METHODS
 */



static PyObject *
DatetimeArray_convert(DatetimeArrayObject *self,
                      PyObject *args, PyObject *kwds)
{
    DatetimeArrayObject *output=NULL;
    PyObject *freq=NULL;
    char *relation_raw=NULL, *relation_uc;
    conversion_function converterfrom, converterto;
    int fromfreq, tofreq;

    PyArray_Descr *indescr, *outdescr;
    ts_metadata *metafrom, *metato;

    indescr = get_descr(self);
    metafrom = get_metadata_from_descr(indescr);
    fromfreq = _get_unit_from_descr(indescr);

    /* Get the arguments */
    static char *kwlist[] = {"freq", "relation", NULL};
    if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|s", kwlist,
                                      &freq, &relation_raw))
        return NULL;

    /* Check the conversion frequency */
    if ((tofreq = check_freq(freq)) == INT_ERR_CODE)
        return NULL;


    /* Initialize the output */
    outdescr = PyArray_DescrNewFromType(PyArray_INT64);
    if (outdescr == NULL)
        return NULL;
    Py_INCREF(outdescr);
    output = (DatetimeArrayObject *)PyArray_NewFromDescr(&DatetimeArray_Type,
                                                         outdescr,
                                                         asarray(self)->nd,
                                                         asarray(self)->dimensions,
                                                         NULL,
                                                         NULL,
                                                         asarray(self)->flags,
                                                         NULL);
    if (output == NULL)
        return NULL;
    metato = get_metadata_from_descr(outdescr);
    init_metadata_from_unit(metato, tofreq);


    /* Update the convert_to_start from the relational argument */
    if(relation_raw) {
        if (strlen(relation_raw) > 0) {
            if ((relation_uc = str_uppercase(relation_raw)) == NULL)
                return PyErr_NoMemory();
            // 'BEFORE' and 'AFTER' values for this parameter are deprecated
            if ((relation_uc[0] == 'E') || (relation_uc[0] == 'A'))
                metafrom->convert_to_start = 0;
            else if ((relation_uc[0] == 'S') || (relation_uc[0] == 'B'))
                metafrom->convert_to_start = 1;
            else {
                PyErr_SetString(PyExc_ValueError,"Invalid relation specification");
                free(relation_uc);
                return NULL;
            }
            free(relation_uc);
        }
        else {
            metafrom->convert_to_start = 0;
        };
    }
    else {
        metafrom->convert_to_start = 0;
    }
    metato->convert_to_start = metafrom->convert_to_start;
    /* Correction for business days */
    if ((tofreq == FR_BUS) && (fromfreq < FR_DAY))
        metato->convert_to_start = 1;


    PyArrayIterObject *iterfrom, *iterto;
    iterfrom = (PyArrayIterObject *)PyArray_IterNew((PyObject *)get_base(self));
    iterto = (PyArrayIterObject *)PyArray_IterNew((PyObject *)output);
    if ((iterfrom == NULL) || (iterto == NULL)) {
        Py_XDECREF(iterfrom);
        Py_XDECREF(iterto);
        Py_XDECREF(output);
        return NULL;
    }

    PyObject *fromdateobj=NULL, *todateobj=NULL;
    PyArray_GetItemFunc *getitem = indescr->f->getitem;

    if (tofreq == fromfreq) {
        while (iterfrom->index < iterfrom->size) {
            fromdateobj = getitem(iterfrom->dataptr, self);
            PyArray_SETITEM(output, iterto->dataptr, fromdateobj);
            Py_DECREF(fromdateobj);
            PyArray_ITER_NEXT(iterfrom);
            PyArray_ITER_NEXT(iterto);
        }
    }
    else {
        ts_datetime fromdateval, todateval;
        converterfrom = convert_to_mediator(metafrom->unit, metato->unit, 0);
        converterto = convert_from_mediator(metafrom->unit, metato->unit, 0);
        while (iterfrom->index < iterfrom->size) {
            fromdateobj = getitem(iterfrom->dataptr, self);
            fromdateval = PyInt_AsLong(fromdateobj);
            todateval = converterto(converterfrom(fromdateval, metafrom), metato);
            todateobj = PyInt_FromLong(todateval);

            PyArray_SETITEM(output, iterto->dataptr, todateobj);
            Py_DECREF(fromdateobj);
            Py_DECREF(todateobj);

            PyArray_ITER_NEXT(iterfrom);
            PyArray_ITER_NEXT(iterto);
        }
    }
    Py_DECREF(iterfrom);
    Py_DECREF(iterto);
    return (PyObject *)output;

}


static PyObject *
DatetimeArray_fill_missing_dates(DatetimeArrayObject *self,
                                 PyObject *args, PyObject *kwds)
{
    DatetimeArrayObject *output=NULL;
    PyObject *base=NULL, *result=NULL;
    PyArray_Descr *descr = get_descr(self);

    int output_mask = 0;
    static char *kwlist[] = {"output_mask", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i", kwlist,
                                     &output_mask))
        return NULL;

    ts_metadata *meta = get_metadata_from_descr(descr);
    ts_timestatus status = self->status;

    if (status.has_dups) {
        PyErr_SetString(PyExc_ValueError, "duplicated dates are not allowed");
        goto fail;
    }
    if (! status.is_chrono) {
        PyErr_SetString(PyExc_ValueError, "series must be in chronological order");
        goto fail;
    }
    if (! status.has_missing) {
        DEBUGPRINTF("fill_missing_dates : no missing");
        base = PyArray_Copy(asndarray(self));
    }
    else {
        DEBUGPRINTF("fill_missing_dates : w/ missing");
        PyObject *start=NULL, *end=NULL;
        start = PyArray_Min(asndarray(self), MAX_DIMS, NULL);
        end = PyArray_Max(asndarray(self), MAX_DIMS, NULL);
        if ((start == NULL) || (end == NULL)) {
            Py_XDECREF(start);
            Py_XDECREF(end);
            goto fail;
        }
        PyObject *pystep = PyInt_FromLong(meta->timestep);
        base = PyArray_ArangeObj(start, PyNumber_Add(end, pystep), pystep,
                                 descr);
        Py_DECREF(pystep);
        Py_XDECREF(start);
        Py_XDECREF(end);
    }
    DEBUGPRINTF("initializing the output");
    output = (DatetimeArrayObject *)PyArray_View(asarray(base), descr, Py_TYPE(self));
    if (output == NULL) {
        goto fail;
    }
    update_descr_metadata(get_descr(output), meta);
    ts_timestatus newstatus = {0, 0, 1};
    output->status = newstatus;

    if (output_mask) {
        npy_intp size = PyArray_SIZE(output);
        PyObject *mask = PyArray_ZEROS(1, &size, NPY_BOOL, 0);
        DEBUGPRINTF("mask initialized");
        if (status.has_missing) {
            PyArray_CompareFunc *cmprf = ((PyArrayObject *)self)->descr->f->compare;
            PyArrayIterObject *iself=NULL, *ifill=NULL, *imask=NULL;
            iself = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
            ifill = (PyArrayIterObject *)PyArray_IterNew((PyObject *)output);
            imask = (PyArrayIterObject *)PyArray_IterNew((PyObject *)mask);
            while (iself->index < iself->size) {
                while (cmprf(iself->dataptr, ifill->dataptr, self) == 1) {
                    PyArray_SETITEM(mask, imask->dataptr, PyInt_FromLong(1));
                    PyArray_ITER_NEXT(ifill);
                    PyArray_ITER_NEXT(imask);
                }
                PyArray_ITER_NEXT(iself);
                PyArray_ITER_NEXT(ifill);
                PyArray_ITER_NEXT(imask);
            }
            Py_XDECREF(iself);
            Py_XDECREF(ifill);
            Py_XDECREF(imask);
        }
        DEBUGPRINTF("all ok");
        result = Py_BuildValue("(OO)", output, mask);
    }
    else
        result = Py_BuildValue("O", output);

    return result;

 fail:
    Py_XDECREF(output);
    Py_XDECREF(base);
    return NULL;
}

static PyObject *
DatetimeArray_get_missing_dates_mask(DatetimeArrayObject *self)
{
    PyObject *mask = NULL;

    ts_timestatus status = self->status;
    int is_valid = ((! status.has_missing) && (! status.has_dups) && (status.is_chrono));

    if (is_valid) {
        npy_intp size = PyArray_SIZE(self);
        mask = PyArray_ZEROS(1, &size, NPY_BOOL, 0);
        if (mask == NULL)
            goto fail;
    }
    else {
        PyObject *filled = DatetimeArray_fill_missing_dates(self, NULL, NULL);
        if (filled == NULL) {
            DEBUGPRINTF("FSCK");
            Py_XDECREF(filled);
            Py_XDECREF(mask);
            return NULL;
        }
        npy_intp size = PyArray_SIZE(get_base(filled));
        mask = PyArray_ZEROS(1, &size, NPY_BOOL, 0);
        if (mask == NULL)
            goto fail;

        DEBUGPRINTF("mask initialized");
        PyArray_CompareFunc *cmprf = ((PyArrayObject *)self)->descr->f->compare;
        PyArrayIterObject *iself=NULL, *ifill=NULL, *imask=NULL;
        iself = (PyArrayIterObject *)PyArray_IterNew((PyObject *)self);
        ifill = (PyArrayIterObject *)PyArray_IterNew((PyObject *)filled);
        imask = (PyArrayIterObject *)PyArray_IterNew((PyObject *)mask);
        int i=0;
        while (iself->index < iself->size) {
            DEBUGPRINTF("ini got %ld", PyInt_AsLong(PyArray_GETITEM(self,iself->dataptr)));
            DEBUGPRINTF("fld got 1 at %ld", PyInt_AsLong(PyArray_GETITEM(self,ifill->dataptr)));
            DEBUGPRINTF("[%i]", cmprf(iself->dataptr, ifill->dataptr, self));
            while (cmprf(iself->dataptr, ifill->dataptr, self) == 1) {
                i++;
                DEBUGPRINTF("Add 1 at %i", i);
//                DEBUGPRINTF("fld got 1 at %ld", PyInt_AsLong(PyArray_GETITEM(self,ifill->dataptr)));
                PyArray_SETITEM(mask, imask->dataptr, PyInt_FromLong(1));
                PyArray_ITER_NEXT(ifill);
                PyArray_ITER_NEXT(imask);
            }
            PyArray_ITER_NEXT(iself);
            PyArray_ITER_NEXT(ifill);
            PyArray_ITER_NEXT(imask);
        }
        Py_XDECREF(iself);
        Py_XDECREF(ifill);
        Py_XDECREF(imask);
        Py_DECREF(filled);
    }
    return (PyObject *)mask;

 fail:
    DEBUGPRINTF("FSCK");
    Py_XDECREF(mask);
    return NULL;
}


static PyMethodDef DatetimeArray_methods[] = {
    { "__array_finalize__", (PyCFunction)DatetimeArray_finalize, METH_VARARGS,
      ""},
//    {"__getitem__", (PyCFunction)DatetimeArray_getitem, METH_VARARGS, ""},
    {"has_dups", (PyCFunction)DatetimeArray_has_dups, METH_VARARGS, ""},
    {"has_missing", (PyCFunction)DatetimeArray_has_missing, METH_VARARGS, ""},
    {"is_chrono", (PyCFunction)DatetimeArray_is_chrono, METH_VARARGS, ""},
    {"is_full", (PyCFunction)DatetimeArray_is_full, METH_VARARGS, ""},
    {"is_valid", (PyCFunction)DatetimeArray_is_valid, METH_VARARGS, ""},
    {"date_to_index", (PyCFunction)DatetimeArray_date_to_index, METH_VARARGS, ""},
    {"tovalues", (PyCFunction)DatetimeArray_tovalues, METH_VARARGS, ""},
    {"toordinals", (PyCFunction)DatetimeArray_toordinals, METH_VARARGS, ""},
    {"tolist", (PyCFunction)DatetimeArray_tolist, METH_VARARGS, ""},
    {"fill_missing_dates", (PyCFunction)DatetimeArray_fill_missing_dates, METH_KEYWORDS, ""},
    {"get_missing_dates_mask", (PyCFunction)DatetimeArray_get_missing_dates_mask, METH_VARARGS, ""},
    {"convert", (PyCFunction)DatetimeArray_convert, METH_KEYWORDS, ""},
    {0}
};


static PyTypeObject DatetimeArray_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                        /* ob_size */
    "timeseries.DatetimeArray",                      /* tp_name */
    sizeof(DatetimeArrayObject),              /* tp_basicsize */
    0,                                        /* tp_itemsize */
    (destructor)DatetimeArray_dealloc,          /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    &DatetimeArray_as_mapping,                 /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    0,                                        /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "DatetimeArray",                          /* tp_doc */
    0,                                        /* tp_traverse */
    0,                                        /* tp_clear */
    0,     /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    0,                                        /* tp_iter */
    0,                                        /* tp_iternext */
    DatetimeArray_methods,                    /* tp_methods */
    DatetimeArray_members,                    /* tp_members */
    DatetimeArray_getseters,                  /* tp_getset */
    0,                            /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    0,                                        /* tp_init */
    0,                                        /* tp_alloc */
    DatetimeArray_new,                        /* tp_new */
};






/*
 * */
void import_c_datearray(PyObject *m)
{
    import_array();
    PyDateTime_IMPORT;

    DatetimeArray_Type.tp_base = &PyArray_Type;
    if (PyType_Ready(&DatetimeArray_Type) < 0)
        return;
    Py_INCREF(&DatetimeArray_Type);
    PyModule_AddObject(m, "DatetimeArray", (PyObject *)(&DatetimeArray_Type));
    
//    PyArray_TS_DATETIME = PyArray_RegisterDataType(&TS_DATETIME_Descr);
//    if (PyArray_TS_DATETIME < 0) {
//        DEBUGPRINTF("Could not import the TS_DATETIME description.");
//        return;
//    };
//    TS_DATETIME_Descr.ob_type = &PyArrayDescr_Type;
//    Py_INCREF(&TS_DATETIME_Descr);

    // PyModule_AddObject(m, "Datetime", (PyObject *)(&TS_DATETIME_Descr));

}


