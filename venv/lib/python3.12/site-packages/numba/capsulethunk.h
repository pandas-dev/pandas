/**
 
 This is a modified version of capsulethunk.h for use in llvmpy 
 
**/

#ifndef __CAPSULETHUNK_H
#define __CAPSULETHUNK_H

#if (    (PY_VERSION_HEX <  0x02070000) \
     || ((PY_VERSION_HEX >= 0x03000000) \
      && (PY_VERSION_HEX <  0x03010000)) )

//#define Assert(X) do_assert(!!(X), #X, __FILE__, __LINE__)
#define Assert(X)

static
void do_assert(int cond, const char * msg, const char *file, unsigned line){
    if (!cond) {
        fprintf(stderr, "Assertion failed %s:%d\n%s\n", file, line, msg);
        exit(1);
    }
}

typedef void (*PyCapsule_Destructor)(PyObject *);

struct FakePyCapsule_Desc {
    const char *name;
    void *context;
    PyCapsule_Destructor dtor;
    PyObject *parent;

    FakePyCapsule_Desc() : name(0), context(0), dtor(0) {}
};

static
FakePyCapsule_Desc* get_pycobj_desc(PyObject *p){
    void *desc = ((PyCObject*)p)->desc;
    Assert(desc && "No desc in PyCObject");
    return static_cast<FakePyCapsule_Desc*>(desc);
}

static
void pycobject_pycapsule_dtor(void *p, void *desc){
    Assert(desc);
    Assert(p);
    FakePyCapsule_Desc *fpc_desc = static_cast<FakePyCapsule_Desc*>(desc);
    Assert(fpc_desc->parent);
    Assert(PyCObject_Check(fpc_desc->parent));
    fpc_desc->dtor(static_cast<PyObject*>(fpc_desc->parent));
    delete fpc_desc;
}

static
PyObject* PyCapsule_New(void* ptr, const char *name, PyCapsule_Destructor dtor)
{
    FakePyCapsule_Desc *desc = new FakePyCapsule_Desc;
    desc->name = name;
    desc->context = NULL;
    desc->dtor = dtor;
    PyObject *p = PyCObject_FromVoidPtrAndDesc(ptr, desc,
                                               pycobject_pycapsule_dtor);
    desc->parent = p;
    return p;
}

static
int PyCapsule_CheckExact(PyObject *p)
{
    return PyCObject_Check(p);
}

static
void* PyCapsule_GetPointer(PyObject *p, const char *name)
{
    Assert(PyCapsule_CheckExact(p));
    if (strcmp(get_pycobj_desc(p)->name, name) != 0) {
        PyErr_SetString(PyExc_ValueError, "Invalid PyCapsule object");
    }
    return PyCObject_AsVoidPtr(p);
}

static
void* PyCapsule_GetContext(PyObject *p)
{
    Assert(p);
    Assert(PyCapsule_CheckExact(p));
    return get_pycobj_desc(p)->context;
}

static
int PyCapsule_SetContext(PyObject *p, void *context)
{
    Assert(PyCapsule_CheckExact(p));
    get_pycobj_desc(p)->context = context;
    return 0;
}

static
const char * PyCapsule_GetName(PyObject *p)
{
//    Assert(PyCapsule_CheckExact(p));
    return get_pycobj_desc(p)->name;
}

#endif /* #if PY_VERSION_HEX < 0x02070000 */

#endif /* __CAPSULETHUNK_H */
