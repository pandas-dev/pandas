// Numexpr - Fast numerical array expression evaluator for NumPy.
//
//      License: MIT
//      Author:  See AUTHORS.txt
//
//  See LICENSE.txt for details about copyright and rights to use.
//
// module.cpp contains the CPython-specific module exposure.

#define DO_NUMPY_IMPORT_ARRAY

#include "module.hpp"
#include <structmember.h>
#include <vector>

#include <signal.h>

#include "interpreter.hpp"
#include "numexpr_object.hpp"

using namespace std;

// Global state. The file interpreter.hpp also has some global state
// in its 'th_params' variable
global_state gs;
long global_max_threads=DEFAULT_MAX_THREADS;

/* Do the worker job for a certain thread */
void *th_worker(void *tidptr)
{
    int tid = *(int *)tidptr;
    /* Parameters for threads */
    npy_intp start;
    npy_intp vlen;
    npy_intp block_size;
    NpyIter *iter;
    vm_params params;
    int *pc_error;
    int ret;
    int n_inputs;
    int n_constants;
    int n_temps;
    size_t memsize;
    char **mem;
    npy_intp *memsteps;
    npy_intp istart, iend;
    char **errmsg;
    // For output buffering if needed
    vector<char> out_buffer;

    while (1) {

        /* Sentinels have to be initialised yet */
        if (tid == 0) {
            gs.init_sentinels_done = 0;
        }

        /* Meeting point for all threads (wait for initialization) */
        pthread_mutex_lock(&gs.count_threads_mutex);
        if (gs.count_threads < gs.nthreads) {
            gs.count_threads++;
            /* Beware of spurious wakeups. See issue pydata/numexpr#306. */
            do {
                pthread_cond_wait(&gs.count_threads_cv,
                                  &gs.count_threads_mutex);
            } while (!gs.barrier_passed);
        }
        else {
            gs.barrier_passed = 1;
            pthread_cond_broadcast(&gs.count_threads_cv);
        }
        pthread_mutex_unlock(&gs.count_threads_mutex);

        /* Check if thread has been asked to return */
        if (gs.end_threads) {
            return(0);
        }

        /* Get parameters for this thread before entering the main loop */
        start = th_params.start;
        vlen = th_params.vlen;
        block_size = th_params.block_size;
        params = th_params.params;
        pc_error = th_params.pc_error;

        // If output buffering is needed, allocate it
        if (th_params.need_output_buffering) {
            out_buffer.resize(params.memsizes[0] * BLOCK_SIZE1);
            params.out_buffer = &out_buffer[0];
        } else {
            params.out_buffer = NULL;
        }

        /* Populate private data for each thread */
        n_inputs = params.n_inputs;
        n_constants = params.n_constants;
        n_temps = params.n_temps;
        memsize = (1+n_inputs+n_constants+n_temps) * sizeof(char *);
        /* XXX malloc seems thread safe for POSIX, but for Win? */
        mem = (char **)malloc(memsize);
        memcpy(mem, params.mem, memsize);

        errmsg = th_params.errmsg;

        params.mem = mem;

        /* Loop over blocks */
        pthread_mutex_lock(&gs.count_mutex);
        if (!gs.init_sentinels_done) {
            /* Set sentinels and other global variables */
            gs.gindex = start;
            istart = gs.gindex;
            iend = istart + block_size;
            if (iend > vlen) {
                iend = vlen;
            }
            gs.init_sentinels_done = 1;  /* sentinels have been initialised */
            gs.giveup = 0;            /* no giveup initially */
        } else {
            gs.gindex += block_size;
            istart = gs.gindex;
            iend = istart + block_size;
            if (iend > vlen) {
                iend = vlen;
            }
        }
        /* Grab one of the iterators */
        iter = th_params.iter[tid];
        if (iter == NULL) {
            th_params.ret_code = -1;
            gs.giveup = 1;
        }
        memsteps = th_params.memsteps[tid];
        /* Get temporary space for each thread */
        ret = get_temps_space(params, mem, BLOCK_SIZE1);
        if (ret < 0) {
            /* Propagate error to main thread */
            th_params.ret_code = ret;
            gs.giveup = 1;
        }
        pthread_mutex_unlock(&gs.count_mutex);

        while (istart < vlen && !gs.giveup) {
            /* Reset the iterator to the range for this task */
            ret = NpyIter_ResetToIterIndexRange(iter, istart, iend,
                                                errmsg);
            /* Execute the task */
            if (ret >= 0) {
                ret = vm_engine_iter_task(iter, memsteps, params, pc_error, errmsg);
            }

            if (ret < 0) {
                pthread_mutex_lock(&gs.count_mutex);
                gs.giveup = 1;
                /* Propagate error to main thread */
                th_params.ret_code = ret;
                pthread_mutex_unlock(&gs.count_mutex);
                break;
            }

            pthread_mutex_lock(&gs.count_mutex);
            gs.gindex += block_size;
            istart = gs.gindex;
            iend = istart + block_size;
            if (iend > vlen) {
                iend = vlen;
            }
            pthread_mutex_unlock(&gs.count_mutex);
        }

        /* Meeting point for all threads (wait for finalization) */
        pthread_mutex_lock(&gs.count_threads_mutex);
        if (gs.count_threads > 0) {
            gs.count_threads--;
            do {
                pthread_cond_wait(&gs.count_threads_cv,
                                  &gs.count_threads_mutex);
            } while (gs.barrier_passed);
        }
        else {
            gs.barrier_passed = 0;
            pthread_cond_broadcast(&gs.count_threads_cv);
        }
        pthread_mutex_unlock(&gs.count_threads_mutex);

        /* Release resources */
        free_temps_space(params, mem);
        free(mem);

    }  /* closes while(1) */

    /* This should never be reached, but anyway */
    return(0);
}

/* Initialize threads */
int init_threads(void)
{
    int tid, rc;

    if ( !(gs.nthreads > 1 && (!gs.init_threads_done || gs.pid != getpid())) ) {
        /* Thread pool must always be initialized once and once only. */
        return(0);
    }

    /* Initialize mutex and condition variable objects */
    pthread_mutex_init(&gs.count_mutex, NULL);
    pthread_mutex_init(&gs.parallel_mutex, NULL);

    /* Barrier initialization */
    pthread_mutex_init(&gs.count_threads_mutex, NULL);
    pthread_cond_init(&gs.count_threads_cv, NULL);
    gs.count_threads = 0;      /* Reset threads counter */
    gs.barrier_passed = 0;

    /*
     * Our worker threads should not deal with signals from the rest of the
     * application - mask everything temporarily in this thread, so our workers
     * can inherit that mask
     */
    sigset_t sigset_block_all, sigset_restore;
    rc = sigfillset(&sigset_block_all);
    if (rc != 0) {
        fprintf(stderr, "ERROR; failed to block signals: sigfillset: %s",
                strerror(rc));
        exit(-1);
    }
    rc = pthread_sigmask( SIG_BLOCK, &sigset_block_all, &sigset_restore);
    if (rc != 0) {
        fprintf(stderr, "ERROR; failed to block signals: pthread_sigmask: %s",
                strerror(rc));
        exit(-1);
    }

    /* Now create the threads */
    for (tid = 0; tid < gs.nthreads; tid++) {
        gs.tids[tid] = tid;
        rc = pthread_create(&gs.threads[tid], NULL, th_worker,
                            (void *)&gs.tids[tid]);
        if (rc) {
            fprintf(stderr,
                    "ERROR; return code from pthread_create() is %d\n", rc);
            fprintf(stderr, "\tError detail: %s\n", strerror(rc));
            exit(-1);
        }
    }

    /*
     * Restore the signal mask so the main thread can process signals as
     * expected
     */
    rc = pthread_sigmask( SIG_SETMASK, &sigset_restore, NULL);
    if (rc != 0) {
        fprintf(stderr,
                "ERROR: failed to restore signal mask: pthread_sigmask: %s",
                strerror(rc));
        exit(-1);
    }

    gs.init_threads_done = 1;                 /* Initialization done! */
    gs.pid = (int)getpid();                   /* save the PID for this process */

    return(0);
}

/* Set the number of threads in numexpr's VM */
int numexpr_set_nthreads(int nthreads_new)
{
    int nthreads_old = gs.nthreads;
    int t, rc;
    void *status;

    // if (nthreads_new > MAX_THREADS) {
    //     fprintf(stderr,
    //             "Error.  nthreads cannot be larger than MAX_THREADS (%d)",
    //             MAX_THREADS);
    //     return -1;
    // }
    if (nthreads_new > global_max_threads) {
        fprintf(stderr,
                "Error.  nthreads cannot be larger than environment variable \"NUMEXPR_MAX_THREADS\" (%ld)",
                global_max_threads);
        return -1;
    }
    else if (nthreads_new <= 0) {
        fprintf(stderr, "Error.  nthreads must be a positive integer");
        return -1;
    }

    /* Only join threads if they are not initialized or if our PID is
       different from that in pid var (probably means that we are a
       subprocess, and thus threads are non-existent). */
    if (gs.nthreads > 1 && gs.init_threads_done && gs.pid == getpid()) {
        /* Tell all existing threads to finish */
        gs.end_threads = 1;
        pthread_mutex_lock(&gs.count_threads_mutex);
        if (gs.count_threads < gs.nthreads) {
            gs.count_threads++;
            do {
                pthread_cond_wait(&gs.count_threads_cv,
                                  &gs.count_threads_mutex);
            } while (!gs.barrier_passed);
        }
        else {
            gs.barrier_passed = 1;
            pthread_cond_broadcast(&gs.count_threads_cv);
        }
        pthread_mutex_unlock(&gs.count_threads_mutex);

        /* Join exiting threads */
        for (t=0; t<gs.nthreads; t++) {
            rc = pthread_join(gs.threads[t], &status);
            if (rc) {
                fprintf(stderr,
                        "ERROR; return code from pthread_join() is %d\n",
                        rc);
                fprintf(stderr, "\tError detail: %s\n", strerror(rc));
                exit(-1);
            }
        }
        gs.init_threads_done = 0;
        gs.end_threads = 0;
    }

    /* Launch a new pool of threads (if necessary) */
    gs.nthreads = nthreads_new;
    init_threads();

    return nthreads_old;
}


#ifdef USE_VML

static PyObject *
_get_vml_version(PyObject *self, PyObject *args)
{
    int len=198;
    char buf[198];
    mkl_get_version_string(buf, len);
    return Py_BuildValue("s", buf);
}

static PyObject *
_set_vml_accuracy_mode(PyObject *self, PyObject *args)
{
    int mode_in, mode_old;
    if (!PyArg_ParseTuple(args, "i", &mode_in))
    return NULL;
    mode_old = vmlGetMode() & VML_ACCURACY_MASK;
    vmlSetMode((mode_in & VML_ACCURACY_MASK) | VML_ERRMODE_IGNORE );
    return Py_BuildValue("i", mode_old);
}

static PyObject *
_set_vml_num_threads(PyObject *self, PyObject *args)
{
    int max_num_threads;
    if (!PyArg_ParseTuple(args, "i", &max_num_threads))
    return NULL;
    mkl_domain_set_num_threads(max_num_threads, MKL_DOMAIN_VML);
    Py_RETURN_NONE;
}

static PyObject *
_get_vml_num_threads(PyObject *self, PyObject *args)
{
    int max_num_threads = mkl_domain_get_max_threads (MKL_DOMAIN_VML);
    return Py_BuildValue("i", max_num_threads);
}

#endif

static PyObject*
Py_set_num_threads(PyObject *self, PyObject *args)
{
    int num_threads, nthreads_old;
    if (!PyArg_ParseTuple(args, "i", &num_threads))
    return NULL;
    nthreads_old = numexpr_set_nthreads(num_threads);
    return Py_BuildValue("i", nthreads_old);
}

static PyObject*
Py_get_num_threads(PyObject *self, PyObject *args)
{
    int n_thread;
    n_thread = gs.nthreads;
    return Py_BuildValue("i", n_thread);
}

static PyMethodDef module_methods[] = {
#ifdef USE_VML
    {"_get_vml_version", _get_vml_version, METH_VARARGS,
     "Get the VML/MKL library version."},
    {"_set_vml_accuracy_mode", _set_vml_accuracy_mode, METH_VARARGS,
     "Set accuracy mode for VML functions."},
    {"_set_vml_num_threads", _set_vml_num_threads, METH_VARARGS,
     "Suggests a maximum number of threads to be used in VML operations."},
    {"_get_vml_num_threads", _get_vml_num_threads, METH_VARARGS,
     "Gets the maximum number of threads to be used in VML operations."},
#endif
    {"_set_num_threads", Py_set_num_threads, METH_VARARGS,
     "Suggests a maximum number of threads to be used in operations."},
    {"_get_num_threads", Py_get_num_threads, METH_VARARGS,
     "Gets the maximum number of threads currently in use for operations."},
    {NULL}
};

static int
add_symbol(PyObject *d, const char *sname, int name, const char* routine_name)
{
    PyObject *o, *s;
    int r;

    if (!sname) {
        return 0;
    }

    o = PyLong_FromLong(name);
    s = PyBytes_FromString(sname);
    if (!o || !s) {
        PyErr_SetString(PyExc_RuntimeError, routine_name);
        r = -1;
    }
    else {
        r = PyDict_SetItem(d, s, o);
    }
    Py_XDECREF(o);
    Py_XDECREF(s);
    return r;
}

#ifdef __cplusplus
extern "C" {
#endif

/* XXX: handle the "global_state" state via moduledef */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "interpreter",
        NULL,
        -1,                 /* sizeof(struct global_state), */
        module_methods,
        NULL,
        NULL,               /* module_traverse, */
        NULL,               /* module_clear, */
        NULL
};

#define INITERROR return NULL

PyObject *
PyInit_interpreter(void) {
    PyObject *m, *d;


    char *max_thread_str = getenv("NUMEXPR_MAX_THREADS");
    char *end;
    if (max_thread_str != NULL) {
        global_max_threads = strtol(max_thread_str, &end, 10);
    }



    th_params.memsteps    = (npy_intp**)calloc(sizeof(npy_intp*), global_max_threads);
    th_params.iter        = (NpyIter**)calloc(sizeof(NpyIter*),   global_max_threads);
    th_params.reduce_iter = (NpyIter**)calloc(sizeof(NpyIter*),   global_max_threads);
    gs.threads            = (pthread_t*)calloc(sizeof(pthread_t), global_max_threads);
    gs.tids               = (int*)calloc(sizeof(int),             global_max_threads);
    // TODO: for Py3, deallocate: https://docs.python.org/3/c-api/module.html#c.PyModuleDef.m_free
    // For Python 2.7, people have to exit the process to reclaim the memory.

    if (PyType_Ready(&NumExprType) < 0)
        INITERROR;

    m = PyModule_Create(&moduledef);

    if (m == NULL)
        INITERROR;

    #ifdef Py_GIL_DISABLED
        PyUnstable_Module_SetGIL(m, Py_MOD_GIL_NOT_USED);
    #endif

    Py_INCREF(&NumExprType);
    PyModule_AddObject(m, "NumExpr", (PyObject *)&NumExprType);

    import_array();

    d = PyDict_New();
    if (!d) INITERROR;

#define OPCODE(n, name, sname, ...)                              \
    if (add_symbol(d, sname, name, "add_op") < 0) { INITERROR; }
#include "opcodes.hpp"
#undef OPCODE

    if (PyModule_AddObject(m, "opcodes", d) < 0) INITERROR;

    d = PyDict_New();
    if (!d) INITERROR;

#define add_func(name, sname)                           \
    if (add_symbol(d, sname, name, "add_func") < 0) { INITERROR; }
#define FUNC_FF(name, sname, ...)  add_func(name, sname);
#define FUNC_FFF(name, sname, ...) add_func(name, sname);
#define FUNC_DD(name, sname, ...)  add_func(name, sname);
#define FUNC_BF(name, sname, ...)  add_func(name, sname);
#define FUNC_BD(name, sname, ...)  add_func(name, sname);
#define FUNC_BC(name, sname, ...)  add_func(name, sname);
#define FUNC_DDD(name, sname, ...) add_func(name, sname);
#define FUNC_CC(name, sname, ...)  add_func(name, sname);
#define FUNC_CCC(name, sname, ...) add_func(name, sname);
#define FUNC_II(name, sname, ...) add_func(name, sname);
#define FUNC_LL(name, sname, ...) add_func(name, sname);
#include "functions.hpp"
#undef FUNC_LL
#undef FUNC_II
#undef FUNC_CCC
#undef FUNC_CC
#undef FUNC_DDD
#undef FUNC_BC
#undef FUNC_BD
#undef FUNC_BF
#undef FUNC_DD
#undef FUNC_FFF
#undef FUNC_FF
#undef add_func

    if (PyModule_AddObject(m, "funccodes", d) < 0) INITERROR;

    if (PyModule_AddObject(m, "allaxes", PyLong_FromLong(255)) < 0) INITERROR;
    if (PyModule_AddObject(m, "maxdims", PyLong_FromLong(NPY_MAXDIMS)) < 0) INITERROR;

    if(PyModule_AddIntConstant(m, "MAX_THREADS", global_max_threads) < 0) INITERROR;

    // Let's export the block sizes to Python side for benchmarking comparisons
    if(PyModule_AddIntConstant(m, "__BLOCK_SIZE1__", BLOCK_SIZE1) < 0) INITERROR;
    // Export if we are using VML or not
#ifdef USE_VML
    if(PyModule_AddObject(m, "use_vml", Py_True) < 0) INITERROR;
#else
    if(PyModule_AddObject(m, "use_vml", Py_False) < 0) INITERROR;
#endif

    return m;
}

#ifdef __cplusplus
}  // extern "C"
#endif
