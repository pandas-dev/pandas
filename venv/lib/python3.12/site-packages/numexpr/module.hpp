#ifndef NUMEXPR_MODULE_HPP
#define NUMEXPR_MODULE_HPP

// Deal with the clunky numpy import mechanism
// by inverting the logic of the NO_IMPORT_ARRAY symbol.
#define PY_ARRAY_UNIQUE_SYMBOL numexpr_ARRAY_API
#ifndef DO_NUMPY_IMPORT_ARRAY
#  define NO_IMPORT_ARRAY
#endif

#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <numpy/arrayscalars.h>

#include "numexpr_config.hpp"

struct global_state {
    /* Global variables for threads */
    int nthreads;                    /* number of desired threads in pool */
    int init_threads_done;           /* pool of threads initialized? */
    int end_threads;                 /* should exisiting threads end? */
    // pthread_t threads[MAX_THREADS];  /* opaque structure for threads */
    // int tids[MAX_THREADS];           /* ID per each thread */
    /* NOTE: threads and tids are arrays, they MUST be allocated to length
       `global_max_threads` before module load. */
    pthread_t *threads;              /* opaque structure for threads */
    int *tids;                       /* ID per each thread */
    npy_intp gindex;                 /* global index for all threads */
    int init_sentinels_done;         /* sentinels initialized? */
    int giveup;                      /* should parallel code giveup? */
    int force_serial;                /* force serial code instead of parallel? */
    int pid;                         /* the PID for this process */

    /* Synchronization variables for threadpool state */
    pthread_mutex_t count_mutex;
    int count_threads;
    int barrier_passed;         /* indicates if the thread pool's thread barrier
                                   is unlocked and ready for the VM to process.*/
    pthread_mutex_t count_threads_mutex;
    pthread_cond_t count_threads_cv;

    /* Mutual exclusion for access to global thread params (th_params) */
    pthread_mutex_t parallel_mutex;

    global_state() {
        nthreads = 1;
        init_threads_done = 0;
        barrier_passed = 0;
        end_threads = 0;
        pid = 0;
    }
};

extern global_state gs;

int numexpr_set_nthreads(int nthreads_new);

#endif // NUMEXPR_MODULE_HPP
