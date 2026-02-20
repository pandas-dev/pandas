#ifndef NUMEXPR_INTERPRETER_HPP
#define NUMEXPR_INTERPRETER_HPP

#include "numexpr_config.hpp"

// Forward declaration
struct NumExprObject;

enum OpCodes {
#define OPCODE(n, e, ...) e = n,
#include "opcodes.hpp"
#undef OPCODE
};

enum FuncFFCodes {
#define FUNC_FF(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_FF
};

enum FuncBFCodes {
#define FUNC_BF(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_BF
};

enum FuncFFFCodes {
#define FUNC_FFF(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_FFF
};

enum FuncDDCodes {
#define FUNC_DD(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_DD
};

enum FuncBDCodes {
#define FUNC_BD(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_BD
};

enum FuncBCCodes {
#define FUNC_BC(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_BC
};

enum FuncIICodes {
#define FUNC_II(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_II
};

enum FuncLLCodes {
#define FUNC_LL(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_LL
};

enum FuncDDDCodes {
#define FUNC_DDD(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_DDD
};

enum FuncCCCodes {
#define FUNC_CC(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_CC
};

enum FuncCCCCodes {
#define FUNC_CCC(fop, ...) fop,
#include "functions.hpp"
#undef FUNC_CCC
};

struct vm_params {
    int prog_len;
    unsigned char *program;
    int n_inputs;
    int n_constants;
    int n_temps;
    unsigned int r_end;
    char *output;
    char **inputs;
    char **mem;
    npy_intp *memsteps;
    npy_intp *memsizes;
    struct index_data *index_data;
    // Memory for output buffering. If output buffering is unneeded,
    // it contains NULL.
    char *out_buffer;
};

// Structure for parameters in worker threads
struct thread_data {
    npy_intp start;
    npy_intp vlen;
    npy_intp block_size;
    vm_params params;
    int ret_code;
    int *pc_error;
    char **errmsg;
    // NOTE: memsteps, iter, and reduce_iter are arrays, they MUST be allocated
    // to length `global_max_threads` before module load.
    // One memsteps array per thread
    // npy_intp *memsteps[MAX_THREADS];
    npy_intp **memsteps;
    // One iterator per thread */
    // NpyIter *iter[MAX_THREADS];
    NpyIter **iter;
    // When doing nested iteration for a reduction
    // NpyIter *reduce_iter[MAX_THREADS]
    NpyIter **reduce_iter;
    // Flag indicating reduction is the outer loop instead of the inner
    bool reduction_outer_loop;
    // Flag indicating whether output buffering is needed
    bool need_output_buffering;
};

// Global state which holds thread parameters
extern thread_data th_params;

PyObject *NumExpr_run(NumExprObject *self, PyObject *args, PyObject *kwds);

char get_return_sig(PyObject* program);
int check_program(NumExprObject *self);
int get_temps_space(const vm_params& params, char **mem, size_t block_size);
void free_temps_space(const vm_params& params, char **mem);
int vm_engine_iter_task(NpyIter *iter, npy_intp *memsteps,
                    const vm_params& params, int *pc_error, char **errmsg);

#endif // NUMEXPR_INTERPRETER_HPP
