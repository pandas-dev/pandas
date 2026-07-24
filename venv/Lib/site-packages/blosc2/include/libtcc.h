#ifndef LIBTCC_H
#define LIBTCC_H

#ifndef LIBTCCAPI
# define LIBTCCAPI
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************/
/* set custom allocator for all allocations (optional), NULL for default. */
typedef void *TCCReallocFunc(void *ptr, unsigned long size);
LIBTCCAPI void tcc_set_realloc(TCCReallocFunc *my_realloc);

/*****************************/
typedef struct TCCState TCCState;

/* create a new TCC compilation context */
LIBTCCAPI TCCState *tcc_new(void);

/* free a TCC compilation context */
LIBTCCAPI void tcc_delete(TCCState *s);

/* set CONFIG_TCCDIR at runtime */
LIBTCCAPI void tcc_set_lib_path(TCCState *s, const char *path);

/* set error/warning callback (optional) */
typedef void TCCErrorFunc(void *opaque, const char *msg);
LIBTCCAPI void tcc_set_error_func(TCCState *s, void *error_opaque, TCCErrorFunc *error_func);

/* set options as from command line (multiple supported) */
LIBTCCAPI int tcc_set_options(TCCState *s, const char *str);

/*****************************/
/* preprocessor */

/* add include path */
LIBTCCAPI int tcc_add_include_path(TCCState *s, const char *pathname);

/* add in system include path */
LIBTCCAPI int tcc_add_sysinclude_path(TCCState *s, const char *pathname);

/* define preprocessor symbol 'sym'. value can be NULL, sym can be "sym=val" */
LIBTCCAPI void tcc_define_symbol(TCCState *s, const char *sym, const char *value);

/* undefine preprocess symbol 'sym' */
LIBTCCAPI void tcc_undefine_symbol(TCCState *s, const char *sym);

/*****************************/
/* compiling */

/* add a file (C file, dll, object, library, ld script). Return -1 if error. */
LIBTCCAPI int tcc_add_file(TCCState *s, const char *filename);

/* compile a string containing a C source. Return -1 if error. */
LIBTCCAPI int tcc_compile_string(TCCState *s, const char *buf);

/* Tip: to have more specific errors/warnings from tcc_compile_string(),
   you can prefix the string with "#line <num> \"<filename>\"\n" */

/*****************************/
/* linking commands */

/* set output type. MUST BE CALLED before any compilation */
LIBTCCAPI int tcc_set_output_type(TCCState *s, int output_type);
#define TCC_OUTPUT_MEMORY   1 /* output will be run in memory */
#define TCC_OUTPUT_EXE      2 /* executable file */
#define TCC_OUTPUT_DLL      4 /* dynamic library */
#define TCC_OUTPUT_OBJ      3 /* object file */
#define TCC_OUTPUT_PREPROCESS 5 /* only preprocess */

/* equivalent to -Lpath option */
LIBTCCAPI int tcc_add_library_path(TCCState *s, const char *pathname);

/* the library name is the same as the argument of the '-l' option */
LIBTCCAPI int tcc_add_library(TCCState *s, const char *libraryname);

/* add a symbol to the compiled program */
LIBTCCAPI int tcc_add_symbol(TCCState *s, const char *name, const void *val);

/* output an executable, library or object file. DO NOT call
   tcc_relocate() before. */
LIBTCCAPI int tcc_output_file(TCCState *s, const char *filename);

/* link and run main() function and return its value. DO NOT call
   tcc_relocate() before. */
LIBTCCAPI int tcc_run(TCCState *s, int argc, char **argv);

/* do all relocations (needed before using tcc_get_symbol()) */
LIBTCCAPI int tcc_relocate(TCCState *s1);

/* return symbol value or NULL if not found */
LIBTCCAPI void *tcc_get_symbol(TCCState *s, const char *name);

/* list all (global) symbols and their values via 'symbol_cb()' */
LIBTCCAPI void tcc_list_symbols(TCCState *s, void *ctx,
    void (*symbol_cb)(void *ctx, const char *name, const void *val));

/* experimental/advanced section (see libtcc_test_mt.c for an example) */

/* catch runtime exceptions (optionally limit backtraces at top_func),
   when using tcc_set_options("-bt") and when not using tcc_run() */
LIBTCCAPI void *_tcc_setjmp(TCCState *s1, void *jmp_buf, void *top_func, void *longjmp);
#define tcc_setjmp(s1,jb,f) setjmp(_tcc_setjmp(s1, jb, f, longjmp))

/* debugging */
/* For debugging to work you have to enable it with tcc_set_options */

/* compile a string containing a C source. Return -1 if error.
   Write the string to file filename if debug is set. */
LIBTCCAPI int tcc_compile_string_file(TCCState *s, const char *buf, const char *filename);

/* Output object file. This must be done after tcc_relocate.
   It only generates the file if debug is set.
   The filename can be loaded with gdb command add-symbol-file */
LIBTCCAPI int elf_output_obj(TCCState *s1, const char *filename);

/* Set base address for wasm32 data/stack layout (default 1024).
   Call before tcc_output_file(). Only meaningful for TCC_TARGET_WASM32. */
LIBTCCAPI void tcc_set_wasm_data_base(TCCState *s, unsigned int base);

/* custom error printer for runtime exceptions. Returning 0 stops backtrace */
typedef int TCCBtFunc(void *udata, void *pc, const char *file, int line, const char* func, const char *msg);
LIBTCCAPI void tcc_set_backtrace_func(TCCState *s1, void* userdata, TCCBtFunc*);

#ifdef __cplusplus
}
#endif

#endif
