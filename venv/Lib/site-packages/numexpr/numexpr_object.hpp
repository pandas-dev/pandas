#ifndef NUMEXPR_OBJECT_HPP
#define NUMEXPR_OBJECT_HPP
/*********************************************************************
  Numexpr - Fast numerical array expression evaluator for NumPy.

      License: MIT
      Author:  See AUTHORS.txt

  See LICENSE.txt for details about copyright and rights to use.
**********************************************************************/

struct NumExprObject
{
    PyObject_HEAD
    PyObject *signature;    /* a python string */
    PyObject *tempsig;
    PyObject *constsig;
    PyObject *fullsig;
    PyObject *program;      /* a python string */
    PyObject *constants;    /* a tuple of int/float/complex */
    PyObject *input_names;  /* tuple of strings */
    char **mem;             /* pointers to registers */
    char *rawmem;           /* a chunks of raw memory for storing registers */
    npy_intp *memsteps;
    npy_intp *memsizes;
    int  rawmemsize;
    int  n_inputs;
    int  n_constants;
    int  n_temps;
};

extern PyTypeObject NumExprType;

#endif // NUMEXPR_OBJECT_HPP
