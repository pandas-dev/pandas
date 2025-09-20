#include "cuda_fp16.h"

#define FNDEF(fname) __numba_wrapper_ ## fname

#define UNARY_FUNCTION(fname) extern "C" __device__ int\
  FNDEF(fname)(						\
    short* return_value,\
  short x\
)\
{\
    __half retval = fname(__short_as_half (x));\
\
  *return_value = __half_as_short (retval);\
  /* Signal that no Python exception occurred */	\
  return 0;\
}\

extern "C" __device__ int
FNDEF(hdiv)(
  short* return_value,
  short x,
  short y
)
{
  __half retval = __hdiv(__short_as_half (x), __short_as_half (y));
  
  *return_value = __half_as_short (retval);
  // Signal that no Python exception occurred
  return 0;
}

UNARY_FUNCTION(hsin)
UNARY_FUNCTION(hcos)
UNARY_FUNCTION(hlog)
UNARY_FUNCTION(hlog10)
UNARY_FUNCTION(hlog2)
UNARY_FUNCTION(hexp)
UNARY_FUNCTION(hexp10)
UNARY_FUNCTION(hexp2)
UNARY_FUNCTION(hsqrt)
UNARY_FUNCTION(hrsqrt)
UNARY_FUNCTION(hfloor)
UNARY_FUNCTION(hceil)
UNARY_FUNCTION(hrcp)
UNARY_FUNCTION(hrint)
UNARY_FUNCTION(htrunc)

