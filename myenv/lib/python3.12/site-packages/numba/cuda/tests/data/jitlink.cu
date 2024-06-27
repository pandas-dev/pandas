// Compile with:
//
//   nvcc -gencode arch=compute_50,code=compute_50 -rdc true -ptx jitlink.cu
//
// using the oldest supported toolkit version (10.2 at the time of writing).

extern "C" __device__
int bar(int *out, int a)
{
  *out = a * 2;
  return 0;
}


// The out argument is necessary due to Numba's CUDA calling convention, which
// always reserves the first parameter for a pointer to a returned value, even
// if there is no return value.
extern "C" __device__
int array_mutator(void *out, int *a)
{
  a[0] = a[1];
  return 0;
}  
