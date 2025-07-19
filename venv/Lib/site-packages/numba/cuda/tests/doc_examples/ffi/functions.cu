// magictoken.ex_mul_f32_f32.begin
// Foreign function example: multiplication of a pair of floats

extern "C" __device__ int
mul_f32_f32(
  float* return_value,
  float x,
  float y)
{
  // Compute result and store in caller-provided slot
  *return_value = x * y;

  // Signal that no Python exception occurred
  return 0;
}
// magictoken.ex_mul_f32_f32.end


// magictoken.ex_sum_reduce_proto.begin
extern "C"
__device__ int
sum_reduce(
  float* return_value,
  float* array,
  int n
);
// magictoken.ex_sum_reduce_proto.end


// Performs a simple reduction on an array passed by pointer using the
// ffi.from_buffer() method. Implements the prototype above.
extern "C"
__device__ int
sum_reduce(
  float* return_value,
  float* array,
  int n
)
{
  double sum = 0.0;

  for (size_t i = 0; i < n; ++i) {
    sum += array[i];
  }

  *return_value = (float)sum;

  return 0;
}
