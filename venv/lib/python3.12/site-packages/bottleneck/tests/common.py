from hypothesis.extra.numpy import array_shapes
from hypothesis.extra.numpy import arrays as hy_arrays
from hypothesis.extra.numpy import floating_dtypes, integer_dtypes
from hypothesis.strategies import one_of

hy_array_gen = hy_arrays(
    dtype=one_of(integer_dtypes(sizes=(32, 64)), floating_dtypes(sizes=(32, 64))),
    shape=array_shapes(),
)

hy_int_array_gen = hy_arrays(
    dtype=integer_dtypes(sizes=(32, 64)),
    shape=array_shapes(),
)
