from mypyc.ir.deps import LIBRT_RANDOM
from mypyc.ir.ops import ERR_MAGIC, ERR_NEVER
from mypyc.ir.rtypes import float_rprimitive, int64_rprimitive, random_rprimitive
from mypyc.primitives.registry import function_op, method_op

# Random() -- construct with OS entropy
function_op(
    name="librt.random.Random",
    arg_types=[],
    return_type=random_rprimitive,
    c_function_name="LibRTRandom_Random_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Random(seed) -- construct with integer seed
function_op(
    name="librt.random.Random",
    arg_types=[int64_rprimitive],
    return_type=random_rprimitive,
    c_function_name="LibRTRandom_Random_from_seed_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Random.randint(a, b) -- return random integer in [a, b]
method_op(
    name="randint",
    arg_types=[random_rprimitive, int64_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="LibRTRandom_Random_randint_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Random.randrange(stop) -- return random integer in [0, stop)
method_op(
    name="randrange",
    arg_types=[random_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="LibRTRandom_Random_randrange1_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Random.randrange(start, stop) -- return random integer in [start, stop)
method_op(
    name="randrange",
    arg_types=[random_rprimitive, int64_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="LibRTRandom_Random_randrange2_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Random.random() -- return random float in [0.0, 1.0)
method_op(
    name="random",
    arg_types=[random_rprimitive],
    return_type=float_rprimitive,
    c_function_name="LibRTRandom_Random_random_internal",
    error_kind=ERR_NEVER,
    dependencies=[LIBRT_RANDOM],
)

# Module-level random() -- return random float using thread-local RNG
function_op(
    name="librt.random.random",
    arg_types=[],
    return_type=float_rprimitive,
    c_function_name="LibRTRandom_module_random_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Module-level randrange(stop) -- return random integer using thread-local RNG
function_op(
    name="librt.random.randrange",
    arg_types=[int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="LibRTRandom_module_randrange1_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Module-level randrange(start, stop) -- return random integer using thread-local RNG
function_op(
    name="librt.random.randrange",
    arg_types=[int64_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="LibRTRandom_module_randrange2_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)

# Module-level randint(a, b) -- return random integer using thread-local RNG
function_op(
    name="librt.random.randint",
    arg_types=[int64_rprimitive, int64_rprimitive],
    return_type=int64_rprimitive,
    c_function_name="LibRTRandom_module_randint_internal",
    error_kind=ERR_MAGIC,
    dependencies=[LIBRT_RANDOM],
)
