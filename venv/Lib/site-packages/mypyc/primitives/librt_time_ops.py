from mypyc.ir.deps import LIBRT_TIME
from mypyc.ir.ops import ERR_MAGIC_OVERLAPPING
from mypyc.ir.rtypes import float_rprimitive
from mypyc.primitives.registry import function_op

function_op(
    name="librt.time.time",
    arg_types=[],
    return_type=float_rprimitive,
    c_function_name="LibRTTime_time",
    error_kind=ERR_MAGIC_OVERLAPPING,
    dependencies=[LIBRT_TIME],
)
