"""
Implement a rewrite pass on a LLVM module to remove unnecessary 
refcount operations.
"""

from llvmlite.ir.transforms import CallVisitor

from numba.core import types


class _MarkNrtCallVisitor(CallVisitor):
    """
    A pass to mark all NRT_incref and NRT_decref.
    """
    def __init__(self):
        self.marked = set()

    def visit_Call(self, instr):
        if getattr(instr.callee, 'name', '') in _accepted_nrtfns:
            self.marked.add(instr)


def _rewrite_function(function):
    # Mark NRT usage
    markpass = _MarkNrtCallVisitor()
    markpass.visit_Function(function)
    # Remove NRT usage
    for bb in function.basic_blocks:
        for inst in list(bb.instructions):
            if inst in markpass.marked:
                bb.instructions.remove(inst)


_accepted_nrtfns = 'NRT_incref', 'NRT_decref'


def _legalize(module, dmm, fndesc):
    """
    Legalize the code in the module.
    Returns True if the module is legal for the rewrite pass that removes
    unnecessary refcounts.
    """

    def valid_output(ty):
        """
        Valid output are any type that does not need refcount
        """
        model = dmm[ty]
        return not model.contains_nrt_meminfo()

    def valid_input(ty):
        """
        Valid input are any type that does not need refcount except Array.
        """
        return valid_output(ty) or isinstance(ty, types.Array)


    # Ensure no reference to function marked as
    # "numba_args_may_always_need_nrt"
    try:
        nmd = module.get_named_metadata("numba_args_may_always_need_nrt")
    except KeyError:
        # Nothing marked
        pass
    else:
        # Has functions marked as "numba_args_may_always_need_nrt"
        if len(nmd.operands) > 0:
            # The pass is illegal for this compilation unit.
            return False

    # More legalization base on function type
    argtypes = fndesc.argtypes
    restype = fndesc.restype
    calltypes = fndesc.calltypes

    # Legalize function arguments
    for argty in argtypes:
        if not valid_input(argty):
            return False

    # Legalize function return
    if not valid_output(restype):
        return False

    # Legalize all called functions
    for callty in calltypes.values():
        if callty is not None and not valid_output(callty.return_type):
            return False

    # Ensure no allocation
    for fn in module.functions:
        if fn.name.startswith("NRT_"):
            if fn.name not in _accepted_nrtfns:
                return False

    return True


def remove_unnecessary_nrt_usage(function, context, fndesc):
    """
    Remove unnecessary NRT incref/decref in the given LLVM function.
    It uses highlevel type info to determine if the function does not need NRT.
    Such a function does not:

    - return array object(s);
    - take arguments that need refcounting except array;
    - call function(s) that return refcounted object.

    In effect, the function will not capture or create references that extend
    the lifetime of any refcounted objects beyond the lifetime of the function.

    The rewrite is performed in place.
    If rewrite has happened, this function returns True, otherwise, it returns False.
    """
    dmm = context.data_model_manager
    if _legalize(function.module, dmm, fndesc):
        _rewrite_function(function)
        return True
    else:
        return False
