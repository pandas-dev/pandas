"""
Useful options to debug LLVM passes

llvm.set_option("test", "-debug-pass=Details")
llvm.set_option("test", "-debug-pass=Executions")
llvm.set_option("test", "-debug-pass=Arguments")
llvm.set_option("test", "-debug-pass=Structure")
llvm.set_option("test", "-debug-only=loop-vectorize")
llvm.set_option("test", "-help-hidden")

"""

from llvmlite import binding as llvm


def _inlining_threshold(optlevel, sizelevel=0):
    """
    Compute the inlining threshold for the desired optimisation level

    Refer to http://llvm.org/docs/doxygen/html/InlineSimple_8cpp_source.html
    """
    if optlevel > 2:
        return 275

    # -Os
    if sizelevel == 1:
        return 75

    # -Oz
    if sizelevel == 2:
        return 25

    return 225


def create_pass_manager_builder(opt=2, loop_vectorize=False,
                                slp_vectorize=False):
    """
    Create an LLVM pass manager with the desired optimisation level and options.
    """
    pmb = llvm.create_pass_manager_builder()
    pmb.opt_level = opt
    pmb.loop_vectorize = loop_vectorize
    pmb.slp_vectorize = slp_vectorize
    pmb.inlining_threshold = _inlining_threshold(opt)
    return pmb
