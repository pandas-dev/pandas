"""
A collection of analysis utilities
"""

from ctypes import POINTER, c_char_p, c_int

from llvmlite.binding import ffi
from llvmlite.binding.module import parse_assembly


def get_function_cfg(func, show_inst=True):
    """Return a string of the control-flow graph of the function in DOT
    format. If the input `func` is not a materialized function, the module
    containing the function is parsed to create an actual LLVM module.
    The `show_inst` flag controls whether the instructions of each block
    are printed.
    """
    assert func is not None
    from llvmlite import ir
    if isinstance(func, ir.Function):
        mod = parse_assembly(str(func.module))
        func = mod.get_function(func.name)

    # Assume func is a materialized function
    with ffi.OutputString() as dotstr:
        ffi.lib.LLVMPY_WriteCFG(func, dotstr, show_inst)
        return str(dotstr)


def view_dot_graph(graph, filename=None, view=False):
    """
    View the given DOT source.  If view is True, the image is rendered
    and viewed by the default application in the system.  The file path of
    the output is returned.  If view is False, a graphviz.Source object is
    returned.  If view is False and the environment is in a IPython session,
    an IPython image object is returned and can be displayed inline in the
    notebook.

    This function requires the graphviz package.

    Args
    ----
    - graph [str]: a DOT source code
    - filename [str]: optional.  if given and view is True, this specifies
                      the file path for the rendered output to write to.
    - view [bool]: if True, opens the rendered output file.

    """
    # Optionally depends on graphviz package
    import graphviz as gv

    src = gv.Source(graph)
    if view:
        # Returns the output file path
        return src.render(filename, view=view)
    else:
        # Attempts to show the graph in IPython notebook
        try:
            __IPYTHON__
        except NameError:
            return src
        else:
            import IPython.display as display
            format = 'svg'
            return display.SVG(data=src.pipe(format))


# Ctypes binding
ffi.lib.LLVMPY_WriteCFG.argtypes = [ffi.LLVMValueRef, POINTER(c_char_p), c_int]
