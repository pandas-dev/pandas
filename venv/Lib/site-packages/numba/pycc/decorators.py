import re
import warnings

from numba.core import typing, sigutils
from numba.pycc.compiler import ExportEntry

# Registry is okay to be a global because we are using pycc as a standalone
# commandline tool.
export_registry = []


def export(prototype):
    warnings.warn("export() is deprecated, use the numba.pycc.CC API instead",
                  DeprecationWarning, stacklevel=2)

    sym, sig = parse_prototype(prototype)

    def wrappped(func):
        fn_argtys, fn_retty = sigutils.normalize_signature(sig)
        signature = typing.signature(fn_retty, *fn_argtys)
        entry = ExportEntry(symbol=sym, signature=signature, function=func)
        export_registry.append(entry)

    return wrappped


def exportmany(prototypes):
    warnings.warn("exportmany() is deprecated, use the numba.pycc.CC API instead",
                  DeprecationWarning, stacklevel=2)

    def wrapped(func):
        for proto in prototypes:
            export(proto)(func)
    return wrapped


def process_input_files(inputs):
    """
    Read input source files for execution of legacy @export / @exportmany
    decorators.
    """
    for ifile in inputs:
        with open(ifile) as fin:
            exec(compile(fin.read(), ifile, 'exec'))


def clear_export_registry():
    export_registry[:] = []


# --------------------------------- Internal ---------------------------------

re_symbol = re.compile(r'[_a-z][_a-z0-9]*', re.I)


def parse_prototype(text):
    """Separate the symbol and function-type in a a string with
    "symbol function-type" (e.g. "mult float(float, float)")

    Returns
    ---------
    (symbol_string, functype_string)
    """
    m = re_symbol.match(text)
    if not m:
        raise ValueError("Invalid function name for export prototype")
    s = m.start(0)
    e = m.end(0)
    symbol = text[s:e]
    functype = text[e + 1:]
    return symbol, functype

