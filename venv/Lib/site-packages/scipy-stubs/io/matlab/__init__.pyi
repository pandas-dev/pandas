# Deprecated namespaces, to be removed in v2.0.0
from . import (
    byteordercodes as byteordercodes,
    mio as mio,
    mio4 as mio4,
    mio5 as mio5,
    mio5_params as mio5_params,
    mio5_utils as mio5_utils,
    mio_utils as mio_utils,
    miobase as miobase,
    streams as streams,
)
from ._mio import loadmat, savemat, whosmat
from ._mio5 import varmats_from_mat
from ._mio5_params import MatlabFunction, MatlabObject, MatlabOpaque, mat_struct
from ._miobase import MatReadError, MatReadWarning, MatWriteError, MatWriteWarning, matfile_version

__all__ = [
    "MatReadError",
    "MatReadWarning",
    "MatWriteError",
    "MatWriteWarning",
    "MatlabFunction",
    "MatlabObject",
    "MatlabOpaque",
    "loadmat",
    "mat_struct",
    "matfile_version",
    "savemat",
    "varmats_from_mat",
    "whosmat",
]
