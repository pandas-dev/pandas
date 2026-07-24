from _typeshed import Incomplete

from ._monitor import TMonitor as TMonitor, TqdmSynchronisationWarning as TqdmSynchronisationWarning
from ._tqdm_pandas import tqdm_pandas as tqdm_pandas
from .cli import main as main
from .gui import tqdm as tqdm_gui, trange as tgrange
from .notebook import tqdm_notebook as tqdm_notebook_cls
from .std import (
    TqdmDeprecationWarning as TqdmDeprecationWarning,
    TqdmExperimentalWarning as TqdmExperimentalWarning,
    TqdmKeyError as TqdmKeyError,
    TqdmMonitorWarning as TqdmMonitorWarning,
    TqdmTypeError as TqdmTypeError,
    TqdmWarning as TqdmWarning,
    tqdm as tqdm,
    trange as trange,
)
from .version import __version__ as __version__

__all__ = [
    "tqdm",
    "tqdm_gui",
    "trange",
    "tgrange",
    "tqdm_pandas",
    "tqdm_notebook",
    "tnrange",
    "main",
    "TMonitor",
    "TqdmTypeError",
    "TqdmKeyError",
    "TqdmWarning",
    "TqdmDeprecationWarning",
    "TqdmExperimentalWarning",
    "TqdmMonitorWarning",
    "TqdmSynchronisationWarning",
    "__version__",
]

def tqdm_notebook(*args, **kwargs) -> tqdm_notebook_cls[Incomplete]: ...
def tnrange(*args, **kwargs) -> tqdm_notebook_cls[int]: ...
