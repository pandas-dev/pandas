# ruff: noqa: E402, PLC0415

_NPTC_BUILD = 2025_08_18


def _check_numpy_typing_compat() -> None:
    try:
        import numpy_typing_compat as nptc
    except ImportError:
        oh_no_an_import_error = (
            "`optype.numpy` requires `numpy-typing-compat` to be installed, which is "
            "available through the `optype[numpy]` extra. Please install it with "
            "`uv install optype[numpy]`."
        )
        raise ImportError(oh_no_an_import_error) from None

    if not nptc._check_version():  # noqa: SLF001
        from importlib.metadata import requires

        import numpy as np

        np_reqs = requires("numpy-typing-compat")
        assert np_reqs
        assert np_reqs[0].startswith("numpy")
        np_req = np_reqs[0]

        oh_no_an_import_error = (
            f"The installed version of 'numpy-typing-compat' ({nptc.__version__}) "
            f"requires {np_req!r}, but 'numpy=={np.__version__}' is installed. "
            f"Please install the compatible version of `numpy-typing-compat`, e.g. "
            f"with `uv pip install optype[numpy]` or `pixi add optype-numpy`."
        )
        raise ImportError(oh_no_an_import_error)

    v_build, v_major, v_minor = map(int, nptc.__version__.split(".", 2)[:3])

    if v_build < _NPTC_BUILD:
        v_req = f"{_NPTC_BUILD}.{v_major}.{v_minor}"

        oh_no_an_import_error = (
            f"The installed version of 'optype-numpy-compat' ({nptc.__version__}) is "
            f"unsupported. Please upgrade to {v_req}."
        )
        raise ImportError(oh_no_an_import_error)


_check_numpy_typing_compat()

# ruff: noqa: F403
from . import (
    _any_array,
    _any_dtype,
    _array,
    _dtype,
    _is,
    _literals,
    _scalar,
    _sequence_nd,
    _shape,
    _to,
    _ufunc,
    compat,
    ctypeslib,
    random,
)
from ._any_array import *
from ._any_dtype import *
from ._array import *
from ._dtype import *
from ._is import *
from ._literals import *
from ._scalar import *
from ._sequence_nd import *
from ._shape import *
from ._to import *
from ._ufunc import *

__all__ = ["compat", "ctypeslib", "random"]
__all__ += _literals.__all__
__all__ += _array.__all__
__all__ += _any_array.__all__
__all__ += _sequence_nd.__all__
__all__ += _dtype.__all__
__all__ += _any_dtype.__all__
__all__ += _shape.__all__
__all__ += _is.__all__
__all__ += _to.__all__
__all__ += _ufunc.__all__
__all__ += _scalar.__all__


def __dir__() -> list[str]:
    return __all__
