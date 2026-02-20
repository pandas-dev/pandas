from typing import Final, Literal

import numpy as np

_lqmn: Final[np.ufunc] = ...  # undocumented
_lqn: Final[np.ufunc] = ...  # undocumented
_rctj: Final[np.ufunc] = ...  # undocumented
_rcty: Final[np.ufunc] = ...  # undocumented

assoc_legendre_p_all: Final[dict[tuple[bool, Literal[0, 1, 2]], np.ufunc]] = ...
legendre_p_all: Final[tuple[np.ufunc, np.ufunc, np.ufunc]] = ...
sph_harm_y_all: Final[tuple[np.ufunc, np.ufunc, np.ufunc]] = ...
sph_legendre_p_all: Final[tuple[np.ufunc, np.ufunc, np.ufunc]] = ...

def _set_action(code: int, action: int, /) -> None: ...  # undocumented
