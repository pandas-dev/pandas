# defined in scipy/ndimage/src/_ni_label.pyx

from typing import TypeAlias

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

# matches `ctypedef fused data_t`
_data_t: TypeAlias = npc.integer | np.float32 | np.float64  # noqa: PYI042

###

class NeedMoreBits(Exception): ...  # undocumented

def get_nonzero_line(a: onp.ArrayND[_data_t]) -> int: ...  # undocumented
def get_read_line(a: onp.ArrayND[_data_t]) -> int: ...  # undocumented
def get_write_line(a: onp.ArrayND[_data_t]) -> int: ...  # undocumented

#
def _label(
    input: onp.ArrayND[_data_t | np.bool_],
    structure: onp.ArrayND[np.bool_ | npc.integer],
    output: onp.ArrayND[_data_t | np.bool_],
) -> int: ...  # undocumented
