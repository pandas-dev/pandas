import numpy as np
import optype.numpy as onp

def linear_sum_assignment(
    cost_matrix: onp.ToFloat2D, maximize: onp.ToBool = False
) -> tuple[onp.Array1D[np.intp], onp.Array1D[np.intp]]: ...
