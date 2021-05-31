from typing import (
    Dict,
    Optional,
)

import numpy as np

from pandas.compat._optional import import_optional_dependency

from pandas.core.util.numba_ import (
    NUMBA_FUNC_CACHE,
    get_jit_arguments,
)


def generate_online_numba_ewma_func(engine_kwargs: Optional[Dict[str, bool]]):
    """
    Generate a numba jitted groupby ewma function specified by values
    from engine_kwargs.
    Parameters
    ----------
    engine_kwargs : dict
        dictionary of arguments to be passed into numba.jit
    Returns
    -------
    Numba function
    """
    nopython, nogil, parallel = get_jit_arguments(engine_kwargs)

    cache_key = (lambda x: x, "online_ewma")
    if cache_key in NUMBA_FUNC_CACHE:
        return NUMBA_FUNC_CACHE[cache_key]

    numba = import_optional_dependency("numba")

    @numba.jit(nopython=nopython, nogil=nogil, parallel=parallel)
    def online_ewma(
        values: np.ndarray,
        deltas: np.ndarray,
        minimum_periods: int,
        old_wt_factor: float,
        new_wt: float,
        old_wt: float,
        adjust: bool,
        ignore_na: bool,
    ):
        result = np.empty(len(values))

        weighted_avg = values[0]
        nobs = int(not np.isnan(weighted_avg))
        result[0] = weighted_avg if nobs >= minimum_periods else np.nan

        for j in range(1, len(values)):
            cur = values[j]
            is_observation = not np.isnan(cur)
            nobs += is_observation
            if not np.isnan(weighted_avg):

                if is_observation or not ignore_na:

                    # note that len(deltas) = len(vals) - 1 and deltas[i] is to be
                    # used in conjunction with vals[i+1]
                    old_wt *= old_wt_factor ** deltas[j - 1]
                    if is_observation:

                        # avoid numerical errors on constant series
                        if weighted_avg != cur:
                            weighted_avg = (
                                (old_wt * weighted_avg) + (new_wt * cur)
                            ) / (old_wt + new_wt)
                        if adjust:
                            old_wt += new_wt
                        else:
                            old_wt = 1.0
            elif is_observation:
                weighted_avg = cur

            result[j] = weighted_avg if nobs >= minimum_periods else np.nan

        return result, old_wt

    return online_ewma


class EWMeanState:
    def __init__(self, com, adjust, ignore_na):
        alpha = 1.0 / (1.0 + com)
        self.old_wt_factor = 1.0 - alpha
        self.new_wt = 1.0 if adjust else alpha
        self.old_wt = 1.0
        self.adjust = adjust
        self.ignore_na = ignore_na
        self.last_ewm = None

    def run_ewm(self, weighted_avg, deltas, min_periods, ewm_func):
        result, old_wt = ewm_func(
            weighted_avg,
            deltas,
            min_periods,
            self.old_wt_factor,
            self.new_wt,
            self.old_wt,
            self.adjust,
            self.ignore_na,
        )
        self.old_wt = old_wt
        self.last_ewm = result[-1]
        return result

    def reset(self):
        self.old_wt = 1
        self.last_ewm = None
