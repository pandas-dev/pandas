import numpy as np

def calculate_variable_window_bounds(
    num_values: int,     # int64_t
    window_size: int,    # int64_t
    min_periods: object,
    center: object,
    closed: object,
    index: np.ndarray,  # const int64_t[:]
) -> tuple[
    np.ndarray,  # np.ndarray[np.int64]
    np.ndarray,  # np.ndarray[np.int64]
]: ...
