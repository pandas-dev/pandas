import numpy as np

def hash_object_array(
    arr: np.ndarray,  # np.ndarray[object]
    key: str,
    encoding: str = "utf8",
) -> np.ndarray: ...  # np.ndarray[np.uint64]
