from tensorflow import Tensor
from tensorflow._aliases import DTypeLike, TensorCompatible

def hamming_window(
    window_length: TensorCompatible, periodic: bool | TensorCompatible = True, dtype: DTypeLike = ..., name: str | None = None
) -> Tensor: ...
