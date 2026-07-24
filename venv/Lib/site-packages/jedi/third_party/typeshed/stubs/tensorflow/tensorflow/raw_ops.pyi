from typing import Literal

from tensorflow import Operation, Tensor
from tensorflow._aliases import TensorCompatible

def Fingerprint(data: TensorCompatible, method: Literal["farmhash64"], name: str | None = None) -> Tensor: ...
def Snapshot(input: TensorCompatible, name: str | None = None) -> Tensor: ...
def ResourceApplyAdagradV2(
    var: Tensor,
    accum: Tensor,
    lr: TensorCompatible,
    epsilon: TensorCompatible,
    grad: TensorCompatible,
    use_locking: bool = False,
    update_slots: bool = True,
    name: str | None = None,
) -> Operation: ...
def ResourceSparseApplyAdagradV2(
    var: Tensor,
    accum: Tensor,
    lr: TensorCompatible,
    epsilon: TensorCompatible,
    grad: TensorCompatible,
    indices: TensorCompatible,
    use_locking: bool = False,
    update_slots: bool = True,
    name: str | None = None,
) -> Operation: ...
def ResourceApplyAdam(
    var: Tensor,
    m: Tensor,
    v: Tensor,
    beta1_power: TensorCompatible,
    beta2_power: TensorCompatible,
    lr: TensorCompatible,
    beta1: TensorCompatible,
    beta2: TensorCompatible,
    epsilon: TensorCompatible,
    grad: TensorCompatible,
    use_locking: bool = False,
    use_nesterov: bool = False,
    name: str | None = None,
) -> Operation: ...
def __getattr__(name: str): ...  # incomplete module
