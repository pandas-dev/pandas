from _typeshed import Incomplete

from tensorflow._aliases import IntArray, IntDataSequence

Layout = Incomplete

class Mesh:
    def __init__(
        self,
        dim_names: list[str],
        global_device_ids: IntArray | IntDataSequence,
        local_device_ids: list[int],
        local_devices: list[Incomplete | str],
        mesh_name: str = "",
        global_devices: list[Incomplete | str] | None = None,
        use_xla_spmd: bool = False,
    ) -> None: ...

def __getattr__(name: str): ...  # incomplete module
