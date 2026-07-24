from typing import final

from scipy._lib._uarray._backend import _Backend

@final
class NumPyBackend(_Backend): ...

@final
class EchoBackend(_Backend[None]): ...
