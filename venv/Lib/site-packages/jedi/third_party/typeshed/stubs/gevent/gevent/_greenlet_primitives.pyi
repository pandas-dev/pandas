from abc import abstractmethod
from typing import Any, NoReturn
from typing_extensions import disjoint_base

from gevent._types import _Loop
from greenlet import greenlet

class TrackedRawGreenlet(greenlet): ...

@disjoint_base
class SwitchOutGreenletWithLoop(TrackedRawGreenlet):
    @property
    @abstractmethod
    def loop(self) -> _Loop: ...
    @loop.setter
    def loop(self, value: _Loop) -> None: ...
    def switch(self) -> Any: ...
    def switch_out(self) -> NoReturn: ...

__all__ = ["TrackedRawGreenlet", "SwitchOutGreenletWithLoop"]
