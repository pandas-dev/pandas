from typing import Any

from moto.core.common_models import BaseModel
from moto.core.utils import unix_time


class Timeout(BaseModel):
    def __init__(self, obj: Any, timestamp: float, kind: str):
        self.obj = obj
        self.timestamp = timestamp
        self.kind = kind

    @property
    def reached(self) -> bool:
        return unix_time() >= self.timestamp
