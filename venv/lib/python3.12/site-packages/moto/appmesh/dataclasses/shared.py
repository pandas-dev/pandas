from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, Dict, Literal, Optional
from uuid import uuid4

from moto.appmesh.utils.common import clean_dict

Status = Dict[Literal["status"], str]


@dataclass
class Metadata:
    arn: str
    mesh_owner: str
    resource_owner: str
    created_at: datetime = datetime.now()
    last_updated_at: datetime = datetime.now()
    uid: str = uuid4().hex
    version: int = 1

    def update_timestamp(self) -> None:
        self.last_updated_at = datetime.now()


@dataclass
class Duration:
    unit: str
    value: int
    to_dict = asdict


class MissingField:
    def to_dict(self) -> None:
        return


@dataclass
class Timeout:
    idle: Optional[Duration] = field(default=None)
    per_request: Optional[Duration] = field(default=None)

    def to_dict(self) -> Dict[str, Any]:  # type: ignore[misc]
        return clean_dict(
            {
                "idle": (self.idle or MissingField()).to_dict(),
                "perRequest": (self.per_request or MissingField()).to_dict(),
            }
        )
