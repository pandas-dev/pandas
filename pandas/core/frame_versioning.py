# pandas/core/frame_versioning.py
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import uuid
from typing import Dict, Optional

import pandas as pd


def _generate_snapshot_id(name: Optional[str] = None) -> str:
    if name:
        return name
    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%S%fZ")
    uid = uuid.uuid4().hex[:8]
    return f"{ts}-{uid}"


@dataclass
class SnapshotMeta:
    name: str
    created_at: datetime


class DataFrameSnapshotStore:
    """
    Per-DataFrame snapshot store.
    Stores deep copies of DataFrames (safe, simple).
    """

    def __init__(self) -> None:
        # snapshot_id -> DataFrame
        self._snapshots: Dict[str, pd.DataFrame] = {}
        self._meta: Dict[str, SnapshotMeta] = {}

    def snapshot(self, df: pd.DataFrame, name: Optional[str] = None) -> str:
        sid = _generate_snapshot_id(name)
        # deep copy for safety
        self._snapshots[sid] = df.copy(deep=True)
        self._meta[sid] = SnapshotMeta(name=sid, created_at=datetime.utcnow())
        return sid

    def restore(self, name: str) -> pd.DataFrame:
        if name not in self._snapshots:
            raise KeyError(f"Snapshot not found: {name}")
        # return a deep copy so modifications don't change stored snapshot
        return self._snapshots[name].copy(deep=True)

    def list(self) -> list[str]:
        return list(self._snapshots.keys())

    def drop(self, name: str) -> None:
        if name not in self._snapshots:
            raise KeyError(f"Snapshot not found: {name}")
        del self._snapshots[name]
        del self._meta[name]

    def clear(self) -> None:
        self._snapshots.clear()
        self._meta.clear()

    def info(self, name: Optional[str] = None) -> dict:
        if name:
            if name not in self._meta:
                raise KeyError(f"Snapshot not found: {name}")
            meta = self._meta[name]
            return {"name": meta.name, "created_at": meta.created_at.isoformat()}
        return {"count": len(self._snapshots), "snapshots": [m.name for m in self._meta.values()]}
