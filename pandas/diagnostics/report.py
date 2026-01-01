from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Any


def _fmt_bytes(n: int) -> str:
    # human-friendly bytes formatter
    n_f = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(n_f) < 1024.0:
            return f"{n_f:.2f}{unit}"
        n_f /= 1024.0
    return f"{n_f:.2f}EB"


@dataclass(frozen=True)
class PhaseRecord:
    name: str
    seconds: float
    mem_bytes_net: int | None = None


@dataclass(frozen=True)
class DiagnosticsReport:
    total_seconds: float
    phases: list[PhaseRecord]
    counters: dict[str, int]
    metadata: dict[str, Any]
    memory_bytes_net: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "total_seconds": self.total_seconds,
            "memory_bytes_net": self.memory_bytes_net,
            "metadata": self.metadata,
            "counters": self.counters,
            "phases": [
                {
                    "name": p.name,
                    "seconds": p.seconds,
                    "memory_bytes_net": p.mem_bytes_net,
                }
                for p in self.phases
            ],
        }

    def to_text(self) -> str:
        lines: list[str] = []
        op = self.metadata.get("operation", "operation")
        lines.append(f"pandas diagnostics report: {op}")
        lines.append(f"total: {self.total_seconds:.6f}s")

        if self.memory_bytes_net is not None:
            mem = int(self.memory_bytes_net)
            lines.append(f"net tracemalloc bytes: {mem} ({_fmt_bytes(mem)})")

        if self.metadata:
            keys = (
                "how",
                "left_rows",
                "right_rows",
                "left_cols",
                "right_cols",
                "copy_on_write",
            )
            meta_compact = {k: self.metadata.get(k) for k in keys if k in self.metadata}
            if meta_compact:
                lines.append(f"meta: {meta_compact}")

        if self.counters:
            lines.append("counters:")
            for k, v in sorted(self.counters.items()):
                if isinstance(v, int) and k.endswith(("_nbytes", "_bytes")):
                    lines.append(f"  - {k}: {v} ({_fmt_bytes(v)})")
                else:
                    lines.append(f"  - {k}: {v}")

        if self.phases:
            lines.append("phases:")
            phase_lines: list[str] = []
            for p in self.phases:
                if p.mem_bytes_net is None:
                    phase_lines.append(f"  - {p.name}: {p.seconds:.6f}s")
                else:
                    m = int(p.mem_bytes_net)
                    phase_lines.append(
                        f"  - {p.name}: {p.seconds:.6f}s, mem: {m} ({_fmt_bytes(m)})"
                    )
            lines.extend(phase_lines)

        return "\n".join(lines)


def attach_report(obj: Any, report: DiagnosticsReport) -> None:
    """
    Best-effort attach to result object.
    """
    try:
        setattr(obj, "_diagnostics_report", report)
    except Exception:
        return


def get_attached_report(obj: Any) -> DiagnosticsReport | None:
    rep = getattr(obj, "_diagnostics_report", None)
    return rep if isinstance(rep, DiagnosticsReport) else None


def explain_attached(obj: Any, *, format: str = "text") -> Any:
    rep = get_attached_report(obj)
    if rep is None:
        raise ValueError("No diagnostics report is attached to this object.")
    if format == "text":
        return rep.to_text()
    if format == "json":
        return rep.to_dict()
    if format == "json-string":
        return json.dumps(rep.to_dict(), indent=2, sort_keys=True)
    raise ValueError(
        f"Unknown format={format!r}. Expected 'text', 'json', or 'json-string'."
    )
