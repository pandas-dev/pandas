from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Any


def _fmt_bytes(n: int) -> str:
    n_f = float(abs(n))
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if n_f < 1024.0:
            return f"{n_f:.2f}{unit}"
        n_f /= 1024.0
    return f"{n_f:.2f}EB"


def _fmt_bytes_signed(n: int) -> str:
    sign = "+" if n > 0 else ("-" if n < 0 else "")
    return f"{sign}{_fmt_bytes(n)}" if sign else _fmt_bytes(n)


def _fmt_seconds(s: float) -> str:
    if s < 1e-6:
        return f"{s * 1e9:.3f}ns"
    if s < 1e-3:
        return f"{s * 1e6:.3f}µs"
    if s < 1.0:
        return f"{s * 1e3:.3f}ms"
    return f"{s:.6f}s"


def _pct(part: float, whole: float) -> str:
    if whole <= 0:
        return ""
    return f"{(100.0 * part / whole):5.1f}%"


def _truncate_middle(s: str, max_len: int) -> str:
    if max_len <= 0:
        return ""
    if len(s) <= max_len:
        return s
    if max_len <= 3:
        return s[:max_len]
    k = (max_len - 3) // 2
    return s[:k] + "..." + s[-(max_len - 3 - k) :]


def _looks_numeric(s: str) -> bool:
    ss = s.strip()
    if not ss:
        return False
    if ss.endswith("%"):
        ss = ss[:-1].strip()
    # quick numeric-ish heuristic
    try:
        float(ss.replace(",", ""))
        return True
    except Exception:
        return False


def _table(
    headers: list[str],
    rows: list[list[str]],
    *,
    max_width: int = 120,
    right_align_cols: set[int] | None = None,
) -> str:
    if not rows:
        return ""

    cols = len(headers)
    widths = [len(h) for h in headers]
    for r in rows:
        for i in range(cols):
            widths[i] = max(widths[i], len(r[i]))

    sep = "  "
    total = sum(widths) + len(sep) * (cols - 1)

    if total > max_width:
        overflow = total - max_width
        for i in range(cols - 1, -1, -1):
            if overflow <= 0:
                break
            min_w = min(widths[i], 10)
            can_shrink = widths[i] - min_w
            if can_shrink <= 0:
                continue
            shrink = min(can_shrink, overflow)
            widths[i] -= shrink
            overflow -= shrink

    if right_align_cols is None:
        right_align_cols = set()

    def fmt_row(items: list[str]) -> str:
        out: list[str] = []
        for i, item in enumerate(items):
            cell = _truncate_middle(item, widths[i])
            right = (i in right_align_cols) or _looks_numeric(cell)
            out.append(cell.rjust(widths[i]) if right else cell.ljust(widths[i]))
        return sep.join(out).rstrip()

    hdr = fmt_row(headers)
    line = "-" * min(len(hdr), max_width)
    body = "\n".join(fmt_row(r) for r in rows)
    return f"{hdr}\n{line}\n{body}"


@dataclass(frozen=True)
class PhaseRecord:
    name: str
    seconds: float
    mem_bytes_net: int | None = None
    mem_current_bytes: int | None = None
    mem_peak_bytes: int | None = None


@dataclass(frozen=True)
class DiagnosticsReport:
    total_seconds: float

    phases: list[PhaseRecord]
    auto_phases: list[PhaseRecord]

    counters: dict[str, int]
    metadata: dict[str, Any]

    memory_bytes_net: int | None = None
    memory_current_bytes: int | None = None
    memory_peak_bytes: int | None = None

    profile_top: list[dict[str, Any]] | None = None
    allocations_top: list[dict[str, Any]] | None = None

    def filtered(self, *, include_diagnostics_phases: bool = True) -> DiagnosticsReport:
        """
        Return a copy of this report with optional filtering applied.
        Currently used to hide/show diag.* wrapper phases.
        """
        if include_diagnostics_phases:
            return self

        phases = [p for p in self.phases if not p.name.startswith("diag.")]
        return DiagnosticsReport(
            total_seconds=self.total_seconds,
            phases=phases,
            auto_phases=list(self.auto_phases),
            counters=dict(self.counters),
            metadata=dict(self.metadata),
            memory_bytes_net=self.memory_bytes_net,
            memory_current_bytes=self.memory_current_bytes,
            memory_peak_bytes=self.memory_peak_bytes,
            profile_top=list(self.profile_top)
            if self.profile_top is not None
            else None,
            allocations_top=list(self.allocations_top)
            if self.allocations_top is not None
            else None,
        )

    def to_dict(self) -> dict[str, Any]:
        def _phase(p: PhaseRecord) -> dict[str, Any]:
            return {
                "name": p.name,
                "seconds": p.seconds,
                "memory_bytes_net": p.mem_bytes_net,
                "memory_current_bytes": p.mem_current_bytes,
                "memory_peak_bytes": p.mem_peak_bytes,
            }

        return {
            "total_seconds": self.total_seconds,
            "memory_bytes_net": self.memory_bytes_net,
            "memory_current_bytes": self.memory_current_bytes,
            "memory_peak_bytes": self.memory_peak_bytes,
            "metadata": self.metadata,
            "counters": self.counters,
            "phases": [_phase(p) for p in self.phases],
            "auto_phases": [_phase(p) for p in self.auto_phases],
            "profile_top": self.profile_top or [],
            "allocations_top": self.allocations_top or [],
        }

    def to_json_string(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    def to_text(self) -> str:
        op = self.metadata.get("operation", "operation")
        callable_name = self.metadata.get("callable", None)

        call_s = None
        finalize_s = None
        for p in self.phases:
            if p.name == "call":
                call_s = p.seconds
            elif p.name == "finalize":
                finalize_s = p.seconds

        lines: list[str] = []
        lines.append(f"pandas diagnostics report: {op}")
        if callable_name:
            lines.append(f"callable: {callable_name}")

        summary = f"total: {_fmt_seconds(self.total_seconds)}"
        if call_s is not None:
            summary += f" | call: {_fmt_seconds(call_s)}"
        if finalize_s is not None:
            summary += f" | diagnostics overhead (finalize): {_fmt_seconds(finalize_s)}"
        lines.append(summary)

        if (
            self.memory_bytes_net is not None
            or self.memory_current_bytes is not None
            or self.memory_peak_bytes is not None
        ):
            parts = []
            if self.memory_bytes_net is not None:
                parts.append(f"netΔ {_fmt_bytes_signed(int(self.memory_bytes_net))}")
            if self.memory_current_bytes is not None:
                parts.append(f"current {_fmt_bytes(int(self.memory_current_bytes))}")
            if self.memory_peak_bytes is not None:
                parts.append(f"peak {_fmt_bytes(int(self.memory_peak_bytes))}")
            lines.append("memory: " + " | ".join(parts))

        if self.metadata:
            prefer = (
                "copy_on_write",
                "left_rows",
                "left_cols",
                "right_rows",
                "right_cols",
                "how",
                "on",
                "left_on",
                "right_on",
                "left_index",
                "right_index",
                "sort",
                "validate",
                "n_objs",
                "concat_axis",
                "concat_join",
                "concat_ignore_index",
                "sort_by",
                "sort_axis",
                "sort_ascending",
                "agg_spec",
                "csv_file",
                "csv_nrows",
                "csv_engine",
            )
            kv = [f"{k}={self.metadata.get(k)!r}" for k in prefer if k in self.metadata]
            if kv:
                lines.append("meta: " + ", ".join(kv))

        if self.counters:
            rows: list[list[str]] = []
            for k, v in sorted(self.counters.items()):
                if isinstance(v, int) and k.endswith(("_nbytes", "_bytes")):
                    rows.append([k, f"{v} ({_fmt_bytes(v)})"])
                else:
                    rows.append([k, str(v)])
            lines.append("")
            lines.append("counters:")
            lines.append(
                _table(
                    ["name", "value"],
                    rows,
                    max_width=120,
                    right_align_cols={1},
                )
            )

        # -----------------
        # Sequential phases
        # -----------------
        if self.phases:
            call_total = (
                float(call_s)
                if call_s is not None
                else float(self.total_seconds or 0.0)
            )

            any_cur = any(p.mem_current_bytes is not None for p in self.phases)
            any_peak = any(p.mem_peak_bytes is not None for p in self.phases)
            any_net = any(p.mem_bytes_net is not None for p in self.phases)

            headers = ["phase", "time", "%call", "cum"]
            show_delta = any_net or any_cur
            if show_delta:
                headers.append("memΔ")
            if any_cur:
                headers.append("mem_cur")
            if any_peak:
                headers.append("mem_peak")

            rows: list[list[str]] = []
            cum = 0.0
            prev_cur: int | None = None

            for p in self.phases:
                cum += float(p.seconds)

                delta_s: str = ""
                if show_delta:
                    if p.mem_bytes_net is not None:
                        delta_s = _fmt_bytes_signed(int(p.mem_bytes_net))
                    elif p.mem_current_bytes is not None and prev_cur is not None:
                        delta_s = _fmt_bytes_signed(
                            int(p.mem_current_bytes) - int(prev_cur)
                        )
                    else:
                        delta_s = ""

                if p.mem_current_bytes is not None:
                    prev_cur = int(p.mem_current_bytes)

                row = [
                    p.name,
                    _fmt_seconds(float(p.seconds)),
                    _pct(float(p.seconds), call_total),
                    _fmt_seconds(cum),
                ]
                if show_delta:
                    row.append(delta_s)
                if any_cur:
                    row.append(
                        _fmt_bytes(int(p.mem_current_bytes))
                        if p.mem_current_bytes is not None
                        else ""
                    )
                if any_peak:
                    row.append(
                        _fmt_bytes(int(p.mem_peak_bytes))
                        if p.mem_peak_bytes is not None
                        else ""
                    )

                rows.append(row)

            right_cols = set(range(1, len(headers)))

            lines.append("")
            lines.append("phases (sequential):")
            lines.append(
                _table(headers, rows, max_width=140, right_align_cols=right_cols)
            )

        # -----------------
        # Auto phases (bucketed attribution)
        # -----------------
        if self.auto_phases:
            ap = [
                p
                for p in self.auto_phases
                if float(p.seconds) > 0.0 or (p.mem_bytes_net not in (None, 0))
            ]
            if ap:
                call_total = (
                    float(call_s)
                    if call_s is not None
                    else float(self.total_seconds or 0.0)
                )
                rows = []
                for p in ap:
                    memd = (
                        _fmt_bytes_signed(int(p.mem_bytes_net))
                        if p.mem_bytes_net is not None
                        else ""
                    )
                    rows.append(
                        [
                            p.name,
                            _fmt_seconds(float(p.seconds)),
                            _pct(float(p.seconds), call_total),
                            memd,
                        ]
                    )

                lines.append("")
                lines.append("phases (auto buckets):")
                lines.append(
                    _table(
                        ["bucket", "time", "%call", "memΔ"],
                        rows,
                        max_width=120,
                        right_align_cols={1, 2, 3},
                    )
                )

        # -----------------
        # Hotspots
        # -----------------
        if self.profile_top:
            call_total = (
                float(call_s)
                if call_s is not None
                else float(self.total_seconds or 0.0)
            )
            rows = []
            for h in (self.profile_top or [])[:50]:
                cs = float(h.get("cum_seconds", 0.0))
                calls = int(h.get("calls", 0))
                fn = str(h.get("function", ""))
                file = str(h.get("filename", ""))
                line = h.get("line", None)
                loc = f"{file}:{line}" if line is not None else file
                rows.append(
                    [_fmt_seconds(cs), _pct(cs, call_total), str(calls), fn, loc]
                )

            if rows:
                lines.append("")
                lines.append("top pandas hotspots (profile):")
                lines.append(
                    _table(
                        ["cum", "%call", "calls", "function", "location"],
                        rows,
                        max_width=160,
                        right_align_cols={0, 1, 2},
                    )
                )

        # -----------------
        # Allocation diffs
        # -----------------
        if self.allocations_top:
            rows = []
            for a in (self.allocations_top or [])[:50]:
                b = int(a.get("bytes_net", 0))
                c = int(a.get("count_net", 0))
                f = str(a.get("filename", ""))
                rows.append([_fmt_bytes_signed(b), str(b), str(c), f])

            if rows:
                lines.append("")
                lines.append("top file alloc diffs (tracemalloc):")
                lines.append(
                    _table(
                        ["pretty", "bytesΔ", "countΔ", "filename"],
                        rows,
                        max_width=160,
                        right_align_cols={0, 1, 2},
                    )
                )

        return "\n".join(lines)


__all__ = ["DiagnosticsReport", "PhaseRecord"]
