from __future__ import annotations

import os
from typing import Any


def _shape(obj: Any) -> tuple[int | None, int | None]:
    try:
        s = getattr(obj, "shape", None)
        if s is None or len(s) < 2:
            return (None, None)
        return (int(s[0]), int(s[1]))
    except Exception:
        return (None, None)


def _maybe_len(obj: Any) -> int | None:
    try:
        return len(obj)
    except Exception:
        return None


def _callable_name(func: Any) -> str:
    mod = getattr(func, "__module__", "") or ""
    qual = getattr(func, "__qualname__", None) or getattr(func, "__name__", "callable")
    return f"{mod}.{qual}" if mod else str(qual)


def _safe_basename(path: Any) -> str | None:
    try:
        if isinstance(path, str):
            return os.path.basename(path)
    except Exception:
        pass
    return None


def _safe_value(v: Any, *, limit: int = 200) -> str:
    try:
        if v is None:
            return "None"
        if isinstance(v, (bool, int, float)):
            return str(v)
        if isinstance(v, str):
            return v[:limit]
        if isinstance(v, (list, tuple, set)):
            return f"{type(v).__name__}(len={len(v)})"
        if isinstance(v, dict):
            return f"dict(len={len(v)})"
        if callable(v):
            return getattr(v, "__qualname__", getattr(v, "__name__", "callable"))
        return repr(v)[:limit]
    except Exception:
        return "<unrepr>"


def infer_op_metadata(
    func: Any,
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
    *,
    include_io_names: bool = False,
) -> dict[str, Any]:
    """
    Best-effort metadata for 5 sample heavy ops:
      - read_csv
      - merge/join
      - groupby agg/aggregate
      - concat
      - sort_values/sort_index

    Intentionally conservative and low-overhead.
    """
    meta: dict[str, Any] = {}

    callable_full = _callable_name(func)
    meta["callable"] = callable_full

    name = getattr(func, "__name__", "callable")
    self_obj = getattr(func, "__self__", None)

    if self_obj is not None:
        meta["operation"] = f"{type(self_obj).__name__}.{name}"
        lr, lc = _shape(self_obj)
        if lr is not None:
            meta["left_rows"] = lr
        if lc is not None:
            meta["left_cols"] = lc
    else:
        mod = getattr(func, "__module__", "") or ""
        meta["operation"] = callable_full if mod.startswith("pandas") else name

    # merge/join
    if name in ("merge", "join") and self_obj is not None:
        try:
            right = args[0] if args else (kwargs.get("right") or kwargs.get("other"))
            rr, rc = _shape(right)
            if rr is not None:
                meta["right_rows"] = rr
            if rc is not None:
                meta["right_cols"] = rc
            for k in (
                "how",
                "on",
                "left_on",
                "right_on",
                "left_index",
                "right_index",
                "sort",
                "validate",
            ):
                if k in kwargs:
                    meta[k] = _safe_value(kwargs.get(k))
        except Exception:
            pass

    # concat
    if name == "concat":
        try:
            objs = args[0] if args else kwargs.get("objs")
            n = _maybe_len(objs)
            if n is not None:
                meta["n_objs"] = n
            for k in (
                "axis",
                "join",
                "ignore_index",
                "keys",
                "sort",
                "verify_integrity",
                "copy",
            ):
                if k in kwargs:
                    meta[f"concat_{k}"] = _safe_value(kwargs.get(k))
        except Exception:
            pass

    # sort
    if name in ("sort_values", "sort_index"):
        try:
            for k in (
                "by",
                "axis",
                "ascending",
                "kind",
                "na_position",
                "ignore_index",
                "level",
                "sort_remaining",
            ):
                if k in kwargs:
                    meta[f"sort_{k}"] = _safe_value(kwargs.get(k))
        except Exception:
            pass

    # groupby agg
    if name in ("agg", "aggregate"):
        try:
            obj = getattr(self_obj, "obj", None)
            if obj is not None:
                r, c = _shape(obj)
                if r is not None:
                    meta.setdefault("left_rows", r)
                if c is not None:
                    meta.setdefault("left_cols", c)
            if args:
                meta["agg_spec"] = _safe_value(args[0])
        except Exception:
            pass

    # read_csv (optional basename only)
    if name == "read_csv":
        try:
            if include_io_names and args:
                meta["csv_file"] = _safe_basename(args[0])
            for k in (
                "nrows",
                "usecols",
                "dtype",
                "parse_dates",
                "engine",
                "memory_map",
                "low_memory",
            ):
                if k in kwargs:
                    meta[f"csv_{k}"] = _safe_value(kwargs.get(k))
        except Exception:
            pass

    return meta


__all__ = ["infer_op_metadata"]
