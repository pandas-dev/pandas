from __future__ import annotations


def assert_eq(a, b, scheduler="sync"):
    if hasattr(a, "compute"):
        a = a.compute(scheduler=scheduler)
    if hasattr(b, "compute"):
        b = b.compute(scheduler=scheduler)

    assert a == b
