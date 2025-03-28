import pytest

from xarray.util.deprecation_helpers import _deprecate_positional_args


def test_deprecate_positional_args_warns_for_function():
    @_deprecate_positional_args("v0.1")
    def f1(a, b, *, c="c", d="d"):
        return a, b, c, d

    result = f1(1, 2)
    assert result == (1, 2, "c", "d")

    result = f1(1, 2, c=3, d=4)
    assert result == (1, 2, 3, 4)

    with pytest.warns(FutureWarning, match=r".*v0.1"):
        result = f1(1, 2, 3)  # type: ignore[misc]
    assert result == (1, 2, 3, "d")

    with pytest.warns(FutureWarning, match=r"Passing 'c' as positional"):
        result = f1(1, 2, 3)  # type: ignore[misc]
    assert result == (1, 2, 3, "d")

    with pytest.warns(FutureWarning, match=r"Passing 'c, d' as positional"):
        result = f1(1, 2, 3, 4)  # type: ignore[misc]
    assert result == (1, 2, 3, 4)

    @_deprecate_positional_args("v0.1")
    def f2(a="a", *, b="b", c="c", d="d"):
        return a, b, c, d

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = f2(1, 2)  # type: ignore[misc]
    assert result == (1, 2, "c", "d")

    @_deprecate_positional_args("v0.1")
    def f3(a, *, b="b", **kwargs):
        return a, b, kwargs

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = f3(1, 2)  # type: ignore[misc]
    assert result == (1, 2, {})

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = f3(1, 2, f="f")  # type: ignore[misc]
    assert result == (1, 2, {"f": "f"})

    @_deprecate_positional_args("v0.1")
    def f4(a, /, *, b="b", **kwargs):
        return a, b, kwargs

    result = f4(1)
    assert result == (1, "b", {})

    result = f4(1, b=2, f="f")
    assert result == (1, 2, {"f": "f"})

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = f4(1, 2, f="f")  # type: ignore[misc]
    assert result == (1, 2, {"f": "f"})

    with pytest.raises(TypeError, match=r"Keyword-only param without default"):

        @_deprecate_positional_args("v0.1")
        def f5(a, *, b, c=3, **kwargs):
            pass


def test_deprecate_positional_args_warns_for_class():
    class A1:
        @_deprecate_positional_args("v0.1")
        def method(self, a, b, *, c="c", d="d"):
            return a, b, c, d

    result = A1().method(1, 2)
    assert result == (1, 2, "c", "d")

    result = A1().method(1, 2, c=3, d=4)
    assert result == (1, 2, 3, 4)

    with pytest.warns(FutureWarning, match=r".*v0.1"):
        result = A1().method(1, 2, 3)  # type: ignore[misc]
    assert result == (1, 2, 3, "d")

    with pytest.warns(FutureWarning, match=r"Passing 'c' as positional"):
        result = A1().method(1, 2, 3)  # type: ignore[misc]
    assert result == (1, 2, 3, "d")

    with pytest.warns(FutureWarning, match=r"Passing 'c, d' as positional"):
        result = A1().method(1, 2, 3, 4)  # type: ignore[misc]
    assert result == (1, 2, 3, 4)

    class A2:
        @_deprecate_positional_args("v0.1")
        def method(self, a=1, b=1, *, c="c", d="d"):
            return a, b, c, d

    with pytest.warns(FutureWarning, match=r"Passing 'c' as positional"):
        result = A2().method(1, 2, 3)  # type: ignore[misc]
    assert result == (1, 2, 3, "d")

    with pytest.warns(FutureWarning, match=r"Passing 'c, d' as positional"):
        result = A2().method(1, 2, 3, 4)  # type: ignore[misc]
    assert result == (1, 2, 3, 4)

    class A3:
        @_deprecate_positional_args("v0.1")
        def method(self, a, *, b="b", **kwargs):
            return a, b, kwargs

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = A3().method(1, 2)  # type: ignore[misc]
    assert result == (1, 2, {})

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = A3().method(1, 2, f="f")  # type: ignore[misc]
    assert result == (1, 2, {"f": "f"})

    class A4:
        @_deprecate_positional_args("v0.1")
        def method(self, a, /, *, b="b", **kwargs):
            return a, b, kwargs

    result = A4().method(1)
    assert result == (1, "b", {})

    result = A4().method(1, b=2, f="f")
    assert result == (1, 2, {"f": "f"})

    with pytest.warns(FutureWarning, match=r"Passing 'b' as positional"):
        result = A4().method(1, 2, f="f")  # type: ignore[misc]
    assert result == (1, 2, {"f": "f"})

    with pytest.raises(TypeError, match=r"Keyword-only param without default"):

        class A5:
            @_deprecate_positional_args("v0.1")
            def __init__(self, a, *, b, c=3, **kwargs):
                pass
