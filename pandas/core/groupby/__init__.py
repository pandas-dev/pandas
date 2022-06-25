import warnings

from pandas._core import groupby as _groupby

for name in dir(_groupby):
    assert name != "_groupby"
    globals()[name] = getattr(_groupby, name)

del _groupby

warnings.warn(
    "pandas.core.groupby is deprecated and will be moved to pandas._core.groupby "
    "in a future version. Import from the public API instead.",
    category=FutureWarning,
    stacklevel=2,
)
