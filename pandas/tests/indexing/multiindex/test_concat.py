from pandas import (
    DataFrame,
    MultiIndex,
    concat,
)
import pandas._testing as tm

df1 = DataFrame({"col": ["a", "b", "c"]}, index=["1", "2", "2"])
df2 = concat([df1], keys=["X"])

iterables = [["X"], ["1", "2", "2"]]

result = df2.index
expected = MultiIndex.from_product(iterables)

tm.assert_index_equal(result, expected)

assert df2.index.has_duplicates == True
assert df2.index.is_unique == False
