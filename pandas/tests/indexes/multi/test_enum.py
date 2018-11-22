import pandas as pd

import pandas.util.testing as tm
import pandas.util._test_decorators as td


@td.skip_if_no('enum')
def test_enum_in_multiindex():
    # GH 21298
    # Allow use of Enums as one of the factors in a MultiIndex.
    from enum import Enum
    MyEnum = Enum("MyEnum", "A B")
    df = pd.DataFrame(columns=pd.MultiIndex.from_product(iterables=[
        MyEnum,
        [1, 2]
    ]))

    exp_index_0 = pd.Index([MyEnum.A, MyEnum.B], dtype='object')
    tm.assert_index_equal(df.columns.levels[0], exp_index_0)

    expected = df.copy()
    df = df.append({(MyEnum.A, 1): "abc", (MyEnum.B, 2): "xyz"},
                   ignore_index=True)
    expected.loc[0, [(MyEnum.A, 1), (MyEnum.B, 2)]] = 'abc', 'xyz'
    tm.assert_frame_equal(df, expected)
