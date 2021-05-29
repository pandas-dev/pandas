import pandas as pd
import pandas._testing as tm


class TestMap:
    def test_map(self):
        # map test on StringDType, GH#40823
        df1 = pd.DataFrame(
            {"col1": [pd.NA, "foo", "bar"]},
            index=["id1", "id2", "id3"],
            dtype=pd.StringDtype(),
        )

        df2 = pd.DataFrame({"id": ["id4", "id2", "id1"]}, dtype=pd.StringDtype())

        df2["col1"] = df2["id"].map(df1["col1"])

        result = df2
        expected = pd.DataFrame(
            {"id": ["id4", "id2", "id1"], "col1": [pd.NA, "foo", pd.NA]},
            dtype=pd.StringDtype(),
        )

        tm.assert_frame_equal(result, expected)
