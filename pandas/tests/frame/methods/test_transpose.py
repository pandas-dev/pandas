import pandas as pd
import pandas._testing as tm


class TestTranspose:
    def test_transpose_tzaware_1col_single_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")

        df = pd.DataFrame(dti)
        assert (df.dtypes == dti.dtype).all()
        res = df.T
        assert (res.dtypes == dti.dtype).all()

    def test_transpose_tzaware_2col_single_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")

        df3 = pd.DataFrame({"A": dti, "B": dti})
        assert (df3.dtypes == dti.dtype).all()
        res3 = df3.T
        assert (res3.dtypes == dti.dtype).all()

    def test_transpose_tzaware_2col_mixed_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")
        dti2 = dti.tz_convert("US/Pacific")

        df4 = pd.DataFrame({"A": dti, "B": dti2})
        assert (df4.dtypes == [dti.dtype, dti2.dtype]).all()
        assert (df4.T.dtypes == object).all()
        tm.assert_frame_equal(df4.T.T, df4)

    def test_transpose_object_to_tzaware_mixed_tz(self):
        # GH#26825
        dti = pd.date_range("2016-04-05 04:30", periods=3, tz="UTC")
        dti2 = dti.tz_convert("US/Pacific")

        # mixed all-tzaware dtypes
        df2 = pd.DataFrame([dti, dti2])
        assert (df2.dtypes == object).all()
        res2 = df2.T
        assert (res2.dtypes == [dti.dtype, dti2.dtype]).all()
