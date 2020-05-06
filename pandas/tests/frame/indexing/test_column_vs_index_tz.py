import pandas as pd
import pandas._testing as tm


class TestColumnvsIndexTZEquality:
    # https://github.com/pandas-dev/pandas/issues/29463
    def check_for_various_tz(self, tz):
        df = pd.DataFrame(
            {
                "val": range(10),
                "time": pd.date_range(start="2019-01-01", freq="1d", periods=10, tz=tz),
            }
        )
        df_query = df.query('"2019-01-03 00:00:00+00" < time')
        l1 = pd.DataFrame(list(df_query["time"]))

        # # This was earlier raising an exception.
        index_query = df.set_index("time").query('"2019-01-03 00:00:00+00" < time')
        l2 = pd.DataFrame(list(index_query.index))
        tm.assert_frame_equal(l1, l2)

    def test_check_column_vs_index_tz_query(self):
        tz_list = [
            "Africa/Abidjan",
            "Africa/Douala",
            "Africa/Mbabane",
            "America/Argentina/Catamarca",
            "America/Belize",
            "America/Curacao",
            "America/Guatemala",
            "America/Kentucky/Louisville",
            "America/Mexico_City",
            "America/Port-au-Prince",
            "America/Sitka",
            "Antarctica/Casey",
            "Asia/Ashkhabad",
            "Asia/Dubai",
            "Asia/Khandyga",
            "Asia/Qatar",
            "Asia/Tomsk",
            "Atlantic/Reykjavik",
            "Australia/Queensland",
            "Canada/Yukon",
            "Etc/GMT+7",
            "Etc/UCT",
            "Europe/Guernsey",
            "Europe/Paris",
            "Europe/Vienna",
            "Indian/Cocos",
            "NZ",
            "Pacific/Honolulu",
            "Pacific/Samoa",
            "US/Eastern",
        ]

        for tz in tz_list:
            self.check_for_various_tz(tz)
