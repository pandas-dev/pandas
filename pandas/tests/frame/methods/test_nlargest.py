"""
Note: for naming purposes, most tests are title with as e.g. "test_nlargest_foo"
but are implicitly also testing nsmallest_foo.
"""

class TestNLargestNSmallest:

    # ----------------------------------------------------------------------
    # Top / bottom
    @pytest.mark.parametrize(
        "order",
        [
            ["a"],
            ["c"],
            ["a", "b"],
            ["a", "c"],
            ["b", "a"],
            ["b", "c"],
            ["a", "b", "c"],
            ["c", "a", "b"],
            ["c", "b", "a"],
            ["b", "c", "a"],
            ["b", "a", "c"],
            # dups!
            ["b", "c", "c"],
        ],
    )
    @pytest.mark.parametrize("n", range(1, 11))
    def test_n(self, df_strings, nselect_method, n, order):
        # GH 10393
        df = df_strings
        if "b" in order:

            error_msg = (
                f"Column 'b' has dtype object, "
                f"cannot use method '{nselect_method}' with this dtype"
            )
            with pytest.raises(TypeError, match=error_msg):
                getattr(df, nselect_method)(n, order)
        else:
            ascending = nselect_method == "nsmallest"
            result = getattr(df, nselect_method)(n, order)
            expected = df.sort_values(order, ascending=ascending).head(n)
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "columns", [["group", "category_string"], ["group", "string"]]
    )
    def test_n_error(self, df_main_dtypes, nselect_method, columns):
        df = df_main_dtypes
        col = columns[1]
        error_msg = (
            f"Column '{col}' has dtype {df[col].dtype}, "
            f"cannot use method '{nselect_method}' with this dtype"
        )
        # escape some characters that may be in the repr
        error_msg = (
            error_msg.replace("(", "\\(")
            .replace(")", "\\)")
            .replace("[", "\\[")
            .replace("]", "\\]")
        )
        with pytest.raises(TypeError, match=error_msg):
            getattr(df, nselect_method)(2, columns)

    def test_n_all_dtypes(self, df_main_dtypes):
        df = df_main_dtypes
        df.nsmallest(2, list(set(df) - {"category_string", "string"}))
        df.nlargest(2, list(set(df) - {"category_string", "string"}))

    @pytest.mark.parametrize(
        "method,expected",
        [
            (
                "nlargest",
                pd.DataFrame(
                    {"a": [2, 2, 2, 1], "b": [3, 2, 1, 3]}, index=[2, 1, 0, 3]
                ),
            ),
            (
                "nsmallest",
                pd.DataFrame(
                    {"a": [1, 1, 1, 2], "b": [1, 2, 3, 1]}, index=[5, 4, 3, 0]
                ),
            ),
        ],
    )
    def test_duplicates_on_starter_columns(self, method, expected):
        # regression test for #22752

        df = pd.DataFrame({"a": [2, 2, 2, 1, 1, 1], "b": [1, 2, 3, 3, 2, 1]})

        result = getattr(df, method)(4, columns=["a", "b"])
        tm.assert_frame_equal(result, expected)

    def test_n_identical_values(self):
        # GH 15297
        df = pd.DataFrame({"a": [1] * 5, "b": [1, 2, 3, 4, 5]})

        result = df.nlargest(3, "a")
        expected = pd.DataFrame({"a": [1] * 3, "b": [1, 2, 3]}, index=[0, 1, 2])
        tm.assert_frame_equal(result, expected)

        result = df.nsmallest(3, "a")
        expected = pd.DataFrame({"a": [1] * 3, "b": [1, 2, 3]})
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "order",
        [["a", "b", "c"], ["c", "b", "a"], ["a"], ["b"], ["a", "b"], ["c", "b"]],
    )
    @pytest.mark.parametrize("n", range(1, 6))
    def test_n_duplicate_index(self, df_duplicates, n, order):
        # GH 13412

        df = df_duplicates
        result = df.nsmallest(n, order)
        expected = df.sort_values(order).head(n)
        tm.assert_frame_equal(result, expected)

        result = df.nlargest(n, order)
        expected = df.sort_values(order, ascending=False).head(n)
        tm.assert_frame_equal(result, expected)

    def test_duplicate_keep_all_ties(self):
        # GH 16818
        df = pd.DataFrame(
            {"a": [5, 4, 4, 2, 3, 3, 3, 3], "b": [10, 9, 8, 7, 5, 50, 10, 20]}
        )
        result = df.nlargest(4, "a", keep="all")
        expected = pd.DataFrame(
            {
                "a": {0: 5, 1: 4, 2: 4, 4: 3, 5: 3, 6: 3, 7: 3},
                "b": {0: 10, 1: 9, 2: 8, 4: 5, 5: 50, 6: 10, 7: 20},
            }
        )
        tm.assert_frame_equal(result, expected)

        result = df.nsmallest(2, "a", keep="all")
        expected = pd.DataFrame(
            {
                "a": {3: 2, 4: 3, 5: 3, 6: 3, 7: 3},
                "b": {3: 7, 4: 5, 5: 50, 6: 10, 7: 20},
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_series_broadcasting(self):
        # smoke test for numpy warnings
        # GH 16378, GH 16306
        df = DataFrame([1.0, 1.0, 1.0])
        df_nan = DataFrame({"A": [np.nan, 2.0, np.nan]})
        s = Series([1, 1, 1])
        s_nan = Series([np.nan, np.nan, 1])

        with tm.assert_produces_warning(None):
            df_nan.clip(lower=s, axis=0)
            for op in ["lt", "le", "gt", "ge", "eq", "ne"]:
                getattr(df, op)(s_nan, axis=0)

    def test_series_nat_conversion(self):
        # GH 18521
        # Check rank does not mutate DataFrame
        df = DataFrame(np.random.randn(10, 3), dtype="float64")
        expected = df.copy()
        df.rank()
        result = df
        tm.assert_frame_equal(result, expected)

    def test_multiindex_column_lookup(self):
        # Check whether tuples are correctly treated as multi-level lookups.
        # GH 23033
        df = pd.DataFrame(
            columns=pd.MultiIndex.from_product([["x"], ["a", "b"]]),
            data=[[0.33, 0.13], [0.86, 0.25], [0.25, 0.70], [0.85, 0.91]],
        )

        # nsmallest
        result = df.nsmallest(3, ("x", "a"))
        expected = df.iloc[[2, 0, 3]]
        tm.assert_frame_equal(result, expected)

        # nlargest
        result = df.nlargest(3, ("x", "b"))
        expected = df.iloc[[3, 2, 1]]
        tm.assert_frame_equal(result, expected)
