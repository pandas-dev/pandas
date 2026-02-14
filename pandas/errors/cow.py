_chained_assignment_msg = (
    "A value is being set on a copy of a DataFrame or Series "
    "through chained assignment.\n"
    "Such chained assignment never works to update the original DataFrame or "
    "Series, because the intermediate object on which we are setting values "
    "always behaves as a copy (due to Copy-on-Write).\n\n"
    "Try using '.loc[row_indexer, col_indexer] = value' instead, to perform "
    "the assignment in a single step.\n\n"
    "See the documentation for a more detailed explanation: "
    "https://pandas.pydata.org/pandas-docs/stable/user_guide/"
    "copy_on_write.html#chained-assignment"
)


_chained_assignment_method_msg = (
    "A value is being set on a copy of a DataFrame or Series "
    "through chained assignment using an inplace method.\n"
    "Such inplace method never works to update the original DataFrame or Series, "
    "because the intermediate object on which we are setting values always "
    "behaves as a copy (due to Copy-on-Write).\n\n"
    "For example, when doing 'df[col].method(value, inplace=True)', try "
    "using 'df.method({col: value}, inplace=True)' instead, to perform "
    "the operation inplace on the original object, or try to avoid an inplace "
    "operation using 'df[col] = df[col].method(value)'.\n\n"
    "See the documentation for a more detailed explanation: "
    "https://pandas.pydata.org/pandas-docs/stable/user_guide/"
    "copy_on_write.html"
)


_chained_assignment_method_update_msg = (
    "A value is being set on a copy of a DataFrame or Series "
    "through chained assignment using an inplace method.\n"
    "Such inplace method never works to update the original DataFrame or Series, "
    "because the intermediate object on which we are setting values always "
    "behaves as a copy (due to Copy-on-Write).\n\n"
    "For example, when doing 'df[col].update(other)', try "
    "using 'df.update({col: other})' instead, to perform "
    "the operation inplace on the original object.\n\n"
    "See the documentation for a more detailed explanation: "
    "https://pandas.pydata.org/pandas-docs/stable/user_guide/"
    "copy_on_write.html"
)
