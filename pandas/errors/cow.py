from typing import Any

_chained_assignment_msg = (
    "A value is trying to be set on a copy of a DataFrame or Series "
    "through chained assignment.\n"
    "When using the Copy-on-Write mode, such chained assignment never works "
    "to update the original DataFrame or Series, because the intermediate "
    "object on which we are setting values always behaves as a copy.\n\n"
    "Try using '.loc[row_indexer, col_indexer] = value' instead, to perform "
    "the assignment in a single step.\n\n"
    "See the caveats in the documentation: "
    "https://pandas.pydata.org/pandas-docs/stable/user_guide/"
    "indexing.html#returning-a-view-versus-a-copy"
)


_chained_assignment_method_msg = (
    "A value is trying to be set on a copy of a DataFrame or Series "
    "through chained assignment using an inplace method.\n"
    "When using the Copy-on-Write mode, such inplace method never works "
    "to update the original DataFrame or Series, because the intermediate "
    "object on which we are setting values always behaves as a copy.\n\n"
    "For example, when doing 'df[col].method(value, inplace=True)', try "
    "using 'df.method({col: value}, inplace=True)' instead, to perform "
    "the operation inplace on the original object.\n\n"
)


_chained_assignment_warning_msg = (
    "ChainedAssignmentError: behaviour will change in pandas 3.0!\n"
    "You are setting values through chained assignment. Currently this works "
    "in certain cases, but when using Copy-on-Write (which will become the "
    "default behaviour in pandas 3.0) this will never work to update the "
    "original DataFrame or Series, because the intermediate object on which "
    "we are setting values will behave as a copy.\n"
    "A typical example is when you are setting values in a column of a "
    "DataFrame, like:\n\n"
    'df["col"][row_indexer] = value\n\n'
    'Use `df.loc[row_indexer, "col"] = values` instead, to perform the '
    "assignment in a single step and ensure this keeps updating the original `df`.\n\n"
    "See the caveats in the documentation: "
    "https://pandas.pydata.org/pandas-docs/stable/user_guide/"
    "indexing.html#returning-a-view-versus-a-copy\n"
)

_chained_assignment_warning_method_msg = (
    "A value is trying to be set on a copy of a DataFrame or Series "
    "through chained assignment using an inplace method.\n"
    "The behavior will change in pandas 3.0. This inplace method will "
    "never work because the intermediate object on which we are setting "
    "values always behaves as a copy.\n\n"
    "For example, when doing 'df[col].method(value, inplace=True)', try "
    "using 'df.method({col: value}, inplace=True)' or "
    "df[col] = df[col].method(value) instead, to perform "
    "the operation inplace on the original object.\n\n"
)


def _check_cacher(obj: Any) -> bool:
    # This is a mess, selection paths that return a view set the _cacher attribute
    # on the Series; most of them also set _item_cache which adds 1 to our relevant
    # reference count, but iloc does not, so we have to check if we are actually
    # in the item cache
    if hasattr(obj, "_cacher"):
        parent = obj._cacher[1]()
        # parent could be dead
        if parent is None:
            return False
        if hasattr(parent, "_item_cache"):
            if obj._cacher[0] in parent._item_cache:
                # Check if we are actually the item from item_cache, iloc creates a
                # new object
                return obj is parent._item_cache[obj._cacher[0]]
    return False
