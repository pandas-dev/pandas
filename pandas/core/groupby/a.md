Purely integer-location based indexing for selection by position per group.

``.iloc[]`` is primarily integer position based (from ``0`` to
``length-1`` of the axis),

Allowed inputs for the first index are:

- An integer, e.g. ``5``.
- A slice object with ints and positive step, e.g. ``1:``, ``4:-3:2``.

Allowed inputs for the second index are as for DataFrame.iloc, namely:

- An integer, e.g. ``5``.
- A list or array of integers, e.g. ``[4, 3, 0]``.
- A slice object with ints, e.g. ``1:7``.
- A boolean array.
- A ``callable`` function with one argument (the calling Series or
  DataFrame) and that returns valid output for indexing (one of the above).

The output format is the same as GroupBy.head and GroupBy.tail, namely a subset
of the original DataFrame or Series with the index and order preserved.

The effect of ``grouped.iloc[i:j, k:l]`` is similar to

    grouped.apply(lambda x: x.iloc[i:j, k:l])

but very much faster and preserving the original index and order.

The behaviour is different from GroupBy.take:
- Input to iloc is a slice of indexes rather than a list of indexes.
- Output from iloc is:
    - In the same order as the original grouped DataFrame or Series.
    - Has the same index columns as the original grouped DataFrame or Series.
      (GroupBy.take introduces an additional index)
- GroupBy.take is extremely slow when there is a high group count.

The behaviour is different from GroupBy.nth:
- Input to iloc is a slice of indexes rather than a list of indexes.
- Output from iloc is:
    - In the same order as the original grouped DataFrame or Series.
    - Has the same index columns as the original grouped DataFrame or Series.
      (nth removes the grouped index)
- GroupBy.nth is quite fast for a high group count but the processing time
  grows with the length of the list of indexes.

Since GroupBy.take and GroupBy.nth only accept a list of individual indexes
it is not possible to define a slice that ends relative to the last row of
each group.

An important use case for GroupBy.iloc is a multi-indexed DataFrame with a
large primary index (Date, say) and a secondary index sorted to a different
order for each Date.
To reduce the DataFrame to a middle slice of each Date:

    df.groupby("Date").iloc[5:-5]

This returns a subset of df containing just the middle rows for each Date
and with its original order and indexing preserved.
(See test_multiindex() in tests/groupby/test_groupby_iloc.py)

Returns
-------
Series or DataFrame

See Also
--------
DataFrame.iloc : Purely integer-location based indexing for selection by position.
GroupBy.head : Return first n rows of each group.
GroupBy.tail : Return last n rows of each group.
GroupBy.nth : Take the nth row from each group if n is an int, or a subset of rows
if n is a list of ints.
GroupBy.take : Return the elements in the given positional indices along an axis.

Examples
--------
    >>> df = pd.DataFrame([["a", 1], ["a", 2], ["a", 3], ["b", 4], ["b", 5]],
    ...                   columns=["A", "B"])
    >>> df.groupby("A").iloc[1:2]
       A  B
    1  a  2
    4  b  5
    >>> df.groupby("A").iloc[:-1, -1:]
       B
    0  1
    1  2
    3  4