Title: TST/ENH: Raise TypeError in Series.searchsorted for incomparable object-dtype values

Summary
-------
This small change makes Series.searchsorted raise a TypeError when the underlying
values are a numpy ndarray with dtype=object containing elements that are not
mutually comparable with the provided search value (for example, mixing int and
str). This aligns the behavior of `searchsorted` with `sort_values` and
reduces surprising cases where NumPy's `searchsorted` can return an index even
though comparisons between the types would fail.

Files changed
------------
- pandas/core/base.py
  - Add a lightweight runtime comparability check for object-dtype ndarrays in
    IndexOpsMixin.searchsorted. If a simple sample comparison between an array
    element and the search value raises TypeError, we propagate that TypeError.

- pandas/tests/series/methods/test_searchsorted.py
  - Add `test_searchsorted_incomparable_object_raises` which asserts that
    `Series([1, 2, "1"]).searchsorted("1")` raises TypeError.

Rationale
--------
Pandas delegates `searchsorted` to NumPy for ndarray-backed data. NumPy's
behavior on mixed-type object arrays can be surprising: it sometimes finds an
insertion index even when Python comparisons between element types would raise
TypeError (e.g. `1 < "1"`). Other pandas operations (like `sort_values`) raise
in that situation, so this change makes `searchsorted` consistent with the
rest of pandas.

Behavior and trade-offs
----------------------
- The comparability check is deliberately lightweight: it attempts a single
  comparison between the first non-NA array element and the sample search
  value. If that raises TypeError, we re-raise.
- This heuristic catches the common case (mixed ints/strings) without scanning
  the whole array (which would be expensive). It may not detect all
  pathological mixed-type arrays (for example, if the first element is
  comparable but later ones are not). If we want a stricter rule we can
  instead sample more elements or check types across the array, at some
  performance cost.

Testing
------
- New test added (see above). To run locally:

  # install in editable mode if importing from source
  python -m pip install -ve .

  # run the single test
  pytest -q pandas/tests/series/methods/test_searchsorted.py::test_searchsorted_incomparable_object_raises

Compatibility
------------
- Backwards compatible for numeric/datetime/etc. arrays: behavior unchanged.
- For object-dtype arrays with mixed types there is now a TypeError where
  previously NumPy might have silently returned an index. This is intentional
  to make behavior consistent with sorting.

Follow-ups
---------
- If desired, we can strengthen the comparability check (sample multiple
  elements or inspect the set of Python types) and add tests for those
  conditions.

PR checklist
-----------
- [ ] Add release note if desired (small change to searchsorted semantics)
- [ ] Add/adjust tests for stronger heuristics if implemented
