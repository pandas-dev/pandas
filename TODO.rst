DONE
----
- SparseSeries name integration + tests
- Refactor Series.repr

TODO
----
- _consolidate, does it always copy?
- Series.align with fill method. Will have to generate more Cython code

TODO docs
---------

- DONE read_csv / read_table
  - auto-sniff delimiter
  - MultiIndex
  - generally more documentation
- DONE pivot_table
- DONE Set mixed-type values with .ix
- DONE get_dtype_counts / dtypes
- DONE save / load functions
- DONE isnull/notnull as instance methods
- DONE DataFrame.to_string
- DONE IPython tab complete hook
- DONE ignore_index in DataFrame.append
- DONE describe for Series with dtype=object
- DONE as_index=False in groupby
- DONOTWANT is_monotonic
- DONE DataFrame.to_csv: different delimiters
- DONE combine_first
- DONE groupby with level name
- DONE MultiIndex get_level_values
- DONE & and | for intersection / union
- DONE Update to reflect Python 3 support in intro
- DONE Index / MultiIndex names
- DONE Unstack / stack by level name
- DONE name attribute on Series
- DONE Multi-key joining
- DONE Inner join on key
- DONE align functions
- DONE df[col_list]
- DONE Panel.rename_axis

Performance blog
----------------
- Series / Time series data alignment
- DataFrame alignment
- Groupby
- joining
- Take
