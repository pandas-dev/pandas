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

- pivot_table

- DONE Set mixed-type values with .ix
- get_dtype_counts / dtypes
- save / load functions
- combine_first
- describe for Series
- DataFrame.to_string
- Index / MultiIndex names
- Unstack / stack by level name
- ignore_index in DataFrame.append
- Inner join on key
- Multi-key joining
- as_index=False in groupby
- is_monotonic
- isnull/notnull as instance methods
- name attribute on Series
- DataFrame.to_csv: different delimiters?
- groupby with level name
- MultiIndex
  - get_level_values

- Update to reflect Python 3 support in intro
- align functions
- df[col_list]
- Panel.rename_axis
- & and | for intersection / union
- IPython tab complete hook

Performance blog
----------------
- Series / Time series data alignment
- DataFrame alignment
- Groupby
- joining
- Take
