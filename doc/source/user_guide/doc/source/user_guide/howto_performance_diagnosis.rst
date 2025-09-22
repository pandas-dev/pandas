.. _howto_performance_diagnosis:

Diagnosing performance bottlenecks
==================================

This page shows practical steps to *measure first*, then fix the most common
pandas bottlenecks (slow ``apply``, heavy memory usage, joins/groupbys, and
unintended dtype conversions). It complements :doc:`enhancingperf` with quick
checklists you can try in any notebook or script.

Quick checklist
---------------

- Inspect dtypes and memory: ``df.info(memory_usage="deep")`` and
  ``df.memory_usage(deep=True)``.
- Measure time with ``time.perf_counter`` or IPython's ``%timeit``.
- Prefer vectorized operations over ``DataFrame.apply`` / row loops.
- Ensure join/groupby keys have the same dtype; consider categoricals for
  low-cardinality keys.
- When arithmetic-heavy, consider :func:`pandas.eval` or moving work to
  specialized libraries (e.g., Numba) if appropriate.

Measure time
------------

.. code-block:: python

   import time
   start = time.perf_counter()

   # your operation here, e.g. df.groupby("key")["value"].mean()

   elapsed = time.perf_counter() - start
   print(f"{elapsed:.3f}s")  # wall-clock timing

In notebooks, ``%timeit`` provides robust micro-benchmarks (avoid inside docs
examples that execute on build).

Measure memory
--------------

.. code-block:: python

   df.info(memory_usage="deep")
   df.memory_usage(deep=True)

Look for large object-dtype columns; consider converting to ``string``,
nullable integer/float (``Int64`` / ``Float64``), or ``category``.

Vectorize instead of apply
--------------------------

.. code-block:: python

   # Slow: Python-level apply calls len() per row/cell
   s = df["name"]
   slow = s.apply(len)

   # Fast: vectorized string method or map over Python built-in
   fast = s.str.len()
   # or
   fast2 = s.map(len)

Joins and groupbys
------------------

- Align dtypes on keys before ``merge``: ``df1["key"] = df1["key"].astype("int64")``
  to match ``df2["key"]``.
- For low-cardinality keys, try categoricals:

.. code-block:: python

   df["key"] = df["key"].astype("category")
   out = df.groupby("key", observed=True)["value"].mean()

Arithmetic / expressions
------------------------

For column-wise arithmetic and boolean logic, ``pandas.eval`` can reduce
temporary objects and speed up some expressions:

.. code-block:: python

   df = df.eval("z = (x + y) * 2")

When to scale out
-----------------

If a single-machine DataFrame is too large or the workflow is inherently
parallel, consider external tools (e.g., Dask) or algorithmic changes. Keep
this page about *diagnosis*; see :doc:`enhancingperf` for advanced options.

See also
--------

- :doc:`enhancingperf`
- :doc:`categorical`
- :doc:`missing_data`
- :doc:`pyarrow`  # Arrow-backed dtypes and memory behavior

Notes for contributors
----------------------

Examples use ``.. code-block:: python`` to avoid executed doctests. Keep code
snippets small and runnable; prefer idiomatic pandas over micro-optimizations.
