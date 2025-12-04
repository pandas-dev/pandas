pandas.API.typing.expression
============================

.. currentmodule:: pandas.api.typing

.. autoclass:: Expression
   :members:
   :undoc-members:
   :show-inheritance:

**Description**

An ``Expression`` represents a symbolic reference to a DataFrame column.
It can be used inside functions like ``DataFrame.assign`` or ``DataFrame.loc`` to
refer to columns in a declarative way.

Example
-------

.. code-block:: python

   import pandas as pd

   df = pd.DataFrame({"a": [1, 2, 3]})
   df.assign(b=pd.col("a") + 5)

**Supported Operations**

- Arithmetic: ``+``, ``-``, ``*``, ``/``
- Comparison: ``<``, ``>``, ``==``, ``!=``
- Universal functions (ufuncs): can be applied directly
- Series accessors: ``.dt``, ``.str``, ``.cat``
- Series methods: ``.sum()``, ``.mean()``, ``.astype()``, etc.
