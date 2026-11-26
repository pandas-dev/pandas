.. _glossary:

********
Glossary
********

.. glossary::
   :sorted:

   alignment
      Before an operation, pandas lines data up by index or column label, not by position. A label that appears on only one side produces a missing value.

   dtype
      The type of the values stored in a Series or DataFrame column, such as ``int64``, ``float64``, ``datetime64``, or ``category``. It determines how the data is stored and which operations work on it.

   file-like
      Any object with a ``read`` or ``write`` method, such as an open file. pandas I/O functions accept these in place of a file path.

   hashable
      An object whose hash value never changes during its lifetime, so it can be stored in a hash table. This is what lets it be a dictionary key or a set element. Strings, numbers, and tuples of hashable objects qualify, while lists and dicts do not.

   method chaining
      Calling a method directly on the result of the previous one, in succession. Each step usually returns a new object, so a whole pipeline can be written as one expression.

   MultiIndex
      An Index with more than one level, where each label is a tuple with one value per level. It allows rows to be selected by several keys.

   regular expression
      A pattern that describes a set of strings, accepted by string methods such as ``str.contains`` and ``str.extract``.

   reindex
      Building a new object with exactly the labels you ask for. Labels that already exist keep their values, new labels get missing values (or a fill value), and labels you leave out are dropped.

   UDF
   user-defined function
      A function you write yourself and hand to a pandas method such as ``apply``, ``map``, ``agg``, or ``transform`` to run on your data.

   vectorized
      Applied to a whole array at once, with the looping done in compiled code, not by a Python ``for`` loop over each element.