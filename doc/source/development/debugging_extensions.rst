.. _debugging_c_extensions:

{{ header }}

======================
Debugging C extensions
======================

Pandas uses Cython and C/C++ `extension modules <https://docs.python.org/3/extending/extending.html>`_ to optimize performance. Unfortunately, the standard Python debugger does not allow you to step into these extensions. Cython extensions can be debugged with the `Cython debugger <https://docs.cython.org/en/latest/src/userguide/debugging.html>`_ and C/C++ extensions can be debugged using the tools shipped with your platform's compiler.

For Python developers with limited or no C/C++ experience this can seem a daunting task. Core developer Will Ayd has written a 3 part blog series to help guide you from the standard Python debugger into these other tools:

  1. `Fundamental Python Debugging Part 1 - Python <https://willayd.com/fundamental-python-debugging-part-1-python.html>`_
  2. `Fundamental Python Debugging Part 2 - Python Extensions <https://willayd.com/fundamental-python-debugging-part-2-python-extensions.html>`_
  3. `Fundamental Python Debugging Part 3 - Cython Extensions <https://willayd.com/fundamental-python-debugging-part-3-cython-extensions.html>`_
