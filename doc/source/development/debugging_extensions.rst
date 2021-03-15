.. _debugging_c_extensions:

{{ header }}

======================
Debugging C extensions
======================

Pandas uses select C extensions for high performance IO operations. In case you need to debug segfaults or general issues with those extensions, the following steps may be helpful.

First, be sure to compile the extensions with the appropriate flags to generate debug symbols and remove optimizations. This can be achieved as follows:

.. code-block:: sh

   python setup.py build_ext --inplace -j4 --with-debugging-symbols

Using a debugger
================

Assuming you are on a Unix-like operating system, you can use either lldb or gdb to debug. The choice between either is largely dependent on your compilation toolchain - typically you would use lldb if using clang and gdb if using gcc. For macOS users, please note that ``gcc`` is on modern systems an alias for ``clang``, so if using Xcode you usually opt for lldb. Regardless of which debugger you choose, please refer to your operating systems instructions on how to install.

After installing a debugger you can create a script that hits the extension module you are looking to debug. For demonstration purposes, let's assume you have a script called ``debug_testing.py`` with the following contents:

.. code-block:: python

   import pandas as pd

   pd.DataFrame([[1, 2]]).to_json()

Place the ``debug_testing.py`` script in the project root and launch a Python process under your debugger. If using lldb:

.. code-block:: sh

   lldb python

If using gdb:

.. code-block:: sh

   gdb python

Before executing our script, let's set a breakpoint in our JSON serializer in its entry function called ``objToJSON``. The lldb syntax would look as follows:

.. code-block:: sh

   breakpoint set --name objToJSON

Similarly for gdb:

.. code-block:: sh

   break objToJSON

.. note::

   You may get a warning that this breakpoint cannot be resolved in lldb. gdb may give a similar warning and prompt you to make the breakpoint on a future library load, which you should say yes to. This should only happen on the very first invocation as the module you wish to debug has not yet been loaded into memory.

Now go ahead and execute your script:

.. code-block:: sh

   run <the_script>.py

Code execution will halt at the breakpoint defined or at the occurrence of any segfault. LLDB's `GDB to LLDB command map <https://lldb.llvm.org/use/map.html>`_ provides a listing of debugger command that you can execute using either debugger.

Another option to execute the entire test suite under lldb would be to run the following:

.. code-block:: sh

   lldb -- python -m pytest

Or for gdb

.. code-block:: sh

   gdb --args python -m pytest

Once the process launches, simply type ``run`` and the test suite will begin, stopping at any segmentation fault that may occur.

Checking memory leaks with valgrind
===================================

You can use `Valgrind <https://www.valgrind.org>`_ to check for and log memory leaks in extensions. For instance, to check for a memory leak in a test from the suite you can run:

.. code-block:: sh

   PYTHONMALLOC=malloc valgrind --leak-check=yes --track-origins=yes --log-file=valgrind-log.txt python -m pytest <path_to_a_test>

Note that code execution under valgrind will take much longer than usual. While you can run valgrind against extensions compiled with any optimization level, it is suggested to have optimizations turned off from compiled extensions to reduce the amount of false positives. The ``--with-debugging-symbols`` flag passed during package setup will do this for you automatically.

.. note::

   For best results, you should run use a Python installation configured with Valgrind support (--with-valgrind)
