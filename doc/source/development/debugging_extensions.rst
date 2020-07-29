.. _debugging_c_extensions:

{{ header }}

======================
Debugging C extensions
======================

Pandas uses select C extensions for high performance IO operations. In case you need to debug segfaults or general issues with those extensions, the following steps may be helpful. These steps are geared towards using lldb as a debugger, though the steps for gdb will be similar.

First, be sure to compile the extensions with the appropriate flags to generate debug symbols and remove optimizations. This can be achieved as follows:

.. code-block:: sh

   python setup.py build_ext --inplace -j4 --with-debugging-symbols

Using a debugger
================

You can create a script that hits the extension module you are looking to debug and place it in the project root. Thereafter launch a Python process under lldb:

.. code-block:: sh

   lldb python

If desired, set breakpoints at various file locations using the below syntax:

.. code-block:: sh

   breakpoint set --file pandas/_libs/src/ujson/python/objToJSON.c --line 1547

At this point you may get *WARNING:  Unable to resolve breakpoint to any actual locations.*. If you have not yet executed anything it is possible that this module has not been loaded into memory, which is why the location cannot be resolved. You can simply ignore for now as it will bind when we actually execute code.

Finally go ahead and execute your script:

.. code-block:: sh

   run <the_script>.py

Code execution will halt at the breakpoint defined or at the occurance of any segfault. LLDB's `GDB to LLDB command map <https://lldb.llvm.org/use/map.html>`_ provides a listing of debugger command that you can execute using either debugger.

Another option to execute the entire test suite under the debugger would be to run the following:

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
