.. _debugging_c_extensions:

{{ header }}

**********************
Debugging C Extensions
**********************

Pandas uses select C extensions for high performance IO operations. In case you need to debug segfaults or general issues with those extensions, the following steps may be helpful. These steps are geared towards using lldb as a debugger, though the steps for gdb will be similar.

First, be sure to compile the extensions with the appropriate flags to generate debug symbols and remove optimizations. This can be achieved as follows:

.. code-block:: sh

   python setup.py build_ext --inplace -j4 --with-debugging-symbols

Next you can create a script that hits the extension module you are looking to debug and place it in the project root. Thereafter launch a Python process under lldb:

.. code-block:: sh

   lldb run python

If desired, set breakpoints at various file locations using the below syntax:

.. code-block:: sh

   breakpoint set --file pandas/_libs/src/ujson/python/objToJSON.c --line 1547

At this point you may get *WARNING:  Unable to resolve breakpoint to any actual locations.*. If you have not yet executed anything it is possible that this module has not been loaded into memory, which is why the location cannot be resolved. You can simply ignore for now as it will bind when we actually execute code.

Finally go ahead and execute your script:

.. code-block:: sh

   run <the_script>.py

Code execution will halt at the breakpoint defined or at the occurance of any segfault. LLDB's `GDB to LLDB command map <https://lldb.llvm.org/use/map.html>`_ provides a listing of debugger command that you can execute using either debugger.
