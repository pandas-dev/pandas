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

Debugging locally
-----------------

By default building pandas from source will generate a release build. To generate a development build you can type::

    pip install -ve . --no-build-isolation --config-settings=builddir="debug" --config-settings=setup-args="-Dbuildtype=debug"

.. note::

   conda environements update CFLAGS/CPPFLAGS with flags that are geared towards generating releases. If using conda, you may need to set ``CFLAGS="$CFLAGS -O0"`` and ``CPPFLAGS="$CPPFLAGS -O0"`` to ensure optimizations are turned off for debugging

By specifying ``builddir="debug"`` all of the targets will be built and placed in the debug directory relative to the project root. This helps to keep your debug and release artifacts separate; you are of course able to choose a different directory name or omit altogether if you do not care to separate build types.

Using Docker
------------

To simplify the debugging process, pandas has created a Docker image with a debug build of Python and the gdb/Cython debuggers pre-installed. You may either ``docker pull pandas/pandas-debug`` to get access to this image or build it from the ``tooling/debug`` folder locallly.

You can then mount your pandas repository into this image via:

.. code-block:: sh

   docker run --rm -it -w /data -v ${PWD}:/data pandas/pandas-debug

Inside the image, you can use meson to build/install pandas and place the build artifacts into a ``debug`` folder using a command as follows:

.. code-block:: sh

    python -m pip install -ve . --no-build-isolation --config-settings=builddir="debug" --config-settings=setup-args="-Dbuildtype=debug"

If planning to use cygdb, the files required by that application are placed within the build folder. So you have to first ``cd`` to the build folder, then start that application.

.. code-block:: sh

   cd debug
   cygdb

Within the debugger you can use `cygdb commands <https://docs.cython.org/en/latest/src/userguide/debugging.html#using-the-debugger>`_ to navigate cython extensions.

Editor support
--------------

The meson build system generates a `compilation database <https://clang.llvm.org/docs/JSONCompilationDatabase.html>`_ automatically and places it in the build directory. Many language servers and IDEs can use this information to provide code-completion, go-to-definition and error checking support as you type.

How each language server / IDE chooses to look for the compilation database may vary. When in doubt you may want to create a symlink at the root of the project that points to the compilation database in your build directory. Assuming you used *debug* as your directory name, you can run::

    ln -s debug/compile_commands.json .
