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

Generating debug builds
-----------------------

By default building pandas from source will generate a release build. To generate a development build you can type::

    pip install -ve . --no-build-isolation --config-settings=builddir="debug" --config-settings=setup-args="-Dbuildtype=debug"

By specifying ``builddir="debug"`` all of the targets will be built and placed in the debug directory relative to the project root. This helps to keep your debug and release artifacts separate; you are of course able to choose a different directory name or omit altogether if you do not care to separate build types.

Editor support
--------------

The meson build system generates a `compilation database <https://clang.llvm.org/docs/JSONCompilationDatabase.html>`_ automatically and places it in the build directory. Many language servers and IDEs can use this information to provide code-completion, go-to-defintion and error checking support as you type.

How each language server / IDE chooses to look for the compilation database may vary. When in doubt you may want to create a symlink at the root of the project that points to the compilation database in your build directory. Assuming you used *debug* as your directory name, you can run::

    ln -s debug/compile_commands.json .
