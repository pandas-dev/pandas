.. _contributing_environment:

{{ header }}

==================================
Creating a development environment
==================================

To test out code changes, you'll need to build pandas from source, which
requires a C/C++ compiler and Python environment. If you're making documentation
changes, you can skip to :ref:`contributing to the documentation <contributing_documentation>` but if you skip
creating the development environment you won't be able to build the documentation
locally before pushing your changes. It's recommended to also install the :ref:`pre-commit hooks <contributing.pre-commit>`.

.. toctree::
    :maxdepth: 2
    :hidden:


Step 1: install a C compiler
----------------------------

Installing a C compiler will depend on your operating system

Windows
~~~~~~~

You will need `Build Tools for Visual Studio 2026
<https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2026>`_.

.. note::
   If you encounter an error indicating ``cl.exe``,
   reopen the installer and also select the optional component
   **MSVC v142 - VS 2019 C++ x64/x86 build tools** in the right pane for installation.

Alternatively, you can install the necessary components on the commandline using
`vs_BuildTools.exe <https://learn.microsoft.com/en-us/visualstudio/install/use-command-line-parameters-to-install-visual-studio?source=recommendations&view=vs-2022>`_

Alternatively, you could use the `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_
and consult the ``Linux`` instructions below.

MacOS
~~~~~

To use the :ref:`conda <contributing.conda>`-based compilers, you will need to install the
Developer Tools using ``xcode-select --install``.

If you prefer to use a different compiler, general information can be found in the
`Python Developer's Guide <https://devguide.python.org/setup/#macos>`_.

Linux
~~~~~

For Linux-based :ref:`conda <contributing.conda>` installations, you won't have to install any
additional components outside of the conda environment. The instructions
below are only needed if your setup isn't based on conda environments.

Some Linux distributions will come with a pre-installed C compiler. To find out
which compilers (and versions) are installed on your system

.. code-block:: bash

    # for Debian/Ubuntu:
    dpkg --list | grep compiler
    # for Red Hat/RHEL/CentOS/Fedora:
    yum list installed | grep -i --color compiler

`GCC (GNU Compiler Collection) <https://gcc.gnu.org/>`_, is a widely used
compiler, which supports C and a number of other languages. If GCC is listed
as an installed compiler nothing more is required.

If no C compiler is installed, or you wish to upgrade, or you're using a different
Linux distribution, consult your favorite search engine for compiler installation/update
instructions.

Let us know if you have any difficulties by opening an issue or reaching out on our contributor
community :ref:`Slack <community.slack>`.

Step 2: create a virtual environment
------------------------------------

Before proceeding:

* Make sure that you have :ref:`cloned the repository <contributing.forking>`
* ``cd`` to the pandas source directory you just created with the clone command

.. _contributing.conda:

Option 1: using conda (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Install miniforge to get `conda <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_
* Create and activate the ``pandas-dev`` conda environment using the following commands:

.. code-block:: bash

   conda env create --file environment.yml
   conda activate pandas-dev

.. _contributing.pip:

Option 2: using pip
~~~~~~~~~~~~~~~~~~~

You'll need to have at least the :ref:`minimum Python version <install.version>` that pandas supports.

**Unix**/**macOS with virtualenv**

.. code-block:: bash

   # Create a virtual environment
   # Use an ENV_DIR of your choice. We'll use ~/virtualenvs/pandas-dev
   # Any parent directories should already exist
   python3 -m venv ~/virtualenvs/pandas-dev

   # Activate the virtualenv
   . ~/virtualenvs/pandas-dev/bin/activate

   # Install the build dependencies
   python -m pip install -r requirements-dev.txt

**Unix**/**macOS with pyenv**

Consult the `docs <https://github.com/pyenv/pyenv>`_ for setting up pyenv.

.. code-block:: bash

   # Create a virtual environment
   # Use an ENV_DIR of your choice. We'll use ~/Users/<yourname>/.pyenv/versions/pandas-dev
   pyenv virtualenv 3.11 pandas-dev

   # Activate the virtualenv
   pyenv activate pandas-dev

   # Now install the build dependencies in the cloned pandas repo
   python -m pip install -r requirements-dev.txt

**Windows**

Below is a brief overview on how to set-up a virtual environment with Powershell
under Windows. For details please refer to the
`official virtualenv user guide <https://virtualenv.pypa.io/en/latest/user_guide.html#activators>`__.

Use an ENV_DIR of your choice. We'll use ``~\\virtualenvs\\pandas-dev`` where
``~`` is the folder pointed to by either ``$env:USERPROFILE`` (Powershell) or
``%USERPROFILE%`` (cmd.exe) environment variable. Any parent directories
should already exist.

.. code-block:: powershell

   # Create a virtual environment
   python -m venv $env:USERPROFILE\virtualenvs\pandas-dev

   # Activate the virtualenv. Use activate.bat for cmd.exe
   ~\virtualenvs\pandas-dev\Scripts\Activate.ps1

   # Install the build dependencies
   python -m pip install -r requirements-dev.txt

Step 3: build and install pandas
--------------------------------

pandas uses the `Meson <https://mesonbuild.com/>`_ build backend via `PEP 517 <https://peps.python.org/pep-0517/>`_
to build the C extensions and install the library.

To compile and install pandas in editable mode, run:

.. code-block:: bash

   python -m pip install --verbose --editable . --no-build-isolation

Additional Meson options can be passed to the ``pip install`` command to modify the installation.
Helpful options include:

* ``-Ceditable-verbose=true``: Print verbose logs during rebuild, even during ``import``.
* ``-Cbuilddir="your builddir here"``: Specify a different build directory for the C extensions.
* ``-Csetup-args="-Ddebug=true"``: Compile the C extensions with debug symbols.

.. note::
   When pandas is installed in ``--editable`` mode, pandas will automatically rebuild the library upon ``import``,
   and build logs will show if ``-Ceditable-verbose=true`` is passed as well.

Now, pandas has been installed into your virtual environment, and the version number will
reflect that it's a development version with a reference to the latest Git hash from which pandas was built.

.. code-block:: ipython

   In [1]: import pandas

   # Your version will be structured similar, but not match, this example.
   In [2]: print(pandas.__version__)
   3.0.0.dev0+880.g2b9e661fbb.dirty

.. note::
   The version number is pulled from the latest repository tag, so ensure tags were fetched when
   :ref:`cloning the repository <contributing.forking>` from your fork.


Additionally, you can try :ref:`running the test suite <contributing.running_tests>`.
