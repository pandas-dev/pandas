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

    contributing_gitpod.rst

Step 1: install a C compiler
----------------------------

How to do this will depend on your platform. If you choose to use ``Docker`` or ``GitPod``
in the next step, then you can skip this step.

**Windows**

You will need `Build Tools for Visual Studio 2022
<https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022>`_.

.. note::
        You DO NOT need to install Visual Studio 2022.
        You only need "Build Tools for Visual Studio 2022" found by
        scrolling down to "All downloads" -> "Tools for Visual Studio".
        In the installer, select the "Desktop development with C++" Workloads.

Alternatively, you can install the necessary components on the commandline using
`vs_BuildTools.exe <https://learn.microsoft.com/en-us/visualstudio/install/use-command-line-parameters-to-install-visual-studio?source=recommendations&view=vs-2022>`_

Alternatively, you could use the `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_
and consult the ``Linux`` instructions below.

**macOS**

To use the :ref:`mamba <contributing.mamba>`-based compilers, you will need to install the
Developer Tools using ``xcode-select --install``. Otherwise
information about compiler installation can be found here:
https://devguide.python.org/setup/#macos

**Linux**

For Linux-based :ref:`mamba <contributing.mamba>` installations, you won't have to install any
additional components outside of the mamba environment. The instructions
below are only needed if your setup isn't based on mamba environments.

Some Linux distributions will come with a pre-installed C compiler. To find out
which compilers (and versions) are installed on your system::

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

Step 2: create an isolated environment
----------------------------------------

Before we begin, please:

* Make sure that you have :any:`cloned the repository <contributing.forking>`
* ``cd`` to the pandas source directory you just created with the clone command

.. _contributing.mamba:

Option 1: using mamba (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Install `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_
* Make sure your mamba is up to date (``mamba update mamba``)

.. code-block:: none

   # Create and activate the build environment
   mamba env create --file environment.yml
   mamba activate pandas-dev

.. _contributing.pip:

Option 2: using pip
~~~~~~~~~~~~~~~~~~~

You'll need to have at least the :ref:`minimum Python version <install.version>` that pandas supports.
You also need to have ``setuptools`` 51.0.0 or later to build pandas.

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

Consult the docs for setting up pyenv `here <https://github.com/pyenv/pyenv>`__.

.. code-block:: bash

   # Create a virtual environment
   # Use an ENV_DIR of your choice. We'll use ~/Users/<yourname>/.pyenv/versions/pandas-dev
   pyenv virtualenv <version> <name-to-give-it>

   # For instance:
   pyenv virtualenv 3.9.10 pandas-dev

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

Option 3: using Docker
~~~~~~~~~~~~~~~~~~~~~~

pandas provides a ``DockerFile`` in the root directory to build a Docker image
with a full pandas development environment.

**Docker Commands**

Build the Docker image::

    # Build the image
    docker build -t pandas-dev .

Run Container::

    # Run a container and bind your local repo to the container
    # This command assumes you are running from your local repo
    # but if not alter ${PWD} to match your local repo path
    docker run -it --rm -v ${PWD}:/home/pandas pandas-dev

*Even easier, you can integrate Docker with the following IDEs:*

**Visual Studio Code**

You can use the DockerFile to launch a remote session with Visual Studio Code,
a popular free IDE, using the ``.devcontainer.json`` file.
See https://code.visualstudio.com/docs/remote/containers for details.

**PyCharm (Professional)**

Enable Docker support and use the Services tool window to build and manage images as well as
run and interact with containers.
See https://www.jetbrains.com/help/pycharm/docker.html for details.

Option 4: using Gitpod
~~~~~~~~~~~~~~~~~~~~~~

Gitpod is an open-source platform that automatically creates the correct development
environment right in your browser, reducing the need to install local development
environments and deal with incompatible dependencies.

If you are a Windows user, unfamiliar with using the command line or building pandas
for the first time, it is often faster to build with Gitpod. Here are the in-depth instructions
for :ref:`building pandas with GitPod <contributing-gitpod>`.

Step 3: build and install pandas
--------------------------------

There are currently two supported ways of building pandas, pip/meson and setuptools(setup.py).
Historically, pandas has only supported using setuptools to build pandas. However, this method
requires a lot of convoluted code in setup.py and also has many issues in compiling pandas in parallel
due to limitations in setuptools.

The newer build system, invokes the meson backend through pip (via a `PEP 517 <https://peps.python.org/pep-0517/>`_ build).
It automatically uses all available cores on your CPU, and also avoids the need for manual rebuilds by
rebuilding automatically whenever pandas is imported (with an editable install).

For these reasons, you should compile pandas with meson.
Because the meson build system is newer, you may find bugs/minor issues as it matures. You can report these bugs
`here <https://github.com/pandas-dev/pandas/issues/49683>`_.

To compile pandas with meson, run::

   # Build and install pandas
   # By default, this will print verbose output
   # showing the "rebuild" taking place on import (see section below for explanation)
   # If you do not want to see this, omit everything after --no-build-isolation
   python -m pip install -ve . --no-build-isolation --config-settings editable-verbose=true

.. note::
   The version number is pulled from the latest repository tag. Be sure to fetch the latest tags from upstream
   before building::

      # set the upstream repository, if not done already, and fetch the latest tags
      git remote add upstream https://github.com/pandas-dev/pandas.git
      git fetch upstream --tags

**Build options**

It is possible to pass options from the pip frontend to the meson backend if you would like to configure your
install. Occasionally, you'll want to use this to adjust the build directory, and/or toggle debug/optimization levels.

You can pass a build directory to pandas by appending ``--config-settings builddir="your builddir here"`` to your pip command.
This option allows you to configure where meson stores your built C extensions, and allows for fast rebuilds.

Sometimes, it might be useful to compile pandas with debugging symbols, when debugging C extensions.
Appending ``--config-settings setup-args="-Ddebug=true"`` will do the trick.

With pip, it is possible to chain together multiple config settings (for example specifying both a build directory
and building with debug symbols would look like
``--config-settings builddir="your builddir here" --config-settings=setup-args="-Dbuildtype=debug"``.

**Compiling pandas with setup.py**

.. note::
   This method of compiling pandas will be deprecated and removed very soon, as the meson backend matures.

To compile pandas with setuptools, run::

   python setup.py develop

.. note::
   If pandas is already installed (via meson), you have to uninstall it first::

        python -m pip uninstall pandas

This is because python setup.py develop will not uninstall the loader script that ``meson-python``
uses to import the extension from the build folder, which may cause errors such as an
``FileNotFoundError`` to be raised.

.. note::
   You will need to repeat this step each time the C extensions change, for example
   if you modified any file in ``pandas/_libs`` or if you did a fetch and merge from ``upstream/main``.

At this point you should be able to import pandas from your locally built version::

   $ python
   >>> import pandas
   >>> print(pandas.__version__)  # note: the exact output may differ
   2.0.0.dev0+880.g2b9e661fbb.dirty

When building pandas with meson, importing pandas will automatically trigger a rebuild, even when C/Cython files are modified.
By default, no output will be produced by this rebuild (the import will just take longer). If you would like to see meson's
output when importing pandas, you can set the environment variable ``MESONPY_EDTIABLE_VERBOSE``. For example, this would be::

   # On Linux/macOS
   MESONPY_EDITABLE_VERBOSE=1 python

   # Windows
   set MESONPY_EDITABLE_VERBOSE=1 # Only need to set this once per session
   python

If you would like to see this verbose output every time, you can set the ``editable-verbose`` config setting to ``true`` like so::

   python -m pip install -ve . --config-settings editable-verbose=true

.. tip::
   If you ever find yourself wondering whether setuptools or meson was used to build your pandas,
   you can check the value of ``pandas._built_with_meson``, which will be true if meson was used
   to compile pandas.
