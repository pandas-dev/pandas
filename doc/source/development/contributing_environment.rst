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

.. contents:: Table of contents:
   :local:


Option 1: creating an environment without Docker
------------------------------------------------

Installing a C compiler
~~~~~~~~~~~~~~~~~~~~~~~

pandas uses C extensions (mostly written using Cython) to speed up certain
operations. To install pandas from source, you need to compile these C
extensions, which means you need a C compiler. This process depends on which
platform you're using.

If you have setup your environment using :ref:`mamba <contributing.mamba>`, the packages ``c-compiler``
and ``cxx-compiler`` will install a fitting compiler for your platform that is
compatible with the remaining mamba packages. On Windows and macOS, you will
also need to install the SDKs as they have to be distributed separately.
These packages will automatically be installed by using the ``pandas``
``environment.yml`` file.

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
as an installed compiler nothing more is required. If no C compiler is
installed (or you wish to install a newer version) you can install a compiler
(GCC in the example code below) with::

    # for recent Debian/Ubuntu:
    sudo apt install build-essential
    # for Red Had/RHEL/CentOS/Fedora
    yum groupinstall "Development Tools"

For other Linux distributions, consult your favorite search engine for
compiler installation instructions.

Let us know if you have any difficulties by opening an issue or reaching out on our contributor
community :ref:`Slack <community.slack>`.

.. _contributing.mamba:

Option 1a: using mamba (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now create an isolated pandas development environment:

* Install `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_
* Make sure your mamba is up to date (``mamba update mamba``)
* Make sure that you have :any:`cloned the repository <contributing.forking>`
* ``cd`` to the pandas source directory

We'll now kick off a three-step process:

1. Install the build dependencies
2. Build and install pandas
3. Install the optional dependencies

.. code-block:: none

   # Create and activate the build environment
   mamba env create --file environment.yml
   mamba activate pandas-dev

   # Build and install pandas
   python setup.py build_ext -j 4
   python -m pip install -e . --no-build-isolation --no-use-pep517

At this point you should be able to import pandas from your locally built version::

   $ python
   >>> import pandas
   >>> print(pandas.__version__)  # note: the exact output may differ
   1.5.0.dev0+1355.ge65a30e3eb.dirty

This will create the new environment, and not touch any of your existing environments,
nor any existing Python installation.

To return to your root environment::

      mamba deactivate

Option 1b: using pip
~~~~~~~~~~~~~~~~~~~~

If you aren't using mamba for your development environment, follow these instructions.
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

   # Build and install pandas
   python setup.py build_ext -j 4
   python -m pip install -e . --no-build-isolation --no-use-pep517

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

   # Build and install pandas
   python setup.py build_ext -j 4
   python -m pip install -e . --no-build-isolation --no-use-pep517

**Windows**

Below is a brief overview on how to set-up a virtual environment with Powershell
under Windows. For details please refer to the
`official virtualenv user guide <https://virtualenv.pypa.io/en/latest/user_guide.html#activators>`__.

Use an ENV_DIR of your choice. We'll use ~\\virtualenvs\\pandas-dev where
'~' is the folder pointed to by either $env:USERPROFILE (Powershell) or
%USERPROFILE% (cmd.exe) environment variable. Any parent directories
should already exist.

.. code-block:: powershell

   # Create a virtual environment
   python -m venv $env:USERPROFILE\virtualenvs\pandas-dev

   # Activate the virtualenv. Use activate.bat for cmd.exe
   ~\virtualenvs\pandas-dev\Scripts\Activate.ps1

   # Install the build dependencies
   python -m pip install -r requirements-dev.txt

   # Build and install pandas
   python setup.py build_ext -j 4
   python -m pip install -e . --no-build-isolation --no-use-pep517

Option 2: creating an environment using Docker
----------------------------------------------

Instead of manually setting up a development environment, you can use `Docker
<https://docs.docker.com/get-docker/>`_ to automatically create the environment with just several
commands. pandas provides a ``DockerFile`` in the root directory to build a Docker image
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

When inside the running container you can build and install pandas the same way as the other methods

.. code-block:: bash

   python setup.py build_ext -j 4
   python -m pip install -e . --no-build-isolation --no-use-pep517

*Even easier, you can integrate Docker with the following IDEs:*

**Visual Studio Code**

You can use the DockerFile to launch a remote session with Visual Studio Code,
a popular free IDE, using the ``.devcontainer.json`` file.
See https://code.visualstudio.com/docs/remote/containers for details.

**PyCharm (Professional)**

Enable Docker support and use the Services tool window to build and manage images as well as
run and interact with containers.
See https://www.jetbrains.com/help/pycharm/docker.html for details.
