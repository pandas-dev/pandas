.. _conf-reference:

``asv.conf.json`` reference
===========================

The ``asv.conf.json`` file contains information about a particular
benchmarking project.  The following describes each of the keys in
this file and their expected values.

.. note::

    The GitHub repository at `asv_samples
    <https://github.com/airspeed-velocity/asv_samples>`_ serves as a
    comprehensive showcase for integrating Air Speed Velocity (ASV) with a wide
    array of Python project configurations. This includes various build systems
    and advanced benchmarking features like custom parameterizations and ASV
    plugins, aiming to benchmark Python code performance across diverse setups.

    The repository is structured with dedicated branches for each build system
    and feature demonstration, providing insights into the impacts of different
    build systems and ASV's extensible features on performance metrics.

.. only:: not man

   .. contents::

``project``
-----------
The name of the project being benchmarked.

``project_url``
---------------
The URL to the homepage of the project.  This can point to anywhere,
really, as it's only used for the link at the top of the benchmark
results page back to your project.

``repo``
--------
The URL to the repository for the project.

The value can also be a path, relative to the location of the
configuration file. For example, if the benchmarks are stored in the
same repository as the project itself, and the configuration file is
located at ``benchmarks/asv.conf.json`` inside the repository, you can
set ``"repo": ".."`` to use the local repository.

Currently, only ``git`` and ``hg`` repositories are supported, so this must be
a URL that ``git`` or ``hg`` know how to clone from, for example:

   - git@github.com:airspeed-velocity/asv.git

   - https://github.com/airspeed-velocity/asv.git

   - ssh://hg@bitbucket.org/yt_analysis/yt

   - hg+https://bitbucket.org/yt_analysis/yt

The repository may be readonly.

``repo_subdir``
---------------

The relative path to your Python project inside the repository.  This is
where its ``setup.py`` file is located.

If empty or omitted, the project is assumed to be located at the root of
the repository.


``build_command``, ``install_command``, ``uninstall_command``
-------------------------------------------------------------

Airspeed Velocity rebuilds the project as needed, using these commands.

The defaults are::

  "install_command":
  ["in-dir={env_dir} python -mpip install {wheel_file}"],

  "uninstall_command":
  ["return-code=any python -mpip uninstall -y {project}"],

  "build_command":
  ["python setup.py build",
   "PIP_NO_BUILD_ISOLATION=false python -mpip wheel --no-deps --no-index -w {build_cache_dir} {build_dir}"],

.. note::

    .. versionchanged:: 0.6.2

        The default build command now assume network connectivity is not
        prohibited. The ``build_command`` is now::

          "build_command":
          ["python setup.py build",
           "python -mpip wheel -w {build_cache_dir} {build_dir}"],

The install commands should install the project in the active Python
environment (virtualenv/conda), so that it can be used by the
benchmark code.

The uninstall commands should uninstall the project from the
environment.

.. note::

    .. versionchanged:: 0.6.0

        If a build command is not specified in the ``asv.conf.json``, the default
        assumes the build system requirements are defined in a ``setup.py`` file.
        ``pyproject.toml`` is the preferred  file format to define the build  system
        requirements of Python projects (`PEP518
        <https://peps.python.org/pep-0518/>`_), and this approach will be the
        default from ``asv v0.6.0`` onwards.

The build commands can optionally be used to cache build results in the
cache directory ``{build_cache_dir}``, which is commit and
environment-specific.  If the cache directory contains any files after
``build_command`` finishes with exit code 0, ``asv`` assumes it
contains a cached build.  When a cached build is available, ``asv``
will only call ``install_command`` but not ``build_command``. (The
number of cached builds retained at any time is determined by the
``build_cache_size`` configuration option.)

The ``install_command`` and ``build_command`` are by default launched
in ``{build_dir}``. The ``uninstall_command`` is launched in the
environment root directory.

The commands are specified in typical POSIX shell syntax (Python
shlex), but are not run in a shell, so that e.g. ``cd`` has no effect
on subsequent commands, and wildcard or environment variable
expansion is not done. The substituted variables ``{variable_name}``
do not need to be quoted. The commands may contain environment
variable specifications in in form ``VARNAME=value`` at the beginning.
In addition, valid return codes can be specified via
``return-code=0,1,2`` and ``return-code=any``.

The ``in-dir=somedir`` specification changes the working directory
for the command.

The commands can be supplied with the arguments:

- ``{project}``: the project name from the configuration file
- ``{env_name}``: name of the currently active environment
- ``{env_type}``: type of the currently active environment
- ``{env_dir}``: full path to the currently active environment root
- ``{conf_dir}``: full path to the directory where ``asv.conf.json`` is
- ``{build_dir}``: full path to the build directory (checked-out source path + ``repo_subdir``)
- ``{build_cache_dir}``: full path to the build cache directory
- ``{commit}``: commit hash of currently installed project
- ``{wheel_file}``: absolute path to a ``*.whl`` file in ``{build_cache_dir}``
  (defined only if there is exactly one existing wheel file in the directory).

Several :doc:`environment variables <env_vars>` are also defined.


``branches``
------------
Branches to generate benchmark results for.

This controls how the benchmark results are displayed, and what
benchmarks ``asv run ALL`` and ``asv run NEW`` run.

If not provided, "main" (Git) or "default" (Mercurial) is chosen.

``show_commit_url``
-------------------
The base URL to show information about a particular commit.  The
commit hash will be added to the end of this URL and then opened in a
new tab when a data point is clicked on in the web interface.

For example, if using Github to host your repository, the
``show_commit_url`` should be:

    http://github.com/owner/project/commit/

``pythons``
-----------
The versions of Python to run the benchmarks in.  If not provided, it
will to default to the version of Python that the ``asv`` command
(main) is being run under.

If provided, it should be a list of strings.  It may be one of the
following:

- a Python version string, e.g. ``"3.7"``, in which case:

  - if ``conda`` is found, ``conda`` will be used to create an
    environment for that version of Python via a temporary
    environment.yml file

  - if ``virtualenv`` is installed, ``asv`` will search for that
    version of Python on the ``PATH`` and create a new virtual
    environment based on it.  ``asv`` does not handle downloading and
    installing different versions of Python for you.  They must
    already be installed and on the path.  Depending on your platform,
    you can install multiple versions of Python using your package
    manager or using `pyenv <https://github.com/yyuu/pyenv>`_.

- an executable name on the ``PATH`` or an absolute path to an
  executable.  In this case, the environment is assumed to be already
  fully loaded and read-only.  Thus, the benchmarked project must
  already be installed, and it will not be possible to benchmark
  multiple revisions of the project.

``conda_environment_file``
--------------------------
A path to a ``conda`` environment file to use as source for the
dependencies. For example::

    "conda_environment_file": "environment.yml"

The environment file should generally install ``wheel`` and ``pip``,
since those are required by the default ``asv`` build commands.  If there
are packages present in ``matrix``, an additional ``conda env update``
call is used to install them after the environment is created.

.. note::

   .. versionchanged:: 0.6.0

    If an ``environment.yml`` file is present where
    ``asv`` is run, it will be used. To turn off this behavior,
    ``conda_environment_file`` can be set to ``IGNORE``.

This option will cause ``asv`` to ignore the Python version in the
environment creation, which is then assumed to be fixed by the
environment file.

``conda_channels``
------------------
A list of ``conda`` channel names (strings) to use in the provided
order as the source channels for the dependencies. For example::

    "conda_channels": ["conda-forge", "defaults"]

The channels will be parsed by ``asv`` to populate the ``channels``
section of a temporary environment.yml file used to build the
benchmarking environment.

``matrix``
----------
Defines a matrix of third-party dependencies and environment variables
to run the benchmarks with.

If provided, it must be a dictionary, containing some of the keys
"req", "env", "env_nobuild". For example::

    "matrix": {
        "req": {
            "numpy": ["1.25", "1.26"],
            "Cython": []
            "six": ["", null]
        },
        "env": {
            "FOO": "bar"
        }
    }

The keys of the ``"req"`` are the names of dependencies, and the
values are lists of versions (as strings) of that dependency.  An
empty string means the "latest" version of that dependency available
on PyPI. Value of ``null`` means the package will not be installed.

If the list is empty, it is equivalent to ``[""]``, in other words,
the "latest" version.

For example, the following will test with two different versions of
Numpy, the latest version of Cython, and six installed as the latest
version and not installed at all::

    "matrix": {
        "req": {
            "numpy": ["1.25", "1.26"],
            "Cython": []
            "six": ["", null],
        }
    }

The matrix dependencies are installed *before* any dependencies that
the project being benchmarked may specify in its ``setup.py`` file.

.. note::

    At present, this functionality only supports dependencies that are
    installable via ``pip`` or ``conda`` or ``mamba`` (depending on which
    environment is used). If ``conda/mamba`` is specified as
    ``environment_type`` and you wish to install the package via ``pip``, then
    preface the package name with ``pip+``. For example, ``emcee`` is only
    available from ``pip``, so the package name to be used is ``pip+emcee``.

    .. versionadded::0.6.0

      ``pip`` dependencies can now accept local (fully qualified) directories,
      and also take flags (e.g. ``-e``)

    .. versionadded::0.6.1

       ``asv`` can now optionally load dependencies from ``environment.yml`` if
       ``conda`` or ``mamba`` is set as the ``environment_type``. As ``asv``
       dependencies are explicitly mentioned only in the ``asv.conf.json``.
       These specifications in ``environment.yml`` or another (user-defined)
       file will be overridden by the environment matrix.

    .. versionadded::0.6.2

       The ``mamba`` plugin will now take channels and channel priority from the
       ``MAMBARC`` environment variable if it is provided. e.g.
       ``MAMBARC=$HOME/.condarc asv run``. By default user ``.rc`` files are not
       read to enforce isolation.

The ``env`` and ``env_nobuild`` dictionaries can be used to set also
environment variables::

   "matrix": {
       "env": {
           "ENV_VAR_1": ["val1", "val2"],
           "ENV_VAR_2": ["val3", null],
       },
       "env_nobuild": {
           "ENV_VAR_3": ["val4", "val5"],
       }
   }

Variables in "no_build" will be passed to every environment during the test
phase, but will not trigger a new build.
A value of ``null`` means that the variable will not be set for the current
combination.

The above matrix will result in 4 different builds with the following
additional environment variables and values:

  - [("ENV_VAR_1", "val1"), ("ENV_VAR_2", "val3")]
  - [("ENV_VAR_1", "val1")]
  - [("ENV_VAR_1", "val2"), ("ENV_VAR_2", "val3")]
  - [("ENV_VAR_1", "val2")]

It will generate 8 different test environments based on those 4 builds with
the following environment variables and values:

  - [("ENV_VAR_1", "val1"), ("ENV_VAR_2", "val3"), ("ENV_VAR_3", "val4")]
  - [("ENV_VAR_1", "val1"), ("ENV_VAR_2", "val3"), ("ENV_VAR_3", "val5")]
  - [("ENV_VAR_1", "val1"), ("ENV_VAR_3", "val4")]
  - [("ENV_VAR_1", "val1"), ("ENV_VAR_3", "val5")]
  - [("ENV_VAR_1", "val2"), ("ENV_VAR_2", "val3"), ("ENV_VAR_3", "val4")]
  - [("ENV_VAR_1", "val2"), ("ENV_VAR_2", "val3"), ("ENV_VAR_3", "val5")]
  - [("ENV_VAR_1", "val2"), ("ENV_VAR_3", "val4")]
  - [("ENV_VAR_1", "val2"), ("ENV_VAR_3", "val5")]


``exclude``
-----------
Combinations of libraries, Python versions, or platforms to be
excluded from the combination matrix. If provided, must be a list of
dictionaries, each specifying an exclude rule.

An exclude rule consists of key-value pairs, specifying matching rules
``matrix[key] ~ value``. The values are strings containing regular
expressions that should match whole strings.  The exclude rule matches
if all of the items in it match.

Each exclude rule can contain the following keys:

- ``python``: Python version (from ``pythons``)

- ``sys_platform``: Current platform, as in ``sys.platform``.
  Common values are: ``linux2``, ``win32``, ``cygwin``, ``darwin``.

- ``environment_type``: The environment type in use (from ``environment_type``).

- ``req``: dictionary of rules vs. the requirements

- ``env``: dictionary of rules vs. environment variables

- ``env_nobuild``: : dictionary of rules vs. the non-build environment variables

For example::

    "pythons": ["3.8", "3.9"],
    "matrix": {
        "req": {
            "numpy": ["1.25", "1.26"],
            "Cython": ["", null],
            "colorama": ["", null]
        },
        "env": {"FOO": ["1", "2"]},
    },
    "exclude": [
        {"python": "3.8", "req": {"numpy": "1.25"}},
        {"sys_platform": "(?!win32).*", "req": {"colorama": ""}},
        {"sys_platform": "win32", "req": {"colorama": null}},
        {"env": {"FOO": "1"}},
    ]

This will generate all combinations of Python version and items in the
matrix, except those with Python 3.8 and Numpy 3.9. In other words,
the combinations::

    python==3.8 numpy==1.26 Cython==latest (colorama==latest) FOO=2
    python==3.8 numpy==1.26 (colorama==latest) FOO=2
    python==3.9 numpy==1.25 Cython==latest (colorama==latest) FOO=2
    python==3.9 numpy==1.25 (colorama==latest) FOO=2
    python==3.9 numpy==1.26 Cython==latest (colorama==latest) FOO=2
    python==3.9 numpy==1.26 (colorama==latest) FOO=2

The ``colorama`` package will be installed only if the current
platform is Windows.


``include``
-----------
Additional package combinations to be included as environments.

If specified, must be a list of dictionaries, indicating the versions
of packages and other environment configuration to be installed. The
dictionary must also include a ``python`` key specifying the Python
version.

Similarly as for the matrix, the ``"req"``, ``"env"`` and ``"env_nobuild"``
entries specify dictionaries containing requirements and environment variables.
In contrast to the matrix, the values are not lists, but a single value only.

In addition, the following keys can be present: ``sys_platform``,
``environment_type``.  If present, the include rule is active only if
the values match, using same matching rules as explained for
``exclude`` above.

The exclude rules are not applied to includes.

For example::

    "include": [
        {"python": "3.9", "req": {"numpy": "1.26"}, "env": {"FOO": "true"}},
        {"platform": "win32", "environment_type": "conda",
         "req": {"python": "3.12", "libpython": ""}}
    ]

This corresponds to two additional environments. One runs on Python 3.9
and including the specified version of Numpy. The second is active only
for Conda on Windows, and installs the latest version of ``libpython``.

``benchmark_dir``
-----------------
The directory, relative to the current directory, that benchmarks are
stored in.  Should rarely need to be overridden.  If not provided,
defaults to ``"benchmarks"``.

``environment_type``
--------------------
Specifies the tool to use to create environments.  May be "conda",
"virtualenv", "mamba" or another value depending on the plugins in use.  If
missing or the empty string, the tool will be automatically determined
by looking for tools on the ``PATH`` environment variable.

``env_dir``
-----------
The directory, relative to the current directory, to cache the Python
environments in.  If not provided, defaults to ``"env"``.

``results_dir``
---------------
The directory, relative to the current directory, that the raw results
are stored in.  If not provided, defaults to ``"results"``.

``html_dir``
------------
The directory, relative to the current directory, to save the website
content in.  If not provided, defaults to ``"html"``.

``hash_length``
---------------
The number of characters to retain in the commit hashes when displayed
in the web interface.  The default value of 8 should be more than
enough for most projects, but projects with extremely large history
may need to increase this value.  This does not affect the storage of
results, where the full commit hash is always retained.

``plugins``
-----------
A list of modules to import containing asv plugins.

``build_cache_size``
--------------------
The number of builds to cache for each environment.

``regressions_first_commits``
-----------------------------

The commits after which the regression search in :ref:`cmd-asv-publish`
should start looking for regressions.

The value is a dictionary mapping benchmark identifier regexps to
commits after which to look for regressions. The benchmark identifiers
are of the form ``benchmark_name(parameters)@branch``, where
``(parameters)`` is present only for parameterized benchmarks. If the
commit identifier is *null*, regression detection for the matching
benchmark is skipped.  The default is to start from the first commit
with results.

Example::

    "regressions_first_commits": {
        ".*": "v0.1.0",
        "benchmark_1": "80fca08d",
        "benchmark_2@main": null,
    }

In this case, regressions are detected only for commits after tag
``v0.1.0`` for all benchmarks. For ``benchmark_1``, regression
detection is further limited to commits after the commit given, and
for ``benchmark_2``, regression detection is skipped completely in the
``main`` branch.

``regressions_thresholds``
--------------------------

The minimum relative change required before :ref:`cmd-asv-publish` reports a
regression.

The value is a dictionary, similar to ``regressions_first_commits``.
If multiple entries match, the largest threshold is taken.  If no
entry matches, the default threshold is ``0.05`` (iow. 5%).

Example::

    "regressions_thresholds": {
        ".*": 0.01,
        "benchmark_1": 0.2,
    }

In this case, the reporting threshold is 1% for all benchmarks, except
``benchmark_1`` which uses a threshold of 20%.
