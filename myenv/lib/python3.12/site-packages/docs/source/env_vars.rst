ASV environment variables
=========================

Benchmarking and build commands are run with the following environment
variables available:

- ``ASV``: ``true``
- ``ASV_PROJECT``: the project name from the configuration file
- ``ASV_ENV_NAME``: name of the currently active environment
- ``ASV_ENV_TYPE``: type of the currently active environment
- ``ASV_ENV_DIR``: full path to the currently active environment root
- ``ASV_CONF_DIR``: full path to the directory where ``asv.conf.json`` is
- ``ASV_BUILD_DIR``: full path to the build directory (checked-out source path + ``repo_subdir``)
- ``ASV_BUILD_CACHE_DIR``: full path to the build cache directory
- ``ASV_COMMIT``: commit hash of currently installed project

If there is no asv-managed environment, build, or cache directory, or
commit hash, those environment variables are unset.

The following environment variables controlling Python and other
behavior are also set:

- ``PATH``: environment-specific binary directories prepended
- ``PIP_USER``: ``false``
- ``PYTHONNOUSERSITE``: ``True`` (for conda environments only)
- ``PYTHONPATH``: unset (if really needed, can be overridden by setting ``ASV_PYTHONPATH``)

.. note::

    .. versionadded::0.6.0

      ``ASV_RUNNER_PATH`` may be set to provide a local installation of
      ``asv_runner``, mostly used for the CI to ensure changes to ``asv_runner``
      do not break ``asv``


Custom environment variables
----------------------------

You can send custom environment variables to build and benchmarking commands
by configuring the ``matrix`` setting in ``asv.conf.json``.
