#!/bin/bash
# Add 3rd party licenses, like numpy does
for file in $PACKAGE_DIR/LICENSES/*; do
  cat $file >> $PACKAGE_DIR/LICENSE
done

# TODO: Delete when there's a PyPI Cython release that supports free-threaded Python 3.13
FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    python -m pip install -U pip
    # python -m pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple cython
    # TODO: Remove below and uncomment above once https://github.com/cython/cython/pull/6717 no longer breaks tests
    python -m pip install git+https://github.com/cython/cython.git@3276b588720a053c78488e5de788605950f4b136
    python -m pip install ninja meson-python versioneer[toml] numpy
fi
