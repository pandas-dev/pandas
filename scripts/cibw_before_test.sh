# TODO: delete along with enabling build isolation by unsetting
# CIBW_BUILD_FRONTEND when scipy is buildable under free-threaded
# python with a released version of cython
FREE_THREADED_BUILD="$(python -c "import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    python -m pip install -U pip
    python -m pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy
fi
