# TODO: Delete when there's a NumPy Windows wheel for the free-threaded build on PyPI.
FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    python -m pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy
fi
