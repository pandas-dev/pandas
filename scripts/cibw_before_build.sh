# TODO: Delete when there's PyPI NumPy/Cython releases the support Python 3.13.
# If free-threading support is not included in those releases, this script will have
# to whether this runs for a free-threaded build instead.
PYTHON_VERSION="$(python -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')")"
if [[ $PYTHON_VERSION == "313" ]]; then
    python -m pip install -U pip
    python -m pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy cython
    python -m pip install ninja meson-python versioneer[toml]
fi
