#!/bin/bash -e

if [[ $PATTERN == "manylinux2014_aarch64" ]]; then
    docker pull multiarch/qemu-user-static
    docker run --privileged multiarch/qemu-user-static --reset -p yes
fi

docker pull quay.io/pypa/$IMAGE
docker run -v $(pwd):/pandas quay.io/pypa/$IMAGE \
/bin/bash -xc "cd pandas && \
/opt/python/cp37-cp37m/bin/python -m venv ~/virtualenvs/pandas-dev && \
. ~/virtualenvs/pandas-dev/bin/activate && \
python -m pip install --no-deps -U pip wheel setuptools && \
pip install cython numpy python-dateutil pytz pytest pytest-xdist hypothesis pytest-azurepipelines && \
python setup.py build_ext -q -j2 && \
python -m pip install --no-build-isolation -e . && \
pytest -m \"$PATTERN\" pandas --junitxml=test-data.xml"
