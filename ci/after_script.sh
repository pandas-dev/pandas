#!/bin/bash

wget https://raw.github.com/y-p/ScatterCI-CLI/master/scatter_cli.py
chmod u+x scatter_cli.py
echo '' > /tmp/build.log
pip install -I requests==2.1.0
echo "${TRAVIS_PYTHON_VERSION:0:4}"
if [ x"${TRAVIS_PYTHON_VERSION:0:4}" == x"2.6" ]; then
    pip install simplejson;
fi

python scatter_cli.py --xunit-file /tmp/nosetests.xml  --log-file /tmp/build.log  --env-file /tmp/env.json --succeed

true # never fail because bad things happened here
