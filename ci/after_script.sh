#!/bin/bash

#wget https://raw.github.com/y-p/ScatterCI-CLI/master/scatter_cli.py
#chmod u+x scatter_cli.py

pip install -I requests==2.1.0
echo "${TRAVIS_PYTHON_VERSION:0:4}"
if [ x"${TRAVIS_PYTHON_VERSION:0:4}" == x"2.6" ]; then
    pip install simplejson;
fi

# ScatterCI accepts a build log, but currently does nothing with it.
echo '' > /tmp/build.log

# nore exposed in the build logs
#export SCATTERCI_ACCESS_KEY=
#export SCATTERCI_HOST=

# Generate a json file describing system and dep versions
ci/print_versions.py -j /tmp/env.json

# nose ran using "--with-xunit --xunit-file nosetest.xml" and generated /tmp/nosetest.xml
# Will timeout if server not available, and should not fail the build
#python scatter_cli.py --xunit-file /tmp/nosetests.xml  --log-file /tmp/build.log  --env-file /tmp/env.json --build-name "$JOB_NAME" --succeed

true # never fail because bad things happened here
