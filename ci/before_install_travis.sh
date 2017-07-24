#!/bin/bash

echo "inside $0"

if [ "${TRAVIS_OS_NAME}" == "linux" ]; then
   sh -e /etc/init.d/xvfb start
fi

# Never fail because bad things happened here.
true
