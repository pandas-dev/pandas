#!/bin/bash

if [ "${TRAVIS_OS_NAME}" != "linux" ]; then
   echo "not using dbs on non-linux"
   exit 0
fi

echo "installing dbs"
mysql -e 'create database pandas_nosetest;'
psql -c 'create database pandas_nosetest;' -U postgres

echo "done"
exit 0
