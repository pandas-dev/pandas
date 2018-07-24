#!/bin/bash

echo "installing dbs"
# mysql -e 'create database pandas_nosetest;'
psql -d $POSTGRESQL_URL -c 'create database pandas_nosetest;' -U postgres || exit 1

echo "done"
exit 0
