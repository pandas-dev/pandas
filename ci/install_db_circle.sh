#!/bin/bash

echo "installing dbs"
mysql -e 'create database pandas_nosetest;'
psql -c 'create database pandas_nosetest;' -U postgres

echo "done"
exit 0
