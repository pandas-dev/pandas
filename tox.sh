#!/usr/bin/env bash


fast=
parallel=

while getopts fp opt; do
    case $opt in
        f) fast=true;;
        p) parallel=true;;
    esac
done


if [ $fast ]; then
    scripts/use_build_cache.py
fi

if [ $parallel ]; then
    ENVS=$(cat tox.ini | grep envlist | tr  "," " " | cut -d " " -f 3-)
    TOX_INI_PAR="tox.ini"
    echo "[Creating distfile]"
    tox --sdistonly
    export DISTFILE="$(find .tox/dist -type f)"

    echo -e "[Starting tests]\n"
    for e in $ENVS; do
        echo "[launching tox for $e]"
        tox -c "$TOX_INI_PAR" -e "$e" &
    done
    wait
else
    tox
fi
