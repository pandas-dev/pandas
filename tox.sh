#!/usr/bin/env bash


use_build_cache=
slow_tests=
full_deps=

while getopts bsd opt; do
    case $opt in
        b) use_build_cache=true;;
        s) slow_tests=true;;
        d) full_deps=true;;
    esac
done


if [ $use_build_cache ]; then
    scripts/use_build_cache.py
fi


# choose to run all tests or just to run not network and not slow (default)
posargs="-A 'not network and not slow'"
if [ $slow_tests ]; then
    posargs=
fi


# choose full or slim build deps
start_i=3
end_i=9
if [ $full_deps ]; then
    start_i=9
    end_i=
fi

ENVS=$(cat tox.ini | grep envlist | tr ',' ' ' | cut -d " " -f ${start_i}-${end_i})


TOX_INI_PAR="tox.ini"
echo "[Creating distfile]"
tox --sdistonly
export DISTFILE="$(find .tox/dist -type f)"


# run the tests
echo -e "[Starting tests]\n"
for e in "${ENVS}"; do
    echo "[launching tox for $e]"
    tox -c "$TOX_INI_PAR" -e "$e" -- "${posargs}" &
done
wait
