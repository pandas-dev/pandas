#!/usr/bin/env bash
#
# tox has an undocumented (as of 1.4.2) config option called "sdistsrc"
# which can make a run use a pre-prepared sdist file.
# we prepare the sdist once , then launch the tox runs in parallel using it.
#
# currently (tox 1.4.2) We have to skip sdist generation when running in parallel
# or we get a race.
#


ENVS=$(cat tox.ini | grep envlist | tr  "," " " | cut -d " " -f 3-)
TOX_INI_PAR="tox_prll.ini"

echo "[Creating distfile]"
tox --sdistonly
export DISTFILE="$(find .tox/dist -type f )"

echo -e "[Starting tests]\n"
for e in $ENVS; do
    echo "[launching tox for $e]"
    tox -c "$TOX_INI_PAR" -e "$e" &
done
