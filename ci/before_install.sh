#!/bin/bash

echo "inside $0"

# overview
if [ ${TRAVIS_PYTHON_VERSION} == "3.3" ]; then
    sudo add-apt-repository -y ppa:doko/ppa # we get the py3.3 debs from here
fi

sudo apt-get update $APT_ARGS # run apt-get update for all versions

# # hack for broken 3.3 env
# if [ x"$VIRTUAL_ENV" == x"" ]; then
#     VIRTUAL_ENV=~/virtualenv/python$TRAVIS_PYTHON_VERSION_with_system_site_packages;
# fi

# # we only recreate the virtualenv for 3.x
# # since the "Detach bug" only affects python3
# # and travis has numpy preinstalled on 2.x which is quicker
# _VENV=$VIRTUAL_ENV # save it
# if [ ${TRAVIS_PYTHON_VERSION:0:1} == "3" ] ; then
#     deactivate # pop out of any venv
#     sudo pip install virtualenv==1.8.4 --upgrade
#     sudo apt-get install $APT_ARGS python3.3 python3.3-dev
#     sudo rm -Rf $_VENV
#     virtualenv -p python$TRAVIS_PYTHON_VERSION $_VENV --system-site-packages;
#     source $_VENV/bin/activate
# fi
