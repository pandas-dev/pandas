#!/bin/bash

# There are 2 distinct pieces that get zipped and cached
# - The venv site-packages dir including the installed dependencies
# - The pandas build artifacts, using the build cache support via
#   scripts/use_build_cache.py
#
# if the user opted in to use the cache and we're on a whitelisted fork
# - if the server doesn't hold a cached version of venv/pandas build,
#   do things the slow way, and put the results on the cache server
#   for the next time.
# -  if the cache files are available, instal some necessaries via apt
#    (no compiling needed), then directly goto script and collect 200$.
#

echo "inside $0"
# Install Dependencie
SITE_PKG_DIR=$VIRTUAL_ENV/lib/python$TRAVIS_PYTHON_VERSION/site-packages
echo "Using SITE_PKG_DIR: $SITE_PKG_DIR"

# workaround for travis ignoring system_site_packages in travis.yml
rm -f $VIRTUAL_ENV/lib/python$TRAVIS_PYTHON_VERSION/no-global-site-packages.txt

if [ x"$LOCALE_OVERRIDE" != x"" ]; then
    # make sure the locale is available
    # probably useless, since you would need to relogin
    sudo locale-gen "$LOCALE_OVERRIDE"
fi;

#scipy is not included in the cached venv
if [ x"$FULL_DEPS" == x"true" ] ; then
   # for pytables gets the lib as well
   sudo apt-get $APT_ARGS install libhdf5-serial-dev

   if [ ${TRAVIS_PYTHON_VERSION:0:1} == "3" ]; then
     sudo apt-get $APT_ARGS install python3-bs4
   elif [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
     sudo apt-get $APT_ARGS install python-bs4
   fi

   if [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
       sudo apt-get $APT_ARGS install python3-scipy
   elif [ ${TRAVIS_PYTHON_VERSION} == "2.7" ]; then
       sudo apt-get $APT_ARGS install python-scipy
   fi
fi

# Everything installed inside this clause into site-packages
# will get included in the cached venv downloaded from the net
# in PTF mode
if ( ! $VENV_FILE_AVAILABLE ); then
    echo "Running full monty"
    # Hard Deps
    pip install $PIP_ARGS  nose python-dateutil pytz
    pip install $PIP_ARGS  cython

    if [ ${TRAVIS_PYTHON_VERSION} == "3.3" ]; then # should be >=3,3
        pip install $PIP_ARGS numpy==1.7.0
    elif [ ${TRAVIS_PYTHON_VERSION} == "3.2" ]; then
        # sudo apt-get $APT_ARGS install python3-numpy; # 1.6.2 or precise
        pip install $PIP_ARGS numpy==1.6.1
    else
        pip install $PIP_ARGS numpy==1.6.1
    fi

    # Optional Deps
    if [ x"$FULL_DEPS" == x"true" ]; then
        echo "Installing FULL_DEPS"
        pip install $PIP_ARGS  cython
        pip install $PIP_ARGS numexpr

        if [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
            pip install $PIP_ARGS xlwt
            pip install $PIP_ARGS bottleneck
        fi

        pip install $PIP_ARGS tables
        pip install $PIP_ARGS matplotlib
        pip install $PIP_ARGS openpyxl
        pip install $PIP_ARGS xlrd>=0.9.0
        pip install $PIP_ARGS 'http://downloads.sourceforge.net/project/pytseries/scikits.timeseries/0.91.3/scikits.timeseries-0.91.3.tar.gz?r='
        pip install $PIP_ARGS patsy
        pip install $PIP_ARGS html5lib

        if [ ${TRAVIS_PYTHON_VERSION:0:1} == "3" ]; then
          sudo apt-get $APT_ARGS remove python3-lxml
        elif [ ${TRAVIS_PYTHON_VERSION:0:1} == "2" ]; then
          sudo apt-get $APT_ARGS remove python-lxml
        fi

        pip install $PIP_ARGS lxml
        # fool statsmodels into thinking pandas was already installed
        # so it won't refuse to install itself. We want it in the zipped venv

        mkdir  $SITE_PKG_DIR/pandas
        touch $SITE_PKG_DIR/pandas/__init__.py
        echo "version='0.10.0-phony'" >  $SITE_PKG_DIR/pandas/version.py
        pip install $PIP_ARGS git+git://github.com/statsmodels/statsmodels@c9062e43b8a5f7385537ca95#egg=statsmodels

        rm -Rf $SITE_PKG_DIR/pandas # scrub phoney pandas
    fi

    # pack up the venv and cache it
    if [ x"$STORE_KEY" != x"" ] && $UPLOAD && $PLEASE_TRAVIS_FASTER ; then
        VENV_FNAME="venv-$TRAVIS_PYTHON_VERSION.zip"

        zip $ZIP_FLAGS -r "$HOME/$VENV_FNAME" $SITE_PKG_DIR/
        ls -l "$HOME/$VENV_FNAME"
        echo "posting venv"
        # silent, don't expose key
        curl -s --form upload=@"$HOME/$VENV_FNAME" "$CACHE_FILE_STORE_URL/$VENV_FNAME"
    fi

fi;

#build and install pandas
if [ x"$BUILD_CACHE_DIR" != x"" ]; then
    scripts/use_build_cache.py -d
    python setup.py install;
else
   python setup.py build_ext install
fi

# package pandas build artifacts and send them home
# that's everything the build cache (scripts/use_build_cache.py)
# stored during the build (.so, pyx->.c and 2to3)
if  (! $CACHE_FILE_AVAILABLE) ; then
 if [ x"$STORE_KEY" != x"" ] && $UPLOAD && $PLEASE_TRAVIS_FASTER ; then
     echo "Posting artifacts"
     strip  "$BUILD_CACHE_DIR/*" &> /dev/null
     echo  "$BUILD_CACHE_DIR"
     cd "$BUILD_CACHE_DIR"/
     zip -r $ZIP_FLAGS "$HOME/$CYTHON_HASH".zip *
     cd "$TRAVIS_BUILD_DIR"
     pwd
     zip "$HOME/$CYTHON_HASH".zip $(find pandas | grep -P '\.(pyx|pxd)$' | sed -r 's/.(pyx|pxd)$/.c/')

     # silent, don't expose key
     curl   --connect-timeout 5 -s --form upload=@"$HOME/$CYTHON_HASH".zip "$CACHE_FILE_STORE_URL/$CYTHON_HASH.zip"
 fi
fi

true
