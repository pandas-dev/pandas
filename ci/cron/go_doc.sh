#!/bin/bash

# This is a one-command cron job for setting up
# a virtualenv-based, linux-based, py2-based environment
# for building the Pandas documentation.
#
# The first run will install all required deps from pypi
# into the venv including monsters like scipy.
# You may want to set it up yourself to speed up the
# process.
#
# This is meant to be run as a cron job under a dedicated
# user account whose HOME directory contains this script.
# a CI directory will be created under it and all files
# stored within it.
#
# The hardcoded dep versions will gradually become obsolete
# You may need to tweak them
#
# @y-p, Jan/2014

# disto latex is sometimes finicky. Optionall use
# a local texlive install
export PATH=/mnt/debian/texlive/2013/bin/x86_64-linux:$PATH

# Having ccache will speed things up
export PATH=/usr/lib64/ccache/:$PATH

# limit disk usage
ccache -M 200M

BASEDIR="$HOME/CI"
REPO_URL="https://github.com/pydata/pandas"
REPO_LOC="$BASEDIR/pandas"

if [ ! -d $BASEDIR ]; then
	mkdir -p $BASEDIR
	virtualenv $BASEDIR/venv
fi

source $BASEDIR/venv/bin/activate

pip install numpy==1.7.2
pip install cython==0.20.0
pip install python-dateutil==2.2
pip install --pre pytz==2013.9
pip install sphinx==1.1.3
pip install numexpr==2.2.2

pip install matplotlib==1.3.0
pip install lxml==3.2.5
pip install beautifulsoup4==4.3.2
pip install html5lib==0.99

# You'll need R as well
pip install rpy2==2.3.9

pip install tables==3.0.0
pip install bottleneck==0.7.0
pip install ipython==0.13.2

# only if you have too
pip install scipy==0.13.2

pip install openpyxl==1.6.2
pip install xlrd==0.9.2
pip install xlwt==0.7.5
pip install xlsxwriter==0.5.1
pip install sqlalchemy==0.8.3

if [ ! -d "$REPO_LOC" ]; then
	git clone "$REPO_URL" "$REPO_LOC"
fi

cd "$REPO_LOC"
git reset --hard
git clean -df
git checkout master
git pull origin
make

source $BASEDIR/venv/bin/activate
export PATH="/usr/lib64/ccache/:$PATH"
pip uninstall pandas -yq
pip install "$REPO_LOC"

cd "$REPO_LOC"/doc

python make.py clean
python make.py html
if [ ! $? == 0 ]; then
	exit 1
fi
python make.py zip_html
# usually requires manual intervention
# python make.py latex

# If you have access:
# python make.py upload_dev
