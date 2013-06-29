#!/bin/bash

# This script is meant to run on a mint precise64 VM.
# The generated wheel files should be compatible
# with travis-ci as of 07/2013.
#
# Runtime can be up to an hour or more.

echo "Running build.sh..."
set -x

WHEEL_DIR=/wheelhouse
VERSIONS="2.6 2.7 3.2 3.3"
SCRIPT_FILE="/tmp/run.sh"
PARALLEL=false

export PIP_ARGS=" --download-cache /tmp -w $WHEEL_DIR --use-wheel --find-links=$WHEEL_DIR"

apt-get update
apt-get install python-software-properties git -y
apt-add-repository ppa:fkrull/deadsnakes -y
apt-get update

apt-get install python-pip libfreetype6-dev libpng12-dev -y
pip install virtualenv
apt-get install libhdf5-serial-dev g++ -y


function generate_wheels {
	VER=$1
	set -x

	if [ x"$VIRTUAL_ENV" != x"" ]; then
		deactivate
	fi

	cd ~/
	sudo rm -Rf venv-$VER
	virtualenv -p python$VER venv-$VER
	source venv-$VER/bin/activate

	pip install -I --download-cache /tmp git+https://github.com/pypa/pip@42102e9d#egg=pip
	pip install -I --download-cache  /tmp https://bitbucket.org/pypa/setuptools/downloads/setuptools-0.8b6.tar.gz
	pip install -I --download-cache /tmp wheel

	export INCLUDE_PATH=/usr/include/python$VER/
	export C_INCLUDE_PATH=/usr/include/python$VER/
	pip wheel $PIP_ARGS cython==0.19.1
	pip install --use-wheel --find-links=$WHEEL_DIR cython==0.19.1

	pip wheel $PIP_ARGS numpy==1.6.1
	pip wheel $PIP_ARGS numpy==1.7.1
	pip install --use-wheel --find-links=$WHEEL_DIR numpy==1.7.1
	pip wheel $PIP_ARGS bottleneck==0.6.0

	pip wheel $PIP_ARGS numexpr==1.4.2
	pip install --use-wheel --find-links=$WHEEL_DIR numexpr==1.4.2
	pip wheel $PIP_ARGS tables==2.3.1
	pip wheel $PIP_ARGS tables==2.4.0

	pip uninstall numexpr -y
	pip wheel $PIP_ARGS numexpr==2.1
	pip install --use-wheel --find-links=$WHEEL_DIR numexpr==2.1
	pip wheel $PIP_ARGS tables==3.0.0
	pip uninstall numexpr -y

	pip wheel $PIP_ARGS matplotlib==1.2.1
}


for VER in $VERSIONS ; do
	apt-get install python$VER python$VER-dev -y
done

if $PARALLEL; then
	echo '#!/bin/bash' > $SCRIPT_FILE
	echo "export WHEEL_DIR=$WHEEL_DIR" >> $SCRIPT_FILE
	echo "export PIP_ARGS='$PIP_ARGS'">> $SCRIPT_FILE

	declare -f generate_wheels >>  $SCRIPT_FILE
	echo 'generate_wheels $1' >> $SCRIPT_FILE
	chmod u+x $SCRIPT_FILE

	pip install -I --download-cache /tmp git+https://github.com/pypa/pip@42102e9d#egg=pip
	pip install --download-cache /tmp --no-install wheel
	pip install --download-cache /tmp --no-install https://bitbucket.org/pypa/setuptools/downloads/setuptools-0.8b6.tar.gz

	for VER in 2.6 2.7 3.2 3.3 ; do
		$SCRIPT_FILE  $VER &
	done

	wait

else
	for VER in 2.6 2.7 3.2 3.3 ; do
		generate_wheels $VER
	done
fi
