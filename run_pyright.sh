#!/bin/bash

# run pyright in mypy's pre-commit environment

# find mypy's environment
mypy_path=$(find ~/.cache/pre-commit/ -type f -name "mypy")

# remove "bin/mypy"
mypy_env=$(dirname $(dirname $mypy_path))

env_name=$(basename $mypy_env)
env_folder=$(dirname $mypy_env)

# the environment name needs to be specified in the config file
sed -i "s/tool.pyright\]/tool.pyright\]\nvenv = \"$env_name\"/" pyproject.toml

# the folder can be specified as an argument
RET=0
pyright --venv-path $env_folder
RET=$(($RET + $?))

# revert sed
git checkout pyproject.toml

exit $RET
