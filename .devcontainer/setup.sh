#!/bin/bash

set -e

conda init --all
conda env create -f environment.yml

git submodule update --init
