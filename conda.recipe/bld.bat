@echo off

conda remove jinja2 --quiet
conda install jinja2 --quiet
%PYTHON% setup.py install
