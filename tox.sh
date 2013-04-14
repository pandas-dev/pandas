#!/usr/bin/env bash


if [ x"$1" == x"fast" ]; then
    scripts/use_build_cache.py
fi;

tox
