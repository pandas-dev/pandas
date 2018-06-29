#!/bin/sh
command -v coverage >/dev/null && coverage erase
command -v python-coverage >/dev/null && python-coverage erase
pytest pandas --cov=pandas -r sxX --strict
