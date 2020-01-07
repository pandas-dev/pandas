#!/usr/bin/env bash

lint_files() {
    git diff upstream/master --name-only -u -- "*.py" | xargs -r $1
}

echo "FLAKE8 check: it passes if no errors are reported"
echo ""

lint_files flake8

echo ""
echo "ISORT check: it passes if no files need to be fixed"
echo "If fixes are made, please commit those with your changes"
echo ""

lint_files isort

echo ""
echo "BLACK check: it passes if no files need to be reformatted"
echo "If reformats are done, please commit those with your changes"
echo ""

lint_files black
