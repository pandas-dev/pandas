"""
Check that doc/source/reference/general_utility_functions.rst documents all exceptions and
warnings in  pandas/errors/__init__.py.

This is meant to be run as a pre-commit hook - to run it manually, you can do:

    pre-commit run pandas-errors-documented --all-files
"""
import inspect
import sys

from pandas import errors


API_PATH = "doc/source/reference/general_utility_functions.rst"


def main():
    with open(API_PATH, "r") as f:
        api_docs = f.read()
        missing = []
        for obj in errors.__dict__.values():
            if inspect.isclass(obj) and issubclass(obj, (Exception, Warning)) and obj.__name__ not in api_docs:
                missing.append(obj.__name__)
        if missing:
            sys.stdout.write(f"The follow exceptions and/or warnings are not documented in {API_PATH}: {missing}")
            sys.exit(1)
        sys.exit(0)


if __name__ == "__main__":
    main()
