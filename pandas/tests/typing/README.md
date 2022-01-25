## Purpose of those tests

The tests contained in the `valid` directory are snippets that when
process through a type checker ensure that type annotations and type
stubs from this repository conform to common pandas API use-patterns.

## Running the tests

Tests can be run in following ways:

`pyright pandas/tests/typing`

`mypy pandas/tests/typing`

They'll also be automatically detected and executed by pytest. This
is to ensure that the test code itself is valid.

## Developing the tests

Some tests contain type checker ignore-instructions along with an
error that's supposed to be thrown.

    # error: No overload variant of "to_datetime" matches argument type "DataFrame"
    pd.to_datetime(df)  # type: ignore[call-overload]

All such constructs are placed because of the missing/invalid API
type information. When the API signature becomes valid again type
checker will ask you to remove `type: ignore`. Please remove the
above comment as well.

When adding new tests please use the above solution as well.

## Origins and attribution

The tests come from the [pandas-stubs](https://github.com/VirtusLab/pandas-stubs)
repository originally released under the MIT license.
