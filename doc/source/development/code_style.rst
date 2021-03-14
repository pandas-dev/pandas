.. _code_style:

{{ header }}

=======================
pandas code style guide
=======================

.. contents:: Table of contents:
   :local:

pandas follows the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_
standard and uses `Black <https://black.readthedocs.io/en/stable/>`_
and `Flake8 <https://flake8.pycqa.org/en/latest/>`_ to ensure a
consistent code format throughout the project. We encourage you to use
:ref:`pre-commit <contributing.pre-commit>` to automatically run ``black``,
``flake8``, ``isort``, and related code checks when you make a git commit.

Patterns
========

Using foo.__class__
-------------------


pandas uses 'type(foo)' instead 'foo.__class__' as it is making the code more
readable.
For example:

**Good:**

.. code-block:: python

    foo = "bar"
    type(foo)

**Bad:**

.. code-block:: python

    foo = "bar"
    foo.__class__


String formatting
=================

Concatenated strings
--------------------

Using f-strings
~~~~~~~~~~~~~~~

pandas uses f-strings formatting instead of '%' and '.format()' string formatters.

The convention of using f-strings on a string that is concatenated over several lines,
is to prefix only the lines containing values which need to be interpreted.

For example:

**Good:**

.. code-block:: python

    foo = "old_function"
    bar = "new_function"

    my_warning_message = (
        f"Warning, {foo} is deprecated, "
        "please use the new and way better "
        f"{bar}"
    )

**Bad:**

.. code-block:: python

    foo = "old_function"
    bar = "new_function"

    my_warning_message = (
        f"Warning, {foo} is deprecated, "
        f"please use the new and way better "
        f"{bar}"
    )

White spaces
~~~~~~~~~~~~

Only put white space at the end of the previous line, so
there is no whitespace at the beginning of the concatenated string.

For example:

**Good:**

.. code-block:: python

    example_string = (
        "Some long concatenated string, "
        "with good placement of the "
        "whitespaces"
    )

**Bad:**

.. code-block:: python

    example_string = (
        "Some long concatenated string,"
        " with bad placement of the"
        " whitespaces"
    )

Representation function (aka 'repr()')
--------------------------------------

pandas uses 'repr()' instead of '%r' and '!r'.

The use of 'repr()' will only happen when the value is not an obvious string.

For example:

**Good:**

.. code-block:: python

    value = str
    f"Unknown received value, got: {repr(value)}"

**Good:**

.. code-block:: python

    value = str
    f"Unknown received type, got: '{type(value).__name__}'"


Imports (aim for absolute)
==========================

In Python 3, absolute imports are recommended. Using absolute imports, doing something
like ``import string`` will import the string module rather than ``string.py``
in the same directory. As much as possible, you should try to write out
absolute imports that show the whole import chain from top-level pandas.

Explicit relative imports are also supported in Python 3 but it is not
recommended to use them. Implicit relative imports should never be used
and are removed in Python 3.

For example:

::

    # preferred
    import pandas.core.common as com

    # not preferred
    from .common import test_base

    # wrong
    from common import test_base

Testing
=======

Failing tests
--------------

See https://docs.pytest.org/en/latest/skipping.html for background.

Do not use ``pytest.xfail``
---------------------------

Do not use this method. It has the same behavior as ``pytest.skip``, namely
it immediately stops the test and does not check if the test will fail. If
this is the behavior you desire, use ``pytest.skip`` instead.

Using ``pytest.mark.xfail``
---------------------------

Use this method if a test is known to fail but the manner in which it fails
is not meant to be captured. It is common to use this method for a test that
exhibits buggy behavior or a non-implemented feature. If
the failing test has flaky behavior, use the argument ``strict=False``. This
will make it so pytest does not fail if the test happens to pass.

Prefer the decorator ``@pytest.mark.xfail`` and the argument ``pytest.param``
over usage within a test so that the test is appropriately marked during the
collection phase of pytest. For xfailing a test that involves multiple
parameters, a fixture, or a combination of these, it is only possible to
xfail during the testing phase. To do so, use the ``request`` fixture:

.. code-block:: python

    import pytest

    def test_xfail(request):
        request.node.add_marker(pytest.mark.xfail(reason="Indicate why here"))

xfail is not to be used for tests involving failure due to invalid user arguments.
For these tests, we need to verify the correct exception type and error message
is being raised, using ``pytest.raises`` instead.

Miscellaneous
=============

Reading from a url
------------------

**Good:**

.. code-block:: python

    from pandas.io.common import urlopen

    with urlopen("http://www.google.com") as url:
        raw_text = url.read()
