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

We use a ``flake8`` plugin, `pandas-dev-flaker <https://github.com/pandas-dev/pandas-dev-flaker>`_, to
check our codebase for unwanted patterns. See its ``README`` for the up-to-date list of rules we enforce.

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
        mark = pytest.mark.xfail(raises=TypeError, reason="Indicate why here")
        request.node.add_marker(mark)

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
