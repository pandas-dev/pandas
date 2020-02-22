.. _testing:

{{ header }}

=======
Testing
=======

.. contents:: Table of contents:
   :local:

First off - thank you for writing test cases - they're really important in
developing pandas!

Typical Imports
===============

.. code-block:: python

    import nose
    import unittest
    import pandas.util.testing as tm
    from pandas.util.testing import makeCustomDataframe as mkdf

Making your tests behave well
=============================

``pandas`` committers run test cases after every change (as does Travis), so it's
important that you make your tests well-behaved. Balancing that, it's important
that your test cases cover the functionality of your addition, so that when
others make changes, they can be confident that they aren't introducing errors
in your code.  This includes:

1. marking network-using test cases with ``@network`` (see below).
2. marking slow tests with ``@slow``
3. using smaller test cases where it makes sense (for example, if you're
   testing a ``numexpr`` evaluation, you can generally just set ``expr._MIN_ELEMENTS = 0``
   and go ahead, rather than needing to test on a frame of at least 10K
   elements).
4. making sure to skip tests (or even test files) if a required import is not
   available.

In addition, stylistically, the preference is to group multiple related tests
under one test function and *not* to use the generator functionality of nose
in order to keep the actual # of tests small.

E.g.:

.. code-block:: python

    @slow
    def test_million_element_arithmetic():
        df = mkdf(100000, 100000)
        tm.assert_frame_equal(df.mod(df) * df * 0, df * 0)

Additional imports
------------------

When creating a subclass of ``unittest.TestCase`` there are useful instance
methods such as ``self.assertEqual(a, b)`` that allow you to test the equality
of two objects. These are not available *as functions* in the Python standard
library. However, these methods are available as functions in the ``nose.tools``
module. To use ``self.assertEqual(a, b)`` in a function you would put
``from nose.tools import assert_equal`` somewhere in the file and then call it
wherever you need it.

**Important**: make sure to document failure conditions (and use the
``assertRaisesRegexp`` where necessary to make it clearer *which* exception
you want to get). Testing for ``Exception`` is strongly discouraged.

Testing using a File
====================

This context manager allows safe read/write access to a temporary file, with
a generated filename (or your filename if provided). The file will be
automatically deleted when the context block is exited.

.. code-block:: python

    with tm.ensure_clean('my_file_path') as path:
        # do something with the path

Testing for Exceptions
======================

Generally, it's not acceptable to just check that something raises ``Exception``,
because that tends to mask a lot of errors. For example, if a function's signature
changes between releases, you could be catching the wrong kind of error altogether.
Going forward, the goal is to have no test cases that pass if ``Exception`` or a
subclass is raised (we're not quite there yet).

Another element that is helpful is to use ``assertRaisesRegexp`` from ``pandas.util.testing``.
It lets you be very explicit with what you expect (and prevents hiding errors like
changing signatures, etc.)

.. code-block:: python

    with tm.assertRaises(ValueError):
        raise ValueError("an error")
    with tm.assertRaisesRegexp(TypeError, 'invalid literal'):
        int('abc')

Handling tests requiring network connectivity
=============================================

**Please run your tests without an internet connection before submitting a PR!** (it's really important that your tests *not* fail when you have no internet connection (i.e., they should skip with out a network connection). In general, network tests are finicky. All tests that involve networking *must* be marked as "network", either by using the ``network`` decorator or the ``with_connectivity_check`` decorator from ``pandas.util.testing``.Unless you *absolutely* need to test that a function/method correctly handles connectivity errors, you should use the ``network`` decorator, which will catch all ``IOError`` s (which includes ``URLError``). If you believe that your test case will only fail if you simply aren't connected to the internet, you can use the ``with_connectivity_test`` to check:

.. code-block:: python

    >>> @with_connectivity_check
    ... def test_my_function():
    ...     urllib2.urlopen("funny://rabbithead")
    >>> test_my_function()
    Traceback (most recent call last)
        ...
    URLError...#some message

If you want to have the decorator always raise errors, just pass ``raise_on_error=True``
to the ``network`` decorator:

.. code-block:: python

    >>> @network(raise_on_error=True)
    ... def test2():
    ...     raise URLError("WRONG!")
    Traceback (most recent call last)
        ...
    URLError: WRONG!

The ``with_connectivity_check`` decorator defaults to checking ``http://www.google.com``
to determine whether it is connected. But if you had a test that depends on yahoo,
it might make sense to check yahoo instead:

.. code-block:: python

    @with_connectivity_check("http://www.yahoo.com")
    def some_test_with_yahoo():
        # do something etc.

It's a good idea to break up network tests into at least two parts:

1. Tests that check that the code works and gracefully handles errors.
2. Tests that really only matter if you have network connectivity (like making
   sure that the current Google Analytics feed is being processed properly).

For (1), you might want to use ``@network(raise_on_error=True)``, because those
tests should *not* fail without connectivity.

For (2), you should definitely suppress network errors, and, particularly if you
have a slow test, you may even want to check for connectivity *first* (so the
test never even runs if there isn't a network connection). You can do that easily
by passing ``check_before_test=True`` to ``with_connectivity_check``:

.. code-block:: python

    @with_connectivity_check("http://www.somespecificsite.com", check_before_test=True)
    def some_test():
        for i in range(1000):
            test_some_really_long_function(i)

Testing for Warnings
====================

To test for warnings, you can use the ``assert_produces_warning`` contextmanager,
which checks that your code produces a warning.

Probably the most common case is just a test case for a DeprecationWarning:

.. code-block:: python

    >>> with assert_produces_warning(DeprecationWarning):
    ...     some_function_that_raises_deprecation_warning()

With no arguments, it checks that any warning is raised.

.. code-block:: python

    >>> import warnings
    >>> with assert_produces_warning():
    ...     warnings.warn(UserWarning())
    ...

When passed False, it checks that *no* warnings are raised.

.. code-block:: python

    >>> with assert_produces_warning(False):
    ...     warnings.warn(RuntimeWarning())
    ...
    Traceback (most recent call last):
        ...
    AssertionError: Caused unexpected warning(s): ['RuntimeWarning'].

Finally, if you pass it a warning class, it will check that the *specific* 
class of warning was raised and no other.

.. code-block:: python

    >>> with assert_produces_warning(UserWarning):
    ...     warnings.warn(RuntimeWarning())
    Traceback (most recent call last):
        ...
    AssertionError: Did not see expected warning of class 'UserWarning'.

Reading from either a URL or zip file
=====================================

Reading from a url
------------------

.. code-block:: python

    from pandas.io.common import urlopen
    with urlopen('http://www.google.com') as url:
        raw_text = url.read()


Reading a file named ``file.txt`` that's inside of a zip file named ``file.zip``
--------------------------------------------------------------------------------

.. code-block:: python

    from pandas.io.common import ZipFile
    with ZipFile('file.zip') as zf:
        raw_text = zf.read('file.txt')

Hook up travis-ci
=================

We use travis for testings the entire library across various python versions.
If you [hook up your fork](http://about.travis-ci.org/docs/user/getting-started/)
to run travis, then it is displayed prominently whether your pull request passes
or fails the testing suite. This is incredibly helpful.

If it shows that it passes, great! We can consider merging. If there's a failure,
this let's you and us know there is something wrong, and needs some attention
before it can be considered for merging.

Sometimes Travis will say a change failed for reasons unrelated to your pull
request. For example there could be a build error or network error. To get Travis
to retest your pull request, do the following:

.. code-block:: shell

    git commit --amend -C HEAD
    git push origin <yourbranch> -f`
