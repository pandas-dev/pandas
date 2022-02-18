.. _code_style:

{{ header }}

=======================
pandas code style guide
=======================

.. contents:: Table of contents:
   :local:

Patterns
========

We use a ``flake8`` plugin, `pandas-dev-flaker <https://github.com/pandas-dev/pandas-dev-flaker>`_, to
check our codebase for unwanted patterns. See its ``README`` for the up-to-date list of rules we enforce.

Miscellaneous
=============

Reading from a url
------------------

**Good:**

.. code-block:: python

    from pandas.io.common import urlopen

    with urlopen("http://www.google.com") as url:
        raw_text = url.read()
