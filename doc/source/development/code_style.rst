.. _code_style:

{{ header }}

=======================
pandas code style guide
=======================

.. contents:: Table of contents:
   :local:

Patterns
========

foo.__class__
-------------

*pandas* uses 'type(foo)' instead 'foo.__class__' as it is making the code more
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

f-strings
~~~~~~~~~

*pandas* uses f-strings formatting instead of '%' and '.format()' string formatters.

The convention of using f-strings on a string that is concatenated over serveral lines,
is to prefix only the lines containing the value needs to be interpeted.

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

Putting the white space only at the end of the previous line, so
there is no whitespace at the beggining of the concatenated string.

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

*pandas* uses 'repr()' instead of '%r' and '!r'.

The use of 'repr()' will only happend when the value is not an obvious string.

For example:

**Good:**

.. code-block:: python

    value = str
    f"Unknown recived value, got: {repr(value)}"

**Good:**

.. code-block:: python

    value = str
    f"Unknown recived type, got: '{type(value).__name__}'"
