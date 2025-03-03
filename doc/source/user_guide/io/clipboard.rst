.. _io.clipboard:

=========
Clipboard
=========

A handy way to grab data is to use the :meth:`~DataFrame.read_clipboard` method,
which takes the contents of the clipboard buffer and passes them to the
``read_csv`` method. For instance, you can copy the following text to the
clipboard (CTRL-C on many operating systems):

.. code-block:: console

     A B C
   x 1 4 p
   y 2 5 q
   z 3 6 r

And then import the data directly to a ``DataFrame`` by calling:

.. code-block:: python

    >>> clipdf = pd.read_clipboard()
    >>> clipdf
      A B C
    x 1 4 p
    y 2 5 q
    z 3 6 r

The ``to_clipboard`` method can be used to write the contents of a ``DataFrame`` to
the clipboard. Following which you can paste the clipboard contents into other
applications (CTRL-V on many operating systems). Here we illustrate writing a
``DataFrame`` into clipboard and reading it back.

.. code-block:: python

    >>> df = pd.DataFrame(
    ...     {"A": [1, 2, 3], "B": [4, 5, 6], "C": ["p", "q", "r"]}, index=["x", "y", "z"]
    ... )

    >>> df
      A B C
    x 1 4 p
    y 2 5 q
    z 3 6 r
    >>> df.to_clipboard()
    >>> pd.read_clipboard()
      A B C
    x 1 4 p
    y 2 5 q
    z 3 6 r

We can see that we got the same content back, which we had earlier written to the clipboard.

.. note::

   You may need to install xclip or xsel (with PyQt5, PyQt4 or qtpy) on Linux to use these methods.
