.. _io.latex:

=====
LaTeX
=====

.. versionadded:: 1.3.0

Currently there are no methods to read from LaTeX, only output methods.

Writing to LaTeX files
''''''''''''''''''''''

.. note::

   DataFrame *and* Styler objects currently have a ``to_latex`` method. We recommend
   using the :func:`Styler.to_latex` method over :func:`DataFrame.to_latex` due to the
   former's greater flexibility with conditional styling, and the latter's possible
   future deprecation.

Review the documentation for :func:`Styler.to_latex`, which gives examples of
conditional styling and explains the operation of its keyword arguments.

For simple application the following pattern is sufficient.

.. ipython:: python

   df = pd.DataFrame([[1, 2], [3, 4]], index=["a", "b"], columns=["c", "d"])
   print(df.style.to_latex())

To format values before output, chain the :func:`Styler.format` method.

.. ipython:: python

   print(df.style.format("€ {}").to_latex())
