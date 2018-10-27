.. _docstring:

======================
pandas docstring guide
======================

.. note::
  `Video tutorial: Pandas docstring guide
  <https://www.youtube.com/watch?v=EOA0lUeW4NI>`_ by Frank Akogun.

About docstrings and standards
------------------------------

A Python docstring is a string used to document a Python module, class,
function or method, so programmers can understand what it does without having
to read the details of the implementation.

Also, it is a common practice to generate online (html) documentation
automatically from docstrings. `Sphinx <http://www.sphinx-doc.org>`_ serves
this purpose.

Next example gives an idea on how a docstring looks like:

.. code-block:: python

    def add(num1, num2):
    """
    Add up two integer numbers.

    This function simply wraps the `+` operator, and does not
    do anything interesting, except for illustrating what is
    the docstring of a very simple function.

    Parameters
    ----------
    num1 : int
        First number to add
    num2 : int
        Second number to add

    Returns
    -------
    int
        The sum of `num1` and `num2`

    See Also
    --------
    subtract : Subtract one integer from another

    Examples
    --------
    >>> add(2, 2)
    4
    >>> add(25, 0)
    25
    >>> add(10, -10)
    0
    """
    return num1 + num2

Some standards exist about docstrings, so they are easier to read, and they can
be exported to other formats such as html or pdf.

The first conventions every Python docstring should follow are defined in
`PEP-257 <https://www.python.org/dev/peps/pep-0257/>`_.

As PEP-257 is quite open, and some other standards exist on top of it. In the
case of pandas, the numpy docstring convention is followed. The conventions is
explained in this document:

* `numpydoc docstring guide <http://numpydoc.readthedocs.io/en/latest/format.html>`_
  (which is based in the original `Guide to NumPy/SciPy documentation
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_)

numpydoc is a Sphinx extension to support the numpy docstring convention.

The standard uses reStructuredText (reST). reStructuredText is a markup
language that allows encoding styles in plain text files. Documentation
about reStructuredText can be found in:

* `Sphinx reStructuredText primer <http://www.sphinx-doc.org/en/stable/rest.html>`_
* `Quick reStructuredText reference <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
* `Full reStructuredText specification <http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html>`_

Pandas has some helpers for sharing docstrings between related classes, see
:ref:`docstring.sharing`.

The rest of this document will summarize all the above guides, and will
provide additional convention specific to the pandas project.

.. _docstring.tutorial:

Writing a docstring
-------------------

.. _docstring.general:

General rules
~~~~~~~~~~~~~

Docstrings must be defined with three double-quotes. No blank lines should be
left before or after the docstring. The text starts in the next line after the
opening quotes. The closing quotes have their own line
(meaning that they are not at the end of the last sentence).

In rare occasions reST styles like bold text or italics will be used in
docstrings, but is it common to have inline code, which is presented between
backticks. It is considered inline code:

* The name of a parameter
* Python code, a module, function, built-in, type, literal... (e.g. ``os``,
  ``list``, ``numpy.abs``, ``datetime.date``, ``True``)
* A pandas class (in the form ``:class:`pandas.Series```)
* A pandas method (in the form ``:meth:`pandas.Series.sum```)
* A pandas function (in the form ``:func:`pandas.to_datetime```)

.. note::
    To display only the last component of the linked class, method or
    function, prefix it with ``~``. For example, ``:class:`~pandas.Series```
    will link to ``pandas.Series`` but only display the last part, ``Series``
    as the link text. See `Sphinx cross-referencing syntax
    <http://www.sphinx-doc.org/en/stable/domains.html#cross-referencing-syntax>`_
    for details.

**Good:**

.. code-block:: python

    def add_values(arr):
        """
        Add the values in `arr`.

        This is equivalent to Python `sum` of :meth:`pandas.Series.sum`.

        Some sections are omitted here for simplicity.
        """
        return sum(arr)

**Bad:**

.. code-block:: python

    def func():

        """Some function.

        With several mistakes in the docstring.

        It has a blank like after the signature `def func():`.

        The text 'Some function' should go in the line after the
        opening quotes of the docstring, not in the same line.

        There is a blank line between the docstring and the first line
        of code `foo = 1`.

        The closing quotes should be in the next line, not in this one."""

        foo = 1
        bar = 2
        return foo + bar

.. _docstring.short_summary:

Section 1: Short summary
~~~~~~~~~~~~~~~~~~~~~~~~

The short summary is a single sentence that expresses what the function does in
a concise way.

The short summary must start with a capital letter, end with a dot, and fit in
a single line. It needs to express what the object does without providing
details. For functions and methods, the short summary must start with an
infinitive verb.

**Good:**

.. code-block:: python

    def astype(dtype):
        """
        Cast Series type.

        This section will provide further details.
        """
        pass

**Bad:**

.. code-block:: python

    def astype(dtype):    # noqa: F811
        """
        Casts Series type.

        Verb in third-person of the present simple, should be infinitive.
        """
        pass


    def astype(dtype):    # noqa: F811
        """
        Method to cast Series type.

        Does not start with verb.
        """
        pass


    def astype(dtype):    # noqa: F811
        """
        Cast Series type

        Missing dot at the end.
        """
        pass


    def astype(dtype):    # noqa: F811
        """
        Cast Series type from its current type to the new type defined in
        the parameter dtype.

        Summary is too verbose and doesn't fit in a single line.
        """
        pass

.. _docstring.extended_summary:

Section 2: Extended summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The extended summary provides details on what the function does. It should not
go into the details of the parameters, or discuss implementation notes, which
go in other sections.

A blank line is left between the short summary and the extended summary. And
every paragraph in the extended summary is finished by a dot.

The extended summary should provide details on why the function is useful and
their use cases, if it is not too generic.

.. code-block:: python

    def unstack():
        """
        Pivot a row index to columns.

        When using a MultiIndex, a level can be pivoted so each value in
        the index becomes a column. This is especially useful when a subindex
        is repeated for the main index, and data is easier to visualize as a
        pivot table.

        The index level will be automatically removed from the index when added
        as columns.
        """
        pass

.. _docstring.parameters:

Section 3: Parameters
~~~~~~~~~~~~~~~~~~~~~

The details of the parameters will be added in this section. This section has
the title "Parameters", followed by a line with a hyphen under each letter of
the word "Parameters". A blank line is left before the section title, but not
after, and not between the line with the word "Parameters" and the one with
the hyphens.

After the title, each parameter in the signature must be documented, including
`*args` and `**kwargs`, but not `self`.

The parameters are defined by their name, followed by a space, a colon, another
space, and the type (or types). Note that the space between the name and the
colon is important. Types are not defined for `*args` and `**kwargs`, but must
be defined for all other parameters. After the parameter definition, it is
required to have a line with the parameter description, which is indented, and
can have multiple lines. The description must start with a capital letter, and
finish with a dot.

For keyword arguments with a default value, the default will be listed after a
comma at the end of the type. The exact form of the type in this case will be
"int, default 0". In some cases it may be useful to explain what the default
argument means, which can be added after a comma "int, default -1, meaning all
cpus".

In cases where the default value is `None`, meaning that the value will not be
used. Instead of "str, default None", it is preferred to write "str, optional".
When `None` is a value being used, we will keep the form "str, default None".
For example, in `df.to_csv(compression=None)`, `None` is not a value being used,
but means that compression is optional, and no compression is being used if not
provided. In this case we will use `str, optional`. Only in cases like
`func(value=None)` and `None` is being used in the same way as `0` or `foo`
would be used, then we will specify "str, int or None, default None".

**Good:**

.. code-block:: python

    class Series:
        def plot(self, kind, color='blue', **kwargs):
            """
            Generate a plot.

            Render the data in the Series as a matplotlib plot of the
            specified kind.

            Parameters
            ----------
            kind : str
                Kind of matplotlib plot.
            color : str, default 'blue'
                Color name or rgb code.
            **kwargs
                These parameters will be passed to the matplotlib plotting
                function.
            """
            pass

**Bad:**

.. code-block:: python

    class Series:
        def plot(self, kind, **kwargs):
            """
            Generate a plot.

            Render the data in the Series as a matplotlib plot of the
            specified kind.

            Note the blank line between the parameters title and the first
            parameter. Also, note that after the name of the parameter `kind`
            and before the colon, a space is missing.

            Also, note that the parameter descriptions do not start with a
            capital letter, and do not finish with a dot.

            Finally, the `**kwargs` parameter is missing.

            Parameters
            ----------

            kind: str
                kind of matplotlib plot
            """
            pass

.. _docstring.parameter_types:

Parameter types
^^^^^^^^^^^^^^^

When specifying the parameter types, Python built-in data types can be used
directly (the Python type is preferred to the more verbose string, integer,
boolean, etc):

* int
* float
* str
* bool

For complex types, define the subtypes. For `dict` and `tuple`, as more than
one type is present, we use the brackets to help read the type (curly brackets
for `dict` and normal brackets for `tuple`):

* list of int
* dict of {str : int}
* tuple of (str, int, int)
* tuple of (str,)
* set of str

In case where there are just a set of values allowed, list them in curly
brackets and separated by commas (followed by a space). If the values are
ordinal and they have an order, list them in this order. Otherwise, list
the default value first, if there is one:

* {0, 10, 25}
* {'simple', 'advanced'}
* {'low', 'medium', 'high'}
* {'cat', 'dog', 'bird'}

If the type is defined in a Python module, the module must be specified:

* datetime.date
* datetime.datetime
* decimal.Decimal

If the type is in a package, the module must be also specified:

* numpy.ndarray
* scipy.sparse.coo_matrix

If the type is a pandas type, also specify pandas except for Series and
DataFrame:

* Series
* DataFrame
* pandas.Index
* pandas.Categorical
* pandas.SparseArray

If the exact type is not relevant, but must be compatible with a numpy
array, array-like can be specified. If Any type that can be iterated is
accepted, iterable can be used:

* array-like
* iterable

If more than one type is accepted, separate them by commas, except the
last two types, that need to be separated by the word 'or':

* int or float
* float, decimal.Decimal or None
* str or list of str

If ``None`` is one of the accepted values, it always needs to be the last in
the list.

For axis, the convention is to use something like:

* axis : {0 or 'index', 1 or 'columns', None}, default None

.. _docstring.returns:

Section 4: Returns or Yields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the method returns a value, it will be documented in this section. Also
if the method yields its output.

The title of the section will be defined in the same way as the "Parameters".
With the names "Returns" or "Yields" followed by a line with as many hyphens
as the letters in the preceding word.

The documentation of the return is also similar to the parameters. But in this
case, no name will be provided, unless the method returns or yields more than
one value (a tuple of values).

The types for "Returns" and "Yields" are the same as the ones for the
"Parameters". Also, the description must finish with a dot.

For example, with a single value:

.. code-block:: python

    def sample():
        """
        Generate and return a random number.

        The value is sampled from a continuous uniform distribution between
        0 and 1.

        Returns
        -------
        float
            Random number generated.
        """
        return random.random()

With more than one value:

.. code-block:: python

    def random_letters():
        """
        Generate and return a sequence of random letters.

        The length of the returned string is also random, and is also
        returned.

        Returns
        -------
        length : int
            Length of the returned string.
        letters : str
            String of random letters.
        """
        length = random.randint(1, 10)
        letters = ''.join(random.choice(string.ascii_lowercase)
                          for i in range(length))
        return length, letters

If the method yields its value:

.. code-block:: python

    def sample_values():
        """
        Generate an infinite sequence of random numbers.

        The values are sampled from a continuous uniform distribution between
        0 and 1.

        Yields
        ------
        float
            Random number generated.
        """
        while True:
            yield random.random()

.. _docstring.see_also:

Section 5: See Also
~~~~~~~~~~~~~~~~~~~

This section is used to let users know about pandas functionality
related to the one being documented. In rare cases, if no related methods
or functions can be found at all, this section can be skipped.

An obvious example would be the `head()` and `tail()` methods. As `tail()` does
the equivalent as `head()` but at the end of the `Series` or `DataFrame`
instead of at the beginning, it is good to let the users know about it.

To give an intuition on what can be considered related, here there are some
examples:

* ``loc`` and ``iloc``, as they do the same, but in one case providing indices
  and in the other positions
* ``max`` and ``min``, as they do the opposite
* ``iterrows``, ``itertuples`` and ``iteritems``, as it is easy that a user
  looking for the method to iterate over columns ends up in the method to
  iterate over rows, and vice-versa
* ``fillna`` and ``dropna``, as both methods are used to handle missing values
* ``read_csv`` and ``to_csv``, as they are complementary
* ``merge`` and ``join``, as one is a generalization of the other
* ``astype`` and ``pandas.to_datetime``, as users may be reading the
  documentation of ``astype`` to know how to cast as a date, and the way to do
  it is with ``pandas.to_datetime``
* ``where`` is related to ``numpy.where``, as its functionality is based on it

When deciding what is related, you should mainly use your common sense and
think about what can be useful for the users reading the documentation,
especially the less experienced ones.

When relating to other libraries (mainly ``numpy``), use the name of the module
first (not an alias like ``np``). If the function is in a module which is not
the main one, like ``scipy.sparse``, list the full module (e.g.
``scipy.sparse.coo_matrix``).

This section, as the previous, also has a header, "See Also" (note the capital
S and A). Also followed by the line with hyphens, and preceded by a blank line.

After the header, we will add a line for each related method or function,
followed by a space, a colon, another space, and a short description that
illustrated what this method or function does, why is it relevant in this
context, and what are the key differences between the documented function and
the one referencing. The description must also finish with a dot.

Note that in "Returns" and "Yields", the description is located in the
following line than the type. But in this section it is located in the same
line, with a colon in between. If the description does not fit in the same
line, it can continue in the next ones, but it has to be indented in them.

For example:

.. code-block:: python

    class Series:
        def head(self):
            """
            Return the first 5 elements of the Series.

            This function is mainly useful to preview the values of the
            Series without displaying the whole of it.

            Returns
            -------
            Series
                Subset of the original series with the 5 first values.

            See Also
            --------
            Series.tail : Return the last 5 elements of the Series.
            Series.iloc : Return a slice of the elements in the Series,
                which can also be used to return the first or last n.
            """
            return self.iloc[:5]

.. _docstring.notes:

Section 6: Notes
~~~~~~~~~~~~~~~~

This is an optional section used for notes about the implementation of the
algorithm. Or to document technical aspects of the function behavior.

Feel free to skip it, unless you are familiar with the implementation of the
algorithm, or you discover some counter-intuitive behavior while writing the
examples for the function.

This section follows the same format as the extended summary section.

.. _docstring.examples:

Section 7: Examples
~~~~~~~~~~~~~~~~~~~

This is one of the most important sections of a docstring, even if it is
placed in the last position. As often, people understand concepts better
with examples, than with accurate explanations.

Examples in docstrings, besides illustrating the usage of the function or
method, must be valid Python code, that in a deterministic way returns the
presented output, and that can be copied and run by users.

They are presented as a session in the Python terminal. `>>>` is used to
present code. `...` is used for code continuing from the previous line.
Output is presented immediately after the last line of code generating the
output (no blank lines in between). Comments describing the examples can
be added with blank lines before and after them.

The way to present examples is as follows:

1. Import required libraries (except ``numpy`` and ``pandas``)

2. Create the data required for the example

3. Show a very basic example that gives an idea of the most common use case

4. Add examples with explanations that illustrate how the parameters can be
   used for extended functionality

A simple example could be:

.. code-block:: python

    class Series:

        def head(self, n=5):
            """
            Return the first elements of the Series.

            This function is mainly useful to preview the values of the
            Series without displaying the whole of it.

            Parameters
            ----------
            n : int
                Number of values to return.

            Return
            ------
            pandas.Series
                Subset of the original series with the n first values.

            See Also
            --------
            tail : Return the last n elements of the Series.

            Examples
            --------
            >>> s = pd.Series(['Ant', 'Bear', 'Cow', 'Dog', 'Falcon',
            ...                'Lion', 'Monkey', 'Rabbit', 'Zebra'])
            >>> s.head()
            0   Ant
            1   Bear
            2   Cow
            3   Dog
            4   Falcon
            dtype: object

            With the `n` parameter, we can change the number of returned rows:

            >>> s.head(n=3)
            0   Ant
            1   Bear
            2   Cow
            dtype: object
            """
            return self.iloc[:n]

The examples should be as concise as possible. In cases where the complexity of
the function requires long examples, is recommended to use blocks with headers
in bold. Use double star ``**`` to make a text bold, like in ``**this example**``.

.. _docstring.example_conventions:

Conventions for the examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Code in examples is assumed to always start with these two lines which are not
shown:

.. code-block:: python

    import numpy as np          # noqa: F401
    import pandas as pd         # noqa: F401

Any other module used in the examples must be explicitly imported, one per line (as
recommended in :pep:`8#imports`)
and avoiding aliases. Avoid excessive imports, but if needed, imports from
the standard library go first, followed by third-party libraries (like
matplotlib).

When illustrating examples with a single ``Series`` use the name ``s``, and if
illustrating with a single ``DataFrame`` use the name ``df``. For indices,
``idx`` is the preferred name. If a set of homogeneous ``Series`` or
``DataFrame`` is used, name them ``s1``, ``s2``, ``s3``...  or ``df1``,
``df2``, ``df3``... If the data is not homogeneous, and more than one structure
is needed, name them with something meaningful, for example ``df_main`` and
``df_to_join``.

Data used in the example should be as compact as possible. The number of rows
is recommended to be around 4, but make it a number that makes sense for the
specific example. For example in the ``head`` method, it requires to be higher
than 5, to show the example with the default values. If doing the ``mean``, we
could use something like ``[1, 2, 3]``, so it is easy to see that the value
returned is the mean.

For more complex examples (grouping for example), avoid using data without
interpretation, like a matrix of random numbers with columns A, B, C, D...
And instead use a meaningful example, which makes it easier to understand the
concept. Unless required by the example, use names of animals, to keep examples
consistent. And numerical properties of them.

When calling the method, keywords arguments ``head(n=3)`` are preferred to
positional arguments ``head(3)``.

**Good:**

.. code-block:: python

    class Series:

        def mean(self):
            """
            Compute the mean of the input.

            Examples
            --------
            >>> s = pd.Series([1, 2, 3])
            >>> s.mean()
            2
            """
            pass


        def fillna(self, value):
            """
            Replace missing values by `value`.

            Examples
            --------
            >>> s = pd.Series([1, np.nan, 3])
            >>> s.fillna(0)
            [1, 0, 3]
            """
            pass

        def groupby_mean(self):
            """
            Group by index and return mean.

            Examples
            --------
            >>> s = pd.Series([380., 370., 24., 26],
            ...               name='max_speed',
            ...               index=['falcon', 'falcon', 'parrot', 'parrot'])
            >>> s.groupby_mean()
            index
            falcon    375.0
            parrot     25.0
            Name: max_speed, dtype: float64
            """
            pass

        def contains(self, pattern, case_sensitive=True, na=numpy.nan):
            """
            Return whether each value contains `pattern`.

            In this case, we are illustrating how to use sections, even
            if the example is simple enough and does not require them.

            Examples
            --------
            >>> s = pd.Series('Antelope', 'Lion', 'Zebra', numpy.nan)
            >>> s.contains(pattern='a')
            0    False
            1    False
            2     True
            3      NaN
            dtype: bool

            **Case sensitivity**

            With `case_sensitive` set to `False` we can match `a` with both
            `a` and `A`:

            >>> s.contains(pattern='a', case_sensitive=False)
            0     True
            1    False
            2     True
            3      NaN
            dtype: bool

            **Missing values**

            We can fill missing values in the output using the `na` parameter:

            >>> s.contains(pattern='a', na=False)
            0    False
            1    False
            2     True
            3    False
            dtype: bool
            """
            pass

**Bad:**

.. code-block:: python

    def method(foo=None, bar=None):
        """
        A sample DataFrame method.

        Do not import numpy and pandas.

        Try to use meaningful data, when it makes the example easier
        to understand.

        Try to avoid positional arguments like in `df.method(1)`. They
        can be all right if previously defined with a meaningful name,
        like in `present_value(interest_rate)`, but avoid them otherwise.

        When presenting the behavior with different parameters, do not place
        all the calls one next to the other. Instead, add a short sentence
        explaining what the example shows.

        Examples
        --------
        >>> import numpy as np
        >>> import pandas as pd
        >>> df = pd.DataFrame(numpy.random.randn(3, 3),
        ...                   columns=('a', 'b', 'c'))
        >>> df.method(1)
        21
        >>> df.method(bar=14)
        123
        """
        pass


.. _docstring.doctest_tips:

Tips for getting your examples pass the doctests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Getting the examples pass the doctests in the validation script can sometimes
be tricky. Here are some attention points:

* Import all needed libraries (except for pandas and numpy, those are already
  imported as ``import pandas as pd`` and ``import numpy as np``) and define
  all variables you use in the example.

* Try to avoid using random data. However random data might be OK in some
  cases, like if the function you are documenting deals with probability
  distributions, or if the amount of data needed to make the function result
  meaningful is too much, such that creating it manually is very cumbersome.
  In those cases, always use a fixed random seed to make the generated examples
  predictable. Example::

    >>> np.random.seed(42)
    >>> df = pd.DataFrame({'normal': np.random.normal(100, 5, 20)})

* If you have a code snippet that wraps multiple lines, you need to use '...'
  on the continued lines: ::

    >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=['a', 'b', 'c'],
    ...                   columns=['A', 'B'])

* If you want to show a case where an exception is raised, you can do::

    >>> pd.to_datetime(["712-01-01"])
    Traceback (most recent call last):
    OutOfBoundsDatetime: Out of bounds nanosecond timestamp: 712-01-01 00:00:00

  It is essential to include the "Traceback (most recent call last):", but for
  the actual error only the error name is sufficient.

* If there is a small part of the result that can vary (e.g. a hash in an object
  representation), you can use ``...`` to represent this part.

  If you want to show that ``s.plot()`` returns a matplotlib AxesSubplot object,
  this will fail the doctest ::

    >>> s.plot()
    <matplotlib.axes._subplots.AxesSubplot at 0x7efd0c0b0690>

  However, you can do (notice the comment that needs to be added) ::

    >>> s.plot()  # doctest: +ELLIPSIS
    <matplotlib.axes._subplots.AxesSubplot at ...>


.. _docstring.example_plots:

Plots in examples
^^^^^^^^^^^^^^^^^

There are some methods in pandas returning plots. To render the plots generated
by the examples in the documentation, the ``.. plot::`` directive exists.

To use it, place the next code after the "Examples" header as shown below. The
plot will be generated automatically when building the documentation.

.. code-block:: python

    class Series:
        def plot(self):
            """
            Generate a plot with the `Series` data.

            Examples
            --------

            .. plot::
                :context: close-figs

                >>> s = pd.Series([1, 2, 3])
                >>> s.plot()
            """
            pass

.. _docstring.sharing:

Sharing Docstrings
------------------

Pandas has a system for sharing docstrings, with slight variations, between
classes. This helps us keep docstrings consistent, while keeping things clear
for the user reading. It comes at the cost of some complexity when writing.

Each shared docstring will have a base template with variables, like
``%(klass)s``. The variables filled in later on using the ``Substitution``
decorator. Finally, docstrings can be appended to with the ``Appender``
decorator.

In this example, we'll create a parent docstring normally (this is like
``pandas.core.generic.NDFrame``. Then we'll have two children (like
``pandas.core.series.Series`` and ``pandas.core.frame.DataFrame``). We'll
substitute the children's class names in this docstring.

.. code-block:: python

   class Parent:
       def my_function(self):
           """Apply my function to %(klass)s."""
           ...


   class ChildA(Parent):
       @Substitution(klass="ChildA")
       @Appender(Parent.my_function.__doc__)
       def my_function(self):
           ...


   class ChildB(Parent):
       @Substitution(klass="ChildB")
       @Appender(Parent.my_function.__doc__)
       def my_function(self):
           ...

The resulting docstrings are

.. code-block:: python

   >>> print(Parent.my_function.__doc__)
   Apply my function to %(klass)s.
   >>> print(ChildA.my_function.__doc__)
   Apply my function to ChildA.
   >>> print(ChildB.my_function.__doc__)
   Apply my function to ChildB.

Notice two things:

1. We "append" the parent docstring to the children docstrings, which are
   initially empty.
2. Python decorators are applied inside out. So the order is Append then
   Substitution, even though Substitution comes first in the file.

Our files will often contain a module-level ``_shared_doc_kwargs`` with some
common substitution values (things like ``klass``, ``axes``, etc).

You can substitute and append in one shot with something like

.. code-block:: python

   @Appender(template % _shared_doc_kwargs)
   def my_function(self):
       ...

where ``template`` may come from a module-level ``_shared_docs`` dictionary
mapping function names to docstrings. Wherever possible, we prefer using
``Appender`` and ``Substitution``, since the docstring-writing processes is
slightly closer to normal.

See ``pandas.core.generic.NDFrame.fillna`` for an example template, and
``pandas.core.series.Series.fillna`` and ``pandas.core.generic.frame.fillna``
for the filled versions.
