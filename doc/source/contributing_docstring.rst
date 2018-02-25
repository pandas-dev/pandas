======================
pandas docstring guide
======================

About docstrings and standards
------------------------------

A Python docstring is a string used to document a Python function or method,
so programmers can understand what it does without having to read the details
of the implementation.

Also, it is a commonn practice to generate online (html) documentation
automatically from docstrings. `Sphinx <http://www.sphinx-doc.org>`_ serves
this purpose.

Next example gives an idea on how a docstring looks like:

.. code-block:: python

    def add(num1, num2):
    """Add up two integer numbers.

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
case of pandas, the numpy docstring convention is followed. There are two main
documents that explain this convention:

- `Guide to NumPy/SciPy documentation <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
- `numpydoc docstring guide <http://numpydoc.readthedocs.io/en/latest/format.html>`_

numpydoc is a Sphinx extension to support the numpy docstring convention.

The standard uses reStructuredText (reST). reStructuredText is a markup
language that allows encoding styles in plain text files. Documentation
about reStructuredText can be found in:

- `Sphinx reStructuredText primer <http://www.sphinx-doc.org/en/stable/rest.html>`_
- `Quick reStructuredText reference <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
- `Full reStructuredText specification <http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html>`_

The rest of this document will summarize all the above guides, and will
provide additional convention specific to the pandas project.

Writing a docstring
-------------------

General rules
~~~~~~~~~~~~~

Docstrings must be defined with three double-quotes. No blank lines should be
left before or after the docstring. The text starts immediately after the
opening quotes (not in the next line). The closing quotes have their own line
(meaning that they are not at the end of the last sentence).

**Good:**

.. code-block:: python

    def func():
        """Some function.

        With a good docstring.
        """
        foo = 1
        bar = 2
        return foo + bar

**Bad:**

.. code-block:: python

    def func():

        """
        Some function.

        With several mistakes in the docstring.
        
        It has a blank like after the signature `def func():`.
        
        The text 'Some function' should go in the same line as the
        opening quotes of the docstring, not in the next line.
        
        There is a blank line between the docstring and the first line
        of code `foo = 1`.
        
        The closing quotes should be in the next line, not in this one."""

        foo = 1
        bar = 2
        return foo + bar

Section 1: Short summary
~~~~~~~~~~~~~~~~~~~~~~~~

The short summary is a single sentence that express what the function does in a
concise way.

The short summary must start with a verb infinitive, end with a dot, and fit in
a single line. It needs to express what the function does without providing
details.

**Good:**

.. code-block:: python

    def astype(dtype):
        """Cast Series type.

        This section will provide further details.
        """
        pass

**Bad:**

.. code-block:: python

    def astype(dtype):
        """Casts Series type.

        Verb in third-person of the present simple, should be infinitive.
        """
        pass

    def astype(dtype):
        """Method to cast Series type.

        Does not start with verb.
        """
        pass

    def astype(dtype):
        """Cast Series type

        Missing dot at the end.
        """
        pass

    def astype(dtype):
        """Cast Series type from its current type to the new type defined in
        the parameter dtype.

        Summary is too verbose and doesn't fit in a single line.
        """
        pass

Section 2: Extended summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The extended summary provides details on what the function does. It should not
go into the details of the parameters, or discuss implementation notes, which
go in other sections.

A blank line is left between the short summary and the extended summary. And
every paragraph in the extended summary is finished by a dot.

.. code-block:: python

    def unstack():
        """Pivot a row index to columns.

        When using a multi-index, a level can be pivoted so each value in
        the index becomes a column. This is especially useful when a subindex
        is repeated for the main index, and data is easier to visualize as a
        pivot table.

        The index level will be automatically removed from the index when added
        as columns.
        """
        pass

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

Keyword arguments with a default value, the default will be listed in brackets
at the end of the description (before the dot). The exact form of the
description in this case would be "Description of the arg (default is X).". In
some cases it may be useful to explain what the default argument means, which
can be added after a comma "Description of the arg (default is -1, which means
all cpus).".

**Good:**

.. code-block:: python

    class Series:
        def plot(self, kind, color='blue', **kwargs):
            """Generate a plot.

            Render the data in the Series as a matplotlib plot of the
            specified kind.

            Parameters
            ----------
            kind : str
                Kind of matplotlib plot.
            color : str
                Color name or rgb code (default is 'blue').
            **kwargs
                These parameters will be passed to the matplotlib plotting
                function.
            """
            pass

**Bad:**

.. code-block:: python

    class Series:
        def plot(self, kind, **kwargs):
            """Generate a plot.

            Render the data in the Series as a matplotlib plot of the
            specified kind.

            Note the blank line between the parameters title and the first
            parameter. Also, not that after the name of the parameter `kind`
            and before the colo, a space is missing.

            Also, note that the parameter descriptions do not start with a
            capital letter, and do not finish with a dot.

            Finally, the `**kwargs` parameter is missing.

            Parameters
            ----------

            kind: str
                kind of matplotlib plot
            """
            pass

Parameter types
^^^^^^^^^^^^^^^

When specifying the parameter types, Python built-in data types can be used
directly:

- int
- float
- str

For complex types, define the subtypes:

- list of [int]
- dict of {str : int}
- tuple of (str, int, int)
- set of {str}

In case there are just a set of values allowed, list them in curly brackets
and separated by commas (followed by a space). If one of them is the default
value of a keyword argument, it should be listed first.:

- {0, 10, 25}
- {'simple', 'advanced'}

If the type is defined in a Python module, the module must be specified:

- datetime.date
- datetime.datetime
- decimal.Decimal

If the type is in a package, the module must be also specified:

- numpy.ndarray
- scipy.sparse.coo_matrix

If the type is a pandas type, also specify pandas except for Series and
DataFrame:

- Series
- DataFrame
- pandas.Index
- pandas.Categorical
- pandas.SparseArray

If the exact type is not relevant, but must be compatible with a numpy
array, array-like can be specified. If Any type that can be iterated is
accepted, iterable can be used:

- array-like
- iterable

If more than one type is accepted, separate them by commas, except the
last two types, that need to be separated by the word 'or':

- int or float
- float, decimal.Decimal or None
- str or list of str

If None is one of the accepted values, it always needs to be the last in
the list.

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
        """Generate and return a random number.

        The value is sampled from a continuos uniform distribution between
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
        """Generate and return a sequence of random letters.

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
        """Generate an infinite sequence of random numbers.

        The values are sampled from a continuos uniform distribution between
        0 and 1.

        Yields
        ------
        float
            Random number generated.
        """
        while True:
            yield random.random()


Section 5: See also
~~~~~~~~~~~~~~~~~~~

This is an optional section, used to let users know about pandas functionality
related to the one being documented. While optional, this section should exist
in most cases, unless no related methods or functions can be found at all.

An obvious example would be the `head()` and `tail()` methods. As `tail()` does
the equivalent as `head()` but at the end of the `Series` or `DataFrame`
instead of at the beginning, it is good to let the users know about it.

To give an intuition on what can be considered related, here there are some
examples:

* `loc` and `iloc`, as they do the same, but in one case providing indices and
  in the other positions
* `max` and `min`, as they do the opposite
* `iterrows`, `itertuples` and `iteritems`, as it is easy that a user looking
  for the method to iterate over columns ends up in the method to iterate
  over rows, and vice-versa
* `fillna` and `dropna`, as both methods are used to handle missing values
* `read_csv` and `to_csv`, as they are complementary
* `merge` and `join`, as one is a generalization of the other
* `astype` and `pandas.to_datetime`, as users may be reading the documentation
  of `astype` to know how to cast as a date, and the way to do it is with
  `pandas.to_datetime`
* `where` is related to `numpy.where`, as its functionality is based on it

When deciding what is related, you should mainly use your common sense and
think about what can be useful for the users reading the documentation,
especially the less experienced ones.

When relating to other methods (mainly `numpy`), use the name of the module
first (not an alias like `np`). If the function is in a module which is not
the main one, like `scipy.sparse`, list the full module (e.g.
`scipy.sparse.coo_matrix`).

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
line, it can continue in the next ones, but it has to be indenteted in them.

For example:

.. code-block:: python

    class Series:
        def head(self):
            """Return the first 5 elements of the Series.

            This function is mainly useful to preview the values of the
            Series without displaying the whole of it.

            Return
            ------
            pandas.Series
                Subset of the original series with the 5 first values.

            See Also
            --------
            tail : Return the last 5 elements of the Series.
            """
            return self.iloc[:5]

Section 6: Notes
~~~~~~~~~~~~~~~~

This is an optional section used for notes about the implementation of the
algorithm. Or to document technical aspects of the function behavior.

Feel free to skip it, unless you are familiar with the implementation of the
algorithm, or you discover some counter-intuitive behavior while writing the
examples for the function.

This section follows the same format as the extended summary section.

Section 7: Examples
~~~~~~~~~~~~~~~~~~~

This is one of the most important sections of a docstring, even if it is
placed in the last position. As often, people understand concepts better
with examples, than with accurate explanations.

Examples in docstrings are also unit tests, and besides illustrating the
usage of the function or method, they need to be valid Python code, that in a
deterministic way returns the presented output.

They are presented as a session in the Python terminal. `>>>` is used to
present code. `...` is used for code continuing from the previous line.
Output is presented immediately after the last line of code generating the
output (no blank lines in between). Comments describing the examples can
be added with blank lines before and after them.

The way to present examples is as follows:

1. Import required libraries (except `numpy` and `pandas`)

2. Create the data required for the example

3. Show a very basic example that gives an idea of the most common use case

4. Add examples with explanations that illustrate how the parameters can be
   used for extended functionality

A simple example could be:

.. code-block:: python

    class Series:
        def head(self, n=5):
            """Return the first elements of the Series.

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

Conventions for the examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Code in examples is assumed to always start with these two lines which are not
shown:

.. code-block:: python

    import numpy as np
    import pandas as pd


Any other module used in the examples must be explicitly imported, one per line (as
recommended in `PEP-8 <https://www.python.org/dev/peps/pep-0008/#imports>`_)
and avoiding aliases. Avoid excessive imports, but if needed, imports from
the standard library go first, followed by third-party libraries (like
matplotlib).

When illustrating examples with a single `Series` use the name `s`, and if
illustrating with a single `DataFrame` use the name `df`. If a set of
homogeneous `Series` or `DataFrame` is used, name them `s1`, `s2`, `s3`...
or `df1`, `df2`, `df3`... If the data is not homogeneous, and more than
one structure is needed, name them with something meaningful, for example
`df_main` and `df_to_join`.

Data used in the example should be as compact as possible. The number of rows
is recommended to be 4, unless the example requires a larger number. As for
example in the head method, where it requires to be higher than 5, to show
the example with the default values.

Avoid using data without interpretation, like a matrix of random numbers
with columns A, B, C, D... And instead use a meaningful example, which makes
it easier to understand the concept. Unless required by the example, use
names of animals, to keep examples consistent. And numerical properties of
them.

When calling the method, keywords arguments `head(n=3)` are preferred to
positional arguments `head(3)`.

**Good:**

.. code-block:: python

    def method():
        """A sample DataFrame method.

        Examples
        --------
        >>> df = pd.DataFrame([389., 24., 80.5, numpy.nan]
        ...                   columns=('max_speed'),
        ...                   index=['falcon', 'parrot', 'lion', 'monkey'])
        """
        pass

**Bad:**

.. code-block:: python

    def method():
        """A sample DataFrame method.

        Examples
        --------
        >>> import numpy as np
        >>> import pandas as pd
        >>> df = pd.DataFrame(numpy.random.randn(3, 3),
        ...                   columns=('a', 'b', 'c'))
        """
        pass

Once you finished the docstring
-------------------------------

When you finished the changes to the docstring, go to the
:ref:`instructions to submit your changes <pandas_pr>` to continue.
