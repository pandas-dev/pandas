.. _style:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows = 15

*******
Styling
*******

.. versionadded:: 0.17.1


You can apply **conditional formatting**, the visual styling of a DataFrame
depending on the data within, by using the ``.style`` property.
This is a property on every object that returns a ``Styler`` object, which has
useful methods for formatting and displaying DataFrames.

The styling is accomplished through CSS_.
You write functions that take DataFrames or Series, and return like-indexed
objects with CSS ``"attribute: value"`` pairs.
You can build up your styles incrementally using method chains, before rending.

.. _CSS: https://developer.mozilla.org/en-US/docs/Web/CSS

Let's create some dummy data to work with

.. ipython:: python

    np.random.seed(24)
    df = pd.DataFrame({'A': np.linspace(1, 10, 10)})
    df = pd.concat([df, pd.DataFrame(np.random.randn(10, 4), columns=list('BCDE'))],
                   axis=1)
    df.iloc[0, 2] = np.nan

.. ipython:: python

    df.style.bar().render()

