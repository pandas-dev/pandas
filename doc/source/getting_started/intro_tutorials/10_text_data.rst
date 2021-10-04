.. _10min_tut_10_text:

{{ header }}

.. ipython:: python

    import pandas as pd

.. raw:: html

    <div class="card gs-data">
        <div class="card-header">
            <div class="gs-data-title">
                Data used for this tutorial:
            </div>
        </div>
        <ul class="list-group list-group-flush">
            <li class="list-group-item">
.. include:: includes/titanic.rst

.. ipython:: python

    titanic = pd.read_csv("data/titanic.csv")
    titanic.head()

.. raw:: html

            </li>
        </ul>
    </div>

How to manipulate textual data?
-------------------------------

.. raw:: html

    <ul class="task-bullet">
        <li>

Make all name characters lowercase.

.. ipython:: python

    titanic["Name"].str.lower()

To make each of the strings in the ``Name`` column lowercase, select the ``Name`` column
(see the :ref:`tutorial on selection of data <10min_tut_03_subset>`), add the ``str`` accessor and
apply the ``lower`` method. As such, each of the strings is converted element-wise.

.. raw:: html

        </li>
    </ul>

Similar to datetime objects in the :ref:`time series tutorial <10min_tut_09_timeseries>`
having a ``dt`` accessor, a number of
specialized string methods are available when using the ``str``
accessor. These methods have in general matching names with the
equivalent built-in string methods for single elements, but are applied
element-wise (remember :ref:`element-wise calculations <10min_tut_05_columns>`?)
on each of the values of the columns.

.. raw:: html

    <ul class="task-bullet">
        <li>

Create a new column ``Surname`` that contains the surname of the passengers by extracting the part before the comma.

.. ipython:: python

    titanic["Name"].str.split(",")

Using the :meth:`Series.str.split` method, each of the values is returned as a list of
2 elements. The first element is the part before the comma and the
second element is the part after the comma.

.. ipython:: python

    titanic["Surname"] = titanic["Name"].str.split(",").str.get(0)
    titanic["Surname"]

As we are only interested in the first part representing the surname
(element 0), we can again use the ``str`` accessor and apply :meth:`Series.str.get` to
extract the relevant part. Indeed, these string functions can be
concatenated to combine multiple functions at once!

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

More information on extracting parts of strings is available in the user guide section on :ref:`splitting and replacing strings <text.split>`.

.. raw:: html

   </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

Extract the passenger data about the countesses on board of the Titanic.

.. ipython:: python

    titanic["Name"].str.contains("Countess")

.. ipython:: python

    titanic[titanic["Name"].str.contains("Countess")]

(*Interested in her story? See* `Wikipedia <https://en.wikipedia.org/wiki/No%C3%ABl_Leslie,_Countess_of_Rothes>`__\ *!*)

The string method :meth:`Series.str.contains` checks for each of the values in the
column ``Name`` if the string contains the word ``Countess`` and returns
for each of the values ``True`` (``Countess`` is part of the name) or
``False`` (``Countess`` is not part of the name). This output can be used
to subselect the data using conditional (boolean) indexing introduced in
the :ref:`subsetting of data tutorial <10min_tut_03_subset>`. As there was
only one countess on the Titanic, we get one row as a result.

.. raw:: html

        </li>
    </ul>

.. note::
    More powerful extractions on strings are supported, as the
    :meth:`Series.str.contains` and :meth:`Series.str.extract` methods accept `regular
    expressions <https://docs.python.org/3/library/re.html>`__, but out of
    scope of this tutorial.

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

More information on extracting parts of strings is available in the user guide section on :ref:`string matching and extracting <text.extract>`.

.. raw:: html

   </div>

.. raw:: html

    <ul class="task-bullet">
        <li>

Which passenger of the Titanic has the longest name?

.. ipython:: python

    titanic["Name"].str.len()

To get the longest name we first have to get the lengths of each of the
names in the ``Name`` column. By using pandas string methods, the
:meth:`Series.str.len` function is applied to each of the names individually
(element-wise).

.. ipython:: python

    titanic["Name"].str.len().idxmax()

Next, we need to get the corresponding location, preferably the index
label, in the table for which the name length is the largest. The
:meth:`~Series.idxmax` method does exactly that. It is not a string method and is
applied to integers, so no ``str`` is used.

.. ipython:: python

    titanic.loc[titanic["Name"].str.len().idxmax(), "Name"]

Based on the index name of the row (``307``) and the column (``Name``),
we can do a selection using the ``loc`` operator, introduced in the
`tutorial on subsetting <3_subset_data.ipynb>`__.

.. raw:: html

        </li>
    </ul>

.. raw:: html

    <ul class="task-bullet">
        <li>

In the "Sex" column, replace values of "male" by "M" and values of "female" by "F".

.. ipython:: python

    titanic["Sex_short"] = titanic["Sex"].replace({"male": "M", "female": "F"})
    titanic["Sex_short"]

Whereas :meth:`~Series.replace` is not a string method, it provides a convenient way
to use mappings or vocabularies to translate certain values. It requires
a ``dictionary`` to define the mapping ``{from : to}``.

.. raw:: html

        </li>
    </ul>

.. warning::
    There is also a :meth:`~Series.str.replace` method available to replace a
    specific set of characters. However, when having a mapping of multiple
    values, this would become:

    ::

        titanic["Sex_short"] = titanic["Sex"].str.replace("female", "F")
        titanic["Sex_short"] = titanic["Sex_short"].str.replace("male", "M")

    This would become cumbersome and easily lead to mistakes. Just think (or
    try out yourself) what would happen if those two statements are applied
    in the opposite order…

.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>REMEMBER</h4>

-  String methods are available using the ``str`` accessor.
-  String methods work element-wise and can be used for conditional
   indexing.
-  The ``replace`` method is a convenient method to convert values
   according to a given dictionary.

.. raw:: html

   </div>

.. raw:: html

    <div class="d-flex flex-row gs-torefguide">
        <span class="badge badge-info">To user guide</span>

A full overview is provided in the user guide pages on :ref:`working with text data <text>`.

.. raw:: html

   </div>
