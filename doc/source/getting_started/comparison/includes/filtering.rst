:orphan:

DataFrames can be filtered in multiple ways; the most intuitive of which is using
:ref:`boolean indexing <indexing.boolean>`

.. ipython:: python
   :suppress:

   # ensure tips is defined when scanning with flake8-rst
   if 'tips' not in vars():
       tips = {}

.. ipython:: python

   tips[tips["total_bill"] > 10]

The above statement is simply passing a ``Series`` of ``True``/``False`` objects to the DataFrame,
returning all rows with ``True``.

.. ipython:: python

    is_dinner = tips["time"] == "Dinner"
    is_dinner
    is_dinner.value_counts()
    tips[is_dinner]
