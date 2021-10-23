The equivalent pandas methods are :meth:`Series.str.upper`, :meth:`Series.str.lower`, and
:meth:`Series.str.title`.

.. ipython:: python

   firstlast = pd.DataFrame({"string": ["John Smith", "Jane Cook"]})
   firstlast["upper"] = firstlast["string"].str.upper()
   firstlast["lower"] = firstlast["string"].str.lower()
   firstlast["title"] = firstlast["string"].str.title()
   firstlast
