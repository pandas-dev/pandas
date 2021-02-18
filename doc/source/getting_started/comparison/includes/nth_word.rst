The simplest way to extract words in pandas is to split the strings by spaces, then reference the
word by index. Note there are more powerful approaches should you need them.

.. ipython:: python

   firstlast = pd.DataFrame({"String": ["John Smith", "Jane Cook"]})
   firstlast["First_Name"] = firstlast["String"].str.split(" ", expand=True)[0]
   firstlast["Last_Name"] = firstlast["String"].str.rsplit(" ", expand=True)[0]
   firstlast
