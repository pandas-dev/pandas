Python extracts a substring from a string based on its text
by using regular expressions. There are much more powerful
approaches, but this just shows a simple approach.

.. ipython:: python

   firstlast = pd.DataFrame({"String": ["John Smith", "Jane Cook"]})
   firstlast["First_Name"] = firstlast["String"].str.split(" ", expand=True)[0]
   firstlast["Last_Name"] = firstlast["String"].str.rsplit(" ", expand=True)[0]
   firstlast
