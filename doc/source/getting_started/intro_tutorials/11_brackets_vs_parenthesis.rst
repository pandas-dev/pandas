In both Python and pandas, it’s important to understand the difference between parentheses ``()`` and square brackets ``[]``:

- **Parentheses** are used to call functions and methods. For example, ``df.mean()`` calculates the mean of a DataFrame, and parentheses are also used to group expressions, such as ``(a + b) * c``, or to create tuples: ``(1, 2, 3)``.
- **Square brackets** are used for indexing, selecting data, and defining lists. In pandas, you use square brackets to select columns or filter data—for instance, ``df['A']`` selects column ``A``, while ``df[0:3]`` selects rows by position. In standard Python, square brackets are also used to define a list, as in ``my_list = [1, 2, 3]``.

Remember:  
**Use [] for selection or indexing, and () for calling functions or grouping expressions.**
