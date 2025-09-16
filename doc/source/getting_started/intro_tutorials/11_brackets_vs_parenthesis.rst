.. _10min_tut_11_brackets_vs_parenthesis:

Parentheses vs. Square Brackets in Python and pandas
====================================================

In Python and pandas, parentheses ``()`` and square brackets ``[]`` have distinct uses, which can sometimes be confusing for new users. Parentheses are used to call functions and methods (for example, ``df.mean()``), to group expressions in calculations (such as ``(a + b) * c``), and to create tuples. Square brackets, on the other hand, are used for defining lists and for indexing or selecting data—in both core Python and pandas. For example, ``df['A']`` selects column ``A`` from a DataFrame, and ``df[0:3]`` selects rows by position. In pandas, square brackets are always used when you want to select or filter data, while parentheses are used any time you are calling a method or function. Remember: use ``[]`` for selection or indexing, and ``()`` for function calls and grouping.
