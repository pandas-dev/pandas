from __future__ import annotations
from pandas.core.frame import DataFrame
import pandas as pd
import numpy as np
from itertools import cycle


def generate_fake_dataframe(
        size : int,
        cols : str,
        col_names : list | None = None,
        intervals: dict | list | None = None,
        seed: int | None = None
) -> DataFrame:
    """
    Generate a DataFrame with dummy data in a controlled manner.
    
    Supported datatypes for the columns are int, float, date and categorical.
    For categorical columns, we can choose to get data from these families:
    'names', 'animals', 'cities' and 'colors'.

    Parameters
    ----------
    size : int
        The number of rows in the DataFrame.
    cols : str
        String whose characters represent the types of the columns generated.
        'i' = int, 'f' = float, 'd' = date and 'c' = categorical.
        i.e cols = "iifccd" (2 int cols, 1 float, 2 categorical and 1 date)
    col_names : list or None
        List of column names for the generated DataFrame.
        Must be of the same length as cols.
        If None, column names are assigned as 'column_n_dtype' where n is an
        integer representing the nth column.
    intervals: dict, list or None
        Intervals define how the data generated. 'i', 'f' and 'd' intervals are set
        as tuples of the form (start,end), where data is sampled between those bounds.
        'c' intervals can be defined as:
            a) a tuple of the form (family, n) where n is the max number of different
               elements from the family used to populate the column.
            b) a list of objects used to populate the column.
        The default intervals set are (0, 10) for integers,
        (0,100) for floats,  ('2020-01-01', '2020-12-31') for dates and
        ('names', 5) for categorical columns.
        `intervals` can receive a `dict` of the form: {char : interval ...}
        which updates the default intervals for specific column types.
        i.e intervals = {'i' : (-5,5), 'd' : ("1999-01-01","2021-12-05")}
        A list can be passed instead with intervals for each column. This allows control
        over how each column is generated (this list must have the same length as `cols`
        i.e cols = "fcc", intervals = [(2.5, 9.5), ('colors', 3), ['custom', 'family']].
        If intervals = None, the default intervals are used. None can also be placed
        inside the `intervals` list to use the default interval for that column type.
        i.e cols = "ifc", intervals = [(2,3), None, ('colors', 3)]

    seed: int
        Seed used for the random sampling. If the same parameters and seed are used,
        the generated dataframe is guaranteed to be the same.

    Returns
    -------
    DataFrame

    See Also
    --------
    makeDataFrame() : Generates a (30, 4) DataFrame with random float values.
    makeMissingDataframe() : Generates a (30, 4) DataFrame with some NANs.
    makeMixedDataFrame() : Generates a predefined (5, 4) DataFrame.

    Examples
    --------
    >>> pd.generate_fake_dataframe(1000, "iiccd")
    >>> pd.generate_fake_dataframe(size = 31,
                           cols = "cid",
                           col_names = ["student", "grade", "birthday"])
    >>> pd.generate_fake_dataframe(20, "cccd", intervals = {"c": ("colors",7)})
    >>> pd.generate_fake_dataframe(
                        size = 10,
                        cols = "cccfd",
                        col_names = ["name", "pet", "city","height", "birthday"],
                        intervals = {"f" : (1.72,1.95),
                                     "d" : ("1996-01-01","1996-12-31")},
                        seed = 1)
    >>> pd.generate_fake_dataframe(250,
                        cols = "cncccd",
                        intervals = [("names", 1),
                                     (2.4,3.4),
                                     ["my", "custom", "list", "of", "categories"],
                                     [True,False],
                                     None
                                    ],
                        seed = 42)
    >>> pd.generate_fake_dataframe(
                        size = 10,
                        cols = "cccfd",
                        col_names=["name", "pet", "city","height", "birthday"],
                        intervals = {"f" : (1.72,1.95),
                                     "d" : ("1996-01-01","1996-12-31")},
                        )
    >>> pd.generate_fake_dataframe(
                        size = 222,
                        cols = "ccffd",
                        intervals = {"c" : [1,0, "-1", 3.14, np.e],
                                    "d" : ("1996-01-01","1996-12-31")}
                        )
    >>> pd.generate_fake_dataframe(
                        3000,
                        "ififi",
                        intervals = [None, None, (-5,5), None, (-1.2,2.4), None]
                        )
    """

    categories_dict = {
        'animals': [
            'cow',
            'rabbit',
            'duck',
            'shrimp',
            'pig',
            'goat',
            'crab',
            'deer',
            'bee',
            'sheep',
            'fish',
            'turkey',
            'dove',
            'chicken',
            'horse'],
        'names' : [
            'James',
            'Mary',
            'Robert',
            'Patricia',
            'John',
            'Jennifer',
            'Michael',
            'Linda',
            'William',
            'Elizabeth',
            'Ahmed',
            'Barbara',
            'Richard',
            'Susan',
            'Salomon',
            'Juan Luis'],
        'cities' : [
            'Stockholm',
            'Denver',
            'Moscow',
            'Marseille',
            'Palermo',
            'Tokyo',
            'Lisbon',
            'Oslo',
            'Nairobi',
            'Río de Janeiro',
            'Berlin',
            'Bogotá',
            'Manila',
            'Madrid',
            'Milwaukee'],
        'colors' : [
            'red',
            'orange',
            'yellow',
            'green',
            'blue',
            'indigo',
            'purple',
            'pink',
            'silver',
            'gold',
            'beige',
            'brown',
            'grey',
            'black',
            'white']}
    default_intervals = {
        "i" : (0, 10),
        "f" : (0, 100),
        "c" : ("names", 5),
        "d" : ("2020-01-01", "2020-12-31")}

    rng = np.random.default_rng(seed)

    first_c = default_intervals["c"][0]
    categories_names = cycle([first_c] +
                             [c for c in categories_dict.keys() if c != first_c])
    default_intervals["c"] = (categories_names, default_intervals["c"][1])

    if isinstance(col_names, list):
        assert len(col_names) == len(cols), \
            f"The fake DataFrame should have {len(cols)} columns \
            but col_names is a list with {len(col_names)} elements"
    elif col_names is None:
        suffix = {"c" : "cat", "i" : "int", "f" : "float", "d" : "date"}
        col_names = [f"column_{str(i)}_{suffix.get(col)}" for i, col in enumerate(cols)]

    if isinstance(intervals, list):
        assert len(intervals) == len(cols), \
            f"The fake DataFrame should have {len(cols)} columns \
            but intervals is a list with {len(intervals)} elements"
    else:
        if isinstance(intervals, dict):
            assert len(set(intervals.keys()) - set(default_intervals.keys())) == 0, \
                "The intervals parameter has invalid keys"
            default_intervals.update(intervals)
        intervals = [default_intervals[col] for col in cols]

    df = DataFrame()
    for col, col_name, interval in zip(cols, col_names, intervals):
        if interval is None:
            interval = default_intervals[col]
        assert (len(interval) == 2 and isinstance(interval, tuple)) \
            or isinstance(interval, list), \
            f"This interval {interval} is neither a tuple of two elements nor a list."
        if col in ("i", "f"):
            start, end = interval
            if col == "i":
                df[col_name] = rng.integers(start, end, size)
            elif col == "f":
                df[col_name] = rng.uniform(start, end, size)
        elif col == "c":
            if isinstance(interval, list):
                categories = np.array(interval)
            else:
                cat_family, length = interval
                if isinstance(cat_family, cycle):
                    cat_family = next(cat_family)
                assert cat_family in categories_dict.keys(), \
                    f"There are no categories for family '{cat_family}'. \
                    Consider passing a list with categories \
                    or use one of the available families: {categories_dict.keys()}"
                categories = rng.choice(
                    categories_dict[cat_family],
                    length,
                    replace=False,
                    shuffle=True)
            df[col_name] = rng.choice(categories, size, shuffle=True)
        elif col == "d":
            start, end = interval
            df[col_name] = rng.choice(pd.date_range(start, end), size)
    return df
