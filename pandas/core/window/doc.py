from textwrap import dedent

from pandas.core.shared_docs import _shared_docs

_shared_docs = dict(**_shared_docs)

doc_template = dedent(
    """
    "Calculate the {window_method} {aggregation_description}"

    Parameters
    ----------
    {parameters}

    {numpy_args_kwargs}

    Returns
    -------
    Series or DataFrame
        Return type is the same as the original object.

    See Also
    --------
    pandas.Series.{window_method} : Calling {window_method} with Series data.
    pandas.DataFrame.{window_method} : Calling {window_method} with DataFrames.
    pandas.Series.{agg_method} : {agg_method} for Series.
    pandas.DataFrame.{agg_method} : {agg_method} for DataFrame.
    {other_see_also}
    """
)

numpy_args_kwargs = dedent(
    """
    *args, **kwargs
        For Numpy compatibility and has no effect on the computed value.
    """
)
