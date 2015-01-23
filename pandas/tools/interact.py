from pandas.core.api import Int64Index, DataFrame, Series, Panel

try:  # tabview optional
    from tabview import tabview
except ImportError:
    pass


def interact(data, index=None, header=None, start=None,
             fixed_index=None, fixed_header=None, **kwargs):
    """
    Visualize the contents of ``data`` interactively.

    Parameters
    ----------
    data : DataFrame, Panel, Series or list
        Object containing the data
    index : bool
        Show the index. When None, detect if an index exists.
    header : bool
        Show the column names or Series name. When None, detect if column names
        or a Series name has been set.
    start : Y or (Y,X) tuple
        Start the viewer at the indicated location.
    fixed_index : bool
        Instruct the viewer to keep the index fixed
    fixed_header : bool
        Instruct the viewer to keep the header fixed
    **kwargs : dict
        Any parameter supported by the underlying viewer.
    """
    if isinstance(data, (Panel, DataFrame)):
        if isinstance(data, Panel):
            data = data.to_frame()

        # detect if data is using the built-in index for the labels/columns
        if index is None:
            index = type(data.index) is not Int64Index
        if header is None:
            header = type(data.index) is not Int64Index

        if index:
            data = data.reset_index()
        buf = []
        if header:
            buf += [data.columns.tolist()]
        buf += data.values.tolist()
        data = buf

    elif isinstance(data, Series):
        if index is None:
            index = type(data.index) is not Int64Index
        if header is None:
            header = data.name is not None and len(data.name)
        buf = []
        if index:
            if header:
                buf += [[data.index.name, data.name]]
            buf += data.reset_index().values.tolist()
        else:
            if header:
                buf += [[data.name]]
            buf += [[x] for x in data.values]
        data = buf

    else:
        # try to convert to a simple list
        data = list(data)

    # defaults
    if fixed_index is None:
        fixed_index = index
    if fixed_header is None:
        fixed_header = header
    if type(start) is not tuple:
        start = (start, 0)

    interact_list(data, **kwargs)


def interact_list(data, start=None, fixed_header=None,
                  fixed_index=None, **kwargs):
    tabview.view(data)
