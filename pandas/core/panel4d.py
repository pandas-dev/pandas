""" Panel4D: a 4-d dict like collection of panels """

import warnings
from pandas.core.generic import NDFrame
from pandas.core.panelnd import create_nd_panel_factory
from pandas.core.panel import Panel
from pandas.util._validators import validate_axis_style_args


Panel4D = create_nd_panel_factory(klass_name='Panel4D',
                                  orders=['labels', 'items', 'major_axis',
                                          'minor_axis'],
                                  slices={'labels': 'labels',
                                          'items': 'items',
                                          'major_axis': 'major_axis',
                                          'minor_axis': 'minor_axis'},
                                  slicer=Panel,
                                  aliases={'major': 'major_axis',
                                           'minor': 'minor_axis'}, stat_axis=2,
                                  ns=dict(__doc__="""
    Panel4D is a 4-Dimensional named container very much like a Panel, but
    having 4 named dimensions. It is intended as a test bed for more
    N-Dimensional named containers.

    .. deprecated:: 0.19.0
       The recommended way to represent these types of n-dimensional data
       are with the `xarray package <http://xarray.pydata.org/en/stable/>`__.
       Pandas provides a `.to_xarray()` method to automate this conversion.

    Parameters
    ----------
    data : ndarray (labels x items x major x minor), or dict of Panels

    labels : Index or array-like : axis=0
    items  : Index or array-like : axis=1
    major_axis : Index or array-like: axis=2
    minor_axis : Index or array-like: axis=3

    dtype : dtype, default None
    Data type to force, otherwise infer
    copy : boolean, default False
    Copy data from inputs. Only affects DataFrame / 2d ndarray input
    """))


def panel4d_init(self, data=None, labels=None, items=None, major_axis=None,
                 minor_axis=None, copy=False, dtype=None):

    # deprecation GH13564
    warnings.warn("\nPanel4D is deprecated and will be removed in a "
                  "future version.\nThe recommended way to represent "
                  "these types of n-dimensional data are with\n"
                  "the `xarray package "
                  "<http://xarray.pydata.org/en/stable/>`__.\n"
                  "Pandas provides a `.to_xarray()` method to help "
                  "automate this conversion.\n",
                  FutureWarning, stacklevel=2)
    self._init_data(data=data, labels=labels, items=items,
                    major_axis=major_axis, minor_axis=minor_axis, copy=copy,
                    dtype=dtype)


def panel4d_reindex(self, labs=None, labels=None, items=None, major_axis=None,
                    minor_axis=None, axis=None, **kwargs):
    # Hack for reindex_axis deprecation
    # Ha, we used labels for two different things
    # I think this will work still.
    if labs is None:
        args = ()
    else:
        args = (labs,)
    kwargs_ = dict(labels=labels,
                   items=items,
                   major_axis=major_axis,
                   minor_axis=minor_axis,
                   axis=axis)
    kwargs_ = {k: v for k, v in kwargs_.items() if v is not None}
    # major = kwargs.pop("major", None)
    # minor = kwargs.pop('minor', None)

    # if major is not None:
    #     if kwargs.get("major_axis"):
    #         raise TypeError("Cannot specify both 'major' and 'major_axis'")
    #     kwargs_['major_axis'] = major
    # if minor is not None:
    #     if kwargs.get("minor_axis"):
    #         raise TypeError("Cannot specify both 'minor' and 'minor_axis'")
    #     kwargs_['minor_axis'] = minor

    if axis is not None:
        kwargs_['axis'] = axis

    axes = validate_axis_style_args(self, args, kwargs_, 'labs', 'reindex')
    kwargs.update(axes)
    return NDFrame.reindex(self, **kwargs)


Panel4D.__init__ = panel4d_init
Panel4D.reindex = panel4d_reindex
