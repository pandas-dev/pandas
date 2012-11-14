""" Panel4D: a 4-d dict like collection of panels """

from pandas.core.panel import Panel
import pandas.lib as lib


class Panel4D(Panel):
    _AXIS_ORDERS  = ['labels','items','major_axis','minor_axis']
    _AXIS_NUMBERS = dict([ (a,i) for i, a in enumerate(_AXIS_ORDERS) ])
    _AXIS_ALIASES = {
        'major' : 'major_axis',
        'minor' : 'minor_axis'
    }
    _AXIS_NAMES   = dict([ (i,a) for i, a in enumerate(_AXIS_ORDERS) ])
    _AXIS_SLICEMAP = {
        'items'      : 'items',
        'major_axis' : 'major_axis',
        'minor_axis' : 'minor_axis'
        }
    _AXIS_LEN     = len(_AXIS_ORDERS)

    # major
    _default_stat_axis = 2

    # info axis
    _het_axis = 0
    _info_axis = _AXIS_ORDERS[_het_axis]

    labels     = lib.AxisProperty(0)
    items      = lib.AxisProperty(1)
    major_axis = lib.AxisProperty(2)
    minor_axis = lib.AxisProperty(3)

    _constructor_sliced = Panel

    def __init__(self, data=None, labels=None, items=None, major_axis=None, minor_axis=None, copy=False, dtype=None):
        """
        Represents a 4 dimensonal structured

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
        """
        self._init_data( data=data, labels=labels, items=items, major_axis=major_axis, minor_axis=minor_axis,
                         copy=copy, dtype=dtype)

    def _get_plane_axes(self, axis):
        axis = self._get_axis_name(axis)

        if axis == 'major_axis':
            items = self.labels
            major = self.items
            minor = self.minor_axis
        elif axis == 'minor_axis':
            items = self.labels
            major = self.items
            minor = self.major_axis
        elif axis == 'items':
            items = self.labels
            major = self.major_axis
            minor = self.minor_axis
        elif axis == 'labels':
            items = self.items
            major = self.major_axis
            minor = self.minor_axis

        return items, major, minor

    def _combine(self, other, func, axis=0):
        if isinstance(other, Panel4D):
            return self._combine_panel4d(other, func)
        return super(Panel4D, self)._combine(other, func, axis=axis)

    def _combine_panel4d(self, other, func):
        labels = self.labels + other.labels
        items  = self.items + other.items
        major  = self.major_axis + other.major_axis
        minor  = self.minor_axis + other.minor_axis

        # could check that everything's the same size, but forget it
        this = self.reindex(labels=labels, items=items, major=major, minor=minor)
        other = other.reindex(labels=labels, items=items, major=major, minor=minor)

        result_values = func(this.values, other.values)

        return self._constructor(result_values, labels, items, major, minor)

    def join(self, other, how='left', lsuffix='', rsuffix=''):
        if isinstance(other, Panel4D):
            join_major, join_minor = self._get_join_index(other, how)
            this = self.reindex(major=join_major, minor=join_minor)
            other = other.reindex(major=join_major, minor=join_minor)
            merged_data = this._data.merge(other._data, lsuffix, rsuffix)
            return self._constructor(merged_data)
        return super(Panel4D, self).join(other=other,how=how,lsuffix=lsuffix,rsuffix=rsuffix)

    ### remove operations ####
    def to_frame(self, *args, **kwargs):
        raise NotImplementedError
    def to_excel(self, *args, **kwargs):
        raise NotImplementedError

