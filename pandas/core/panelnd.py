""" Factory methods to create N-D panels """

import pandas
from pandas.core.panel import Panel
import pandas.lib as lib

def create_nd_panel_factory(klass_name, axis_orders, axis_slices, slicer, axis_aliases = None, stat_axis = 2):
    """ manufacture a n-d class:

        parameters
        ----------
        klass_name  : the klass name
        axis_orders : the names of the axes in order (highest to lowest)
        axis_slices : a dictionary that defines how the axes map to the sliced axis
        slicer      : the class representing a slice of this panel
        axis_aliases: a dictionary defining aliases for various axes 
                        default = { major : major_axis, minor : minor_axis }
        stat_axis   : the default statistic axis
                        default = 2
        het_axis    : the info axis


        returns
        -------
        a class object reprsenting this panel


    """

    # build the klass
    klass = type(klass_name, (slicer,),{}) 

    # add the class variables
    klass._AXIS_ORDERS   = axis_orders
    klass._AXIS_NUMBERS  = dict([ (a,i) for i, a in enumerate(axis_orders) ])
    klass._AXIS_ALIASES  = axis_aliases or dict()
    klass._AXIS_NAMES    = dict([ (i,a) for i, a in enumerate(axis_orders) ])
    klass._AXIS_SLICEMAP = axis_slices
    klass._AXIS_LEN      = len(axis_orders)
    klass._default_stat_axis = stat_axis
    klass._het_axis      = 0
    klass._info_axis     = axis_orders[klass._het_axis]
    klass._constructor_sliced = slicer

    # add the axes
    for i, a in enumerate(axis_orders):
        setattr(klass,a,lib.AxisProperty(i))

    #### define the methods ####
    def __init__(self, *args, **kwargs):
        if not (kwargs.get('data') or len(args)):
            raise Exception("must supply at least a data argument to [%s]" % klass_name)
        if 'copy' not in kwargs:
            kwargs['copy'] = False
        if 'dtype' not in kwargs:
            kwargs['dtype'] = None
        self._init_data( *args, **kwargs)
    klass.__init__ = __init__

    def _get_plane_axes(self, axis):

        axis   = self._get_axis_name(axis)
        index  = self._AXIS_ORDERS.index(axis)

        planes = []
        if index:
            planes.extend(self._AXIS_ORDERS[0:index])
        if index != self._AXIS_LEN:
            planes.extend(self._AXIS_ORDERS[index+1:])

        return [ getattr(self,p) for p in planes ]
    klass._get_plane_axes = _get_plane_axes

    def _combine(self, other, func, axis=0):
        if isinstance(other, klass):
            return self._combine_with_constructor(other, func)
        return super(klass, self)._combine(other, func, axis=axis)
    klass._combine = _combine

    def _combine_with_constructor(self, other, func):

        # combine labels to form new axes
        new_axes = []
        for a in self._AXIS_ORDERS:
            new_axes.append(getattr(self,a) + getattr(other,a))

        # reindex: could check that everything's the same size, but forget it
        d = dict([ (a,ax) for a,ax in zip(self._AXIS_ORDERS,new_axes) ])
        d['copy'] = False
        this = self.reindex(**d)
        other = other.reindex(**d)
    
        result_values = func(this.values, other.values)

        return self._constructor(result_values, **d)
    klass._combine_with_constructor = _combine_with_constructor

    # set as NonImplemented operations which we don't support
    for f in ['to_frame','to_excel','to_sparse','groupby','join','_get_join_index']:
        def func(self, *args, **kwargs):
            raise NotImplementedError
        setattr(klass,f,func)

    return klass


if __name__ == '__main__':

    # create a sample
    from pandas.util import testing
    print pandas.__version__

    # create a 4D
    Panel4DNew = create_nd_panel_factory(
        klass_name   = 'Panel4DNew', 
        axis_orders  = ['labels1','items1','major_axis','minor_axis'], 
        axis_slices  = { 'items1' : 'items', 'major_axis' : 'major_axis', 'minor_axis' : 'minor_axis' },
        slicer       = Panel,
        axis_aliases = { 'major' : 'major_axis', 'minor' : 'minor_axis' },
        stat_axis    = 2)
    
    p4dn = Panel4DNew(dict(L1 = testing.makePanel(), L2 = testing.makePanel()))
    print "creating a 4-D Panel"
    print p4dn, "\n"

    # create a 5D
    Panel5DNew = create_nd_panel_factory(
        klass_name   = 'Panel5DNew', 
        axis_orders  = [ 'cool1', 'labels1','items1','major_axis','minor_axis'], 
        axis_slices  = { 'labels1' : 'labels1', 'items1' : 'items', 'major_axis' : 'major_axis', 'minor_axis' : 'minor_axis' },
        slicer       = Panel4DNew,
        axis_aliases = { 'major' : 'major_axis', 'minor' : 'minor_axis' },
        stat_axis    = 2)
    
    p5dn = Panel5DNew(dict(C1 = p4dn))

    print "creating a 5-D Panel"
    print p5dn, "\n"

    print "Slicing p5dn"
    print p5dn.ix['C1',:,:,0:3,:], "\n"

    print "Transposing p5dn"
    print p5dn.transpose(1,2,3,4,0), "\n"
