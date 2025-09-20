"""Implement common widgets layouts as reusable components"""

import re
from collections import defaultdict

from traitlets import Instance, Bool, Unicode, CUnicode, CaselessStrEnum, Tuple
from traitlets import Integer
from traitlets import HasTraits, TraitError
from traitlets import observe, validate

from .widget import Widget
from .widget_box import GridBox

from .docutils import doc_subst


_doc_snippets = {
    'style_params' : """

    grid_gap : str
        CSS attribute used to set the gap between the grid cells

    justify_content : str, in ['flex-start', 'flex-end', 'center', 'space-between', 'space-around']
        CSS attribute used to align widgets vertically

    align_items : str, in ['top', 'bottom', 'center', 'flex-start', 'flex-end', 'baseline', 'stretch']
        CSS attribute used to align widgets horizontally

    width : str
    height : str
        width and height"""
    }

@doc_subst(_doc_snippets)
class LayoutProperties(HasTraits):
    """Mixin class for layout templates

    This class handles mainly style attributes (height, grid_gap etc.)

    Parameters
    ----------

    {style_params}


    Note
    ----

    This class is only meant to be used in inheritance as mixin with other
    classes. It will not work, unless `self.layout` attribute is defined.

    """

    # style attributes (passed to Layout)
    grid_gap = Unicode(
        None,
        allow_none=True,
        help="The grid-gap CSS attribute.")
    justify_content = CaselessStrEnum(
        ['flex-start', 'flex-end', 'center',
         'space-between', 'space-around'],
        allow_none=True,
        help="The justify-content CSS attribute.")
    align_items = CaselessStrEnum(
        ['top', 'bottom',
         'flex-start', 'flex-end', 'center',
         'baseline', 'stretch'],
        allow_none=True, help="The align-items CSS attribute.")
    width = Unicode(
        None,
        allow_none=True,
        help="The width CSS attribute.")
    height = Unicode(
        None,
        allow_none=True,
        help="The width CSS attribute.")


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._property_rewrite = defaultdict(dict)
        self._property_rewrite['align_items'] = {'top': 'flex-start',
                                                 'bottom': 'flex-end'}
        self._copy_layout_props()
        self._set_observers()

    def _delegate_to_layout(self, change):
        "delegate the trait types to their counterparts in self.layout"
        value, name = change['new'], change['name']
        value = self._property_rewrite[name].get(value, value)
        setattr(self.layout, name, value) # pylint: disable=no-member

    def _set_observers(self):
        "set observers on all layout properties defined in this class"
        _props = LayoutProperties.class_trait_names()
        self.observe(self._delegate_to_layout, _props)

    def _copy_layout_props(self):

        _props = LayoutProperties.class_trait_names()

        for prop in _props:
            value = getattr(self, prop)
            if value:
                value = self._property_rewrite[prop].get(value, value)
                setattr(self.layout, prop, value) #pylint: disable=no-member

@doc_subst(_doc_snippets)
class AppLayout(GridBox, LayoutProperties):
    """ Define an application like layout of widgets.

    Parameters
    ----------

    header: instance of Widget
    left_sidebar: instance of Widget
    center: instance of Widget
    right_sidebar: instance of Widget
    footer: instance of Widget
        widgets to fill the positions in the layout

    merge: bool
        flag to say whether the empty positions should be automatically merged

    pane_widths: list of numbers/strings
        the fraction of the total layout width each of the central panes should occupy
        (left_sidebar,
        center, right_sidebar)

    pane_heights: list of numbers/strings
        the fraction of the width the vertical space that the panes should occupy
         (left_sidebar, center, right_sidebar)

    {style_params}

    Examples
    --------

    """

    # widget positions
    header = Instance(Widget, allow_none=True)
    footer = Instance(Widget, allow_none=True)
    left_sidebar = Instance(Widget, allow_none=True)
    right_sidebar = Instance(Widget, allow_none=True)
    center = Instance(Widget, allow_none=True)

    # extra args
    pane_widths = Tuple(CUnicode(), CUnicode(), CUnicode(),
                        default_value=['1fr', '2fr', '1fr'])
    pane_heights = Tuple(CUnicode(), CUnicode(), CUnicode(),
                         default_value=['1fr', '3fr', '1fr'])

    merge = Bool(default_value=True)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._update_layout()

    @staticmethod
    def _size_to_css(size):
        if re.match(r'\d+\.?\d*(px|fr|%)$', size):
            return size
        if re.match(r'\d+\.?\d*$', size):
            return size + 'fr'

        raise TypeError("the pane sizes must be in one of the following formats: "
                        "'10px', '10fr', 10 (will be converted to '10fr')."
                        "Got '{}'".format(size))

    def _convert_sizes(self, size_list):
        return list(map(self._size_to_css, size_list))

    def _update_layout(self):

        grid_template_areas = [["header", "header", "header"],
                               ["left-sidebar", "center", "right-sidebar"],
                               ["footer", "footer", "footer"]]

        grid_template_columns = self._convert_sizes(self.pane_widths)
        grid_template_rows = self._convert_sizes(self.pane_heights)

        all_children = {'header': self.header,
                        'footer': self.footer,
                        'left-sidebar': self.left_sidebar,
                        'right-sidebar': self.right_sidebar,
                        'center': self.center}

        children = {position : child
                    for position, child in all_children.items()
                    if child is not None}

        if not children:
            return

        for position, child in children.items():
            child.layout.grid_area = position

        if self.merge:

            if len(children) == 1:
                position = list(children.keys())[0]
                grid_template_areas = [[position, position, position],
                                       [position, position, position],
                                       [position, position, position]]

            else:
                if self.center is None:
                    for row in grid_template_areas:
                        del row[1]
                    del grid_template_columns[1]

                if self.left_sidebar is None:
                    grid_template_areas[1][0] = grid_template_areas[1][1]

                if self.right_sidebar is None:
                    grid_template_areas[1][-1] = grid_template_areas[1][-2]

                if (self.left_sidebar is None and
                        self.right_sidebar is None and
                        self.center is None):
                    grid_template_areas = [['header'], ['footer']]
                    grid_template_columns = ['1fr']
                    grid_template_rows = ['1fr', '1fr']

                if self.header is None:
                    del grid_template_areas[0]
                    del grid_template_rows[0]

                if self.footer is None:
                    del grid_template_areas[-1]
                    del grid_template_rows[-1]


        grid_template_areas_css = "\n".join('"{}"'.format(" ".join(line))
                                            for line in grid_template_areas)

        self.layout.grid_template_columns = " ".join(grid_template_columns)
        self.layout.grid_template_rows = " ".join(grid_template_rows)
        self.layout.grid_template_areas = grid_template_areas_css

        self.children = tuple(children.values())

    @observe("footer", "header", "center", "left_sidebar", "right_sidebar", "merge",
             "pane_widths", "pane_heights")
    def _child_changed(self, change): #pylint: disable=unused-argument
        self._update_layout()


@doc_subst(_doc_snippets)
class GridspecLayout(GridBox, LayoutProperties):
    """ Define a N by M grid layout

    Parameters
    ----------

    n_rows : int
        number of rows in the grid

    n_columns : int
        number of columns in the grid

    {style_params}

    Examples
    --------

    >>> from ipywidgets import GridspecLayout, Button, Layout
    >>> layout = GridspecLayout(n_rows=4, n_columns=2, height='200px')
    >>> layout[:3, 0] = Button(layout=Layout(height='auto', width='auto'))
    >>> layout[1:, 1] = Button(layout=Layout(height='auto', width='auto'))
    >>> layout[-1, 0] = Button(layout=Layout(height='auto', width='auto'))
    >>> layout[0, 1] = Button(layout=Layout(height='auto', width='auto'))
    >>> layout
    """

    n_rows = Integer()
    n_columns = Integer()

    def __init__(self, n_rows=None, n_columns=None, **kwargs):
        super().__init__(**kwargs)
        self.n_rows = n_rows
        self.n_columns = n_columns
        self._grid_template_areas = [['.'] * self.n_columns for i in range(self.n_rows)]

        self._grid_template_rows = 'repeat(%d, 1fr)' % (self.n_rows,)
        self._grid_template_columns = 'repeat(%d, 1fr)' % (self.n_columns,)
        self._children = {}
        self._id_count = 0

    @validate('n_rows', 'n_columns')
    def _validate_integer(self, proposal):
        if proposal['value'] > 0:
            return proposal['value']
        raise TraitError('n_rows and n_columns must be positive integer')

    def _get_indices_from_slice(self, row, column):
        "convert a two-dimensional slice to a list of rows and column indices"

        if isinstance(row, slice):
            start, stop, stride = row.indices(self.n_rows)
            rows = range(start, stop, stride)
        else:
            rows = [row]

        if isinstance(column, slice):
            start, stop, stride = column.indices(self.n_columns)
            columns = range(start, stop, stride)
        else:
            columns = [column]

        return rows, columns

    def __setitem__(self, key, value):
        row, column = key
        self._id_count += 1
        obj_id = 'widget%03d' % self._id_count
        value.layout.grid_area = obj_id

        rows, columns = self._get_indices_from_slice(row, column)

        for row in rows:
            for column in columns:
                current_value = self._grid_template_areas[row][column]
                if current_value != '.' and current_value in self._children:
                    del self._children[current_value]
                self._grid_template_areas[row][column] = obj_id

        self._children[obj_id] = value
        self._update_layout()

    def __getitem__(self, key):
        rows, columns = self._get_indices_from_slice(*key)

        obj_id = None
        for row in rows:
            for column in columns:
                new_obj_id = self._grid_template_areas[row][column]
                obj_id = obj_id or new_obj_id
                if obj_id != new_obj_id:
                    raise TypeError('The slice spans several widgets, but '
                                    'only a single widget can be retrieved '
                                    'at a time')

        return self._children[obj_id]

    def _update_layout(self):

        grid_template_areas_css = "\n".join('"{}"'.format(" ".join(line))
                                            for line in self._grid_template_areas)

        self.layout.grid_template_columns = self._grid_template_columns
        self.layout.grid_template_rows = self._grid_template_rows
        self.layout.grid_template_areas = grid_template_areas_css
        self.children = tuple(self._children.values())


@doc_subst(_doc_snippets)
class TwoByTwoLayout(GridBox, LayoutProperties):
    """ Define a layout with 2x2 regular grid.

    Parameters
    ----------

    top_left: instance of Widget
    top_right: instance of Widget
    bottom_left: instance of Widget
    bottom_right: instance of Widget
        widgets to fill the positions in the layout

    merge: bool
        flag to say whether the empty positions should be automatically merged

    {style_params}

    Examples
    --------

    >>> from ipywidgets import TwoByTwoLayout, Button
    >>> TwoByTwoLayout(top_left=Button(description="Top left"),
    ...                top_right=Button(description="Top right"),
    ...                bottom_left=Button(description="Bottom left"),
    ...                bottom_right=Button(description="Bottom right"))

    """

    # widget positions
    top_left = Instance(Widget, allow_none=True)
    top_right = Instance(Widget, allow_none=True)
    bottom_left = Instance(Widget, allow_none=True)
    bottom_right = Instance(Widget, allow_none=True)

    # extra args
    merge = Bool(default_value=True)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._update_layout()

    def _update_layout(self):


        grid_template_areas = [["top-left", "top-right"],
                               ["bottom-left", "bottom-right"]]

        all_children = {'top-left' : self.top_left,
                        'top-right' : self.top_right,
                        'bottom-left' : self.bottom_left,
                        'bottom-right' : self.bottom_right}

        children = {position : child
                    for position, child in all_children.items()
                    if child is not None}

        if not children:
            return

        for position, child in children.items():
            child.layout.grid_area = position

        if self.merge:

            if len(children) == 1:
                position = list(children.keys())[0]
                grid_template_areas = [[position, position],
                                       [position, position]]
            else:
                columns = ['left', 'right']
                for i, column in enumerate(columns):
                    top, bottom = children.get('top-' + column), children.get('bottom-' + column)
                    i_neighbour = (i + 1) % 2
                    if top is None and bottom is None:
                        # merge each cell in this column with the neighbour on the same row
                        grid_template_areas[0][i] = grid_template_areas[0][i_neighbour]
                        grid_template_areas[1][i] = grid_template_areas[1][i_neighbour]
                    elif top is None:
                        # merge with the cell below
                        grid_template_areas[0][i] = grid_template_areas[1][i]
                    elif bottom is None:
                        # merge with the cell above
                        grid_template_areas[1][i] = grid_template_areas[0][i]

        grid_template_areas_css = "\n".join('"{}"'.format(" ".join(line))
                                            for line in grid_template_areas)

        self.layout.grid_template_columns = '1fr 1fr'
        self.layout.grid_template_rows = '1fr 1fr'
        self.layout.grid_template_areas = grid_template_areas_css

        self.children = tuple(children.values())

    @observe("top_left", "bottom_left", "top_right", "bottom_right", "merge")
    def _child_changed(self, change): #pylint: disable=unused-argument
        self._update_layout()
