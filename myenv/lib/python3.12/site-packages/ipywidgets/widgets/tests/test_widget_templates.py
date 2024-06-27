"Testing widget layout templates"

from unittest import TestCase
from unittest import mock

import pytest

import traitlets
import ipywidgets as widgets
from ipywidgets.widgets.widget_templates import LayoutProperties

class TestTwoByTwoLayout(TestCase):
    """test layout templates"""

    def test_merge_cells(self): #pylint: disable=no-self-use
        """test merging cells with missing widgets"""

        button1 = widgets.Button()
        button2 = widgets.Button()
        button3 = widgets.Button()
        button4 = widgets.Button()

        box = widgets.TwoByTwoLayout(top_left=button1,
                                     top_right=button2,
                                     bottom_left=button3,
                                     bottom_right=button4)

        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"bottom-left bottom-right"')
        assert box.top_left.layout.grid_area == 'top-left'
        assert box.top_right.layout.grid_area == 'top-right'
        assert box.bottom_left.layout.grid_area == 'bottom-left'
        assert box.bottom_right.layout.grid_area == 'bottom-right'
        assert len(box.get_state()['children']) == 4

        box = widgets.TwoByTwoLayout(top_left=button1,
                                     top_right=button2,
                                     bottom_left=None,
                                     bottom_right=button4)

        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"top-left bottom-right"')
        assert box.top_left.layout.grid_area == 'top-left'
        assert box.top_right.layout.grid_area == 'top-right'
        assert box.bottom_left is None
        assert box.bottom_right.layout.grid_area == 'bottom-right'
        assert len(box.get_state()['children']) == 3

        box = widgets.TwoByTwoLayout(top_left=None,
                                     top_right=button2,
                                     bottom_left=button3,
                                     bottom_right=button4)

        assert box.layout.grid_template_areas == ('"bottom-left top-right"\n' +
                                                  '"bottom-left bottom-right"')
        assert box.top_left is None
        assert box.top_right.layout.grid_area == 'top-right'
        assert box.bottom_left.layout.grid_area == 'bottom-left'
        assert box.bottom_right.layout.grid_area == 'bottom-right'
        assert len(box.get_state()['children']) == 3

        box = widgets.TwoByTwoLayout(top_left=None,
                                     top_right=button2,
                                     bottom_left=None,
                                     bottom_right=button4)

        assert box.layout.grid_template_areas == ('"top-right top-right"\n' +
                                                  '"bottom-right bottom-right"')
        assert box.top_left is None
        assert box.top_right.layout.grid_area == 'top-right'
        assert box.bottom_left is None
        assert box.bottom_right.layout.grid_area == 'bottom-right'
        assert len(box.get_state()['children']) == 2

        box = widgets.TwoByTwoLayout(top_left=button1,
                                     top_right=None,
                                     bottom_left=button3,
                                     bottom_right=button4)

        assert box.layout.grid_template_areas == ('"top-left bottom-right"\n' +
                                                  '"bottom-left bottom-right"')
        assert box.top_left.layout.grid_area == 'top-left'
        assert box.top_right is None
        assert box.bottom_left.layout.grid_area == 'bottom-left'
        assert box.bottom_right.layout.grid_area == 'bottom-right'
        assert len(box.get_state()['children']) == 3


        box = widgets.TwoByTwoLayout(top_left=button1,
                                     top_right=None,
                                     bottom_left=None,
                                     bottom_right=None)

        assert box.layout.grid_template_areas == ('"top-left top-left"\n' +
                                                  '"top-left top-left"')

        assert box.top_left is button1
        assert box.top_left.layout.grid_area == 'top-left'
        assert box.top_right is None
        assert box.bottom_left is None
        assert box.bottom_right is None
        assert len(box.get_state()['children']) == 1

        box = widgets.TwoByTwoLayout(top_left=None,
                                     top_right=button1,
                                     bottom_left=None,
                                     bottom_right=None)

        assert box.layout.grid_template_areas == ('"top-right top-right"\n' +
                                                  '"top-right top-right"')

        assert box.top_right is button1
        assert box.top_right.layout.grid_area == 'top-right'
        assert box.top_left is None
        assert box.bottom_left is None
        assert box.bottom_right is None
        assert len(box.get_state()['children']) == 1

        box = widgets.TwoByTwoLayout(top_left=None,
                                     top_right=None,
                                     bottom_left=None,
                                     bottom_right=None)

        assert box.layout.grid_template_areas is None
        assert box.top_left is None
        assert box.top_right is None
        assert box.bottom_left is None
        assert box.bottom_right is None
        assert not box.get_state()['children']

        box = widgets.TwoByTwoLayout(top_left=None,
                                     top_right=button1,
                                     bottom_left=None,
                                     bottom_right=None,
                                     merge=False)

        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"bottom-left bottom-right"')

        assert box.top_right is button1
        assert box.top_right.layout.grid_area == 'top-right'
        assert box.top_left is None
        assert box.bottom_left is None
        assert box.bottom_right is None
        assert len(box.get_state()['children']) == 1

    def test_keep_layout_options(self): #pylint: disable=no-self-use
        """test whether layout options are passed down to GridBox"""

        layout = widgets.Layout(align_items="center")
        button1 = widgets.Button()
        button2 = widgets.Button()
        button3 = widgets.Button()
        button4 = widgets.Button()

        box = widgets.TwoByTwoLayout(top_left=button1, top_right=button2,
                                     bottom_left=button3, bottom_right=button4,
                                     layout=layout)

        assert box.layout.align_items == 'center'

    def test_pass_layout_options(self): #pylint: disable=no-self-use
        """test whether the extra layout options of the template class are
           passed down to Layout object"""

        button1 = widgets.Button()
        button2 = widgets.Button()
        button3 = widgets.Button()
        button4 = widgets.Button()

        box = widgets.TwoByTwoLayout(top_left=button1, top_right=button2,
                                     bottom_left=button3, bottom_right=button4,
                                     grid_gap="10px", justify_content="center",
                                     align_items="center")

        assert box.layout.grid_gap == "10px"
        assert box.layout.justify_content == "center"
        assert box.layout.align_items == "center"

        # we still should be able to pass them through layout
        layout = widgets.Layout(grid_gap="10px", justify_content="center",
                                align_items="center")
        box = widgets.TwoByTwoLayout(top_left=button1, top_right=button2,
                                     bottom_left=button3, bottom_right=button4,
                                     layout=layout
                                     )

        assert box.layout.grid_gap == "10px"
        assert box.layout.justify_content == "center"
        assert box.layout.align_items == "center"

        # values passed directly in the constructor should overwrite layout options
        layout = widgets.Layout(grid_gap="10px", justify_content="center",
                                align_items="center")
        box = widgets.TwoByTwoLayout(top_left=button1, top_right=button2,
                                     bottom_left=button3, bottom_right=button4,
                                     layout=layout, grid_gap="30px"
                                     )

        assert box.layout.grid_gap == "30px"
        assert box.layout.justify_content == "center"
        assert box.layout.align_items == "center"


    @mock.patch("ipywidgets.Layout.send_state")
    def test_update_dynamically(self, send_state): #pylint: disable=no-self-use
        """test whether it's possible to add widget outside __init__"""

        button1 = widgets.Button()
        button2 = widgets.Button()
        button3 = widgets.Button()
        button4 = widgets.Button()

        box = widgets.TwoByTwoLayout(top_left=button1, top_right=button3,
                                     bottom_left=None, bottom_right=button4)
        from ipykernel.kernelbase import Kernel

        state = box.get_state()
        assert len(state['children']) == 3
        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"top-left bottom-right"')

        box.layout.comm.kernel = mock.MagicMock(spec=Kernel) #for mocking purposes
        send_state.reset_mock()
        box.bottom_left = button2

        state = box.get_state()
        assert len(state['children']) == 4
        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"bottom-left bottom-right"')
        # check whether frontend was informed
        send_state.assert_called_with(key="grid_template_areas")

        box = widgets.TwoByTwoLayout(top_left=button1, top_right=button3,
                                     bottom_left=None, bottom_right=button4)
        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"top-left bottom-right"')
        box.layout.comm.kernel = mock.MagicMock(spec=Kernel) #for mocking purposes
        send_state.reset_mock()
        box.merge = False
        assert box.layout.grid_template_areas == ('"top-left top-right"\n' +
                                                  '"bottom-left bottom-right"')
        send_state.assert_called_with(key="grid_template_areas")


class TestAppLayout(TestCase):
    """test layout templates"""

    def test_create_with_defaults(self):
        "test creating with default values"

        footer = widgets.Button()
        header = widgets.Button()
        center = widgets.Button()
        left_sidebar = widgets.Button()
        right_sidebar = widgets.Button()

        box = widgets.AppLayout(
            footer=footer,
            header=header,
            center=center,
            left_sidebar=left_sidebar,
            right_sidebar=right_sidebar
        )

        assert box.layout.grid_template_areas == ('"header header header"\n' +
                                                  '"left-sidebar center right-sidebar"\n' +
                                                  '"footer footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.center.layout.grid_area == 'center'
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'

        assert len(box.get_state()['children']) == 5

        # empty layout should produce no effects

        box = widgets.AppLayout()
        assert box.layout.grid_template_areas is None
        assert box.layout.grid_template_columns is None
        assert box.layout.grid_template_rows is None
        assert len(box.get_state()['children']) == 0


    def test_merge_empty_cells(self):
        "test if cells are correctly merged"

        footer = widgets.Button()
        header = widgets.Button()
        center = widgets.Button()
        left_sidebar = widgets.Button()
        right_sidebar = widgets.Button()

        # merge all if only one widget
        box = widgets.AppLayout(
            center=center
        )

        assert box.layout.grid_template_areas == ('"center center center"\n' +
                                                  '"center center center"\n' +
                                                  '"center center center"')
        assert box.center.layout.grid_area == 'center'

        assert len(box.get_state()['children']) == 1

        box = widgets.AppLayout(
            left_sidebar=left_sidebar
        )

        assert box.layout.grid_template_areas == ('"left-sidebar left-sidebar left-sidebar"\n' +
                                                  '"left-sidebar left-sidebar left-sidebar"\n' +
                                                  '"left-sidebar left-sidebar left-sidebar"')
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'

        assert len(box.get_state()['children']) == 1

        # merge left and right sidebars with center

        box = widgets.AppLayout(
            header=header,
            footer=footer,
            left_sidebar=left_sidebar,
            center=center
        )

        assert box.layout.grid_template_areas == ('"header header header"\n' +
                                                  '"left-sidebar center center"\n' +
                                                  '"footer footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.center.layout.grid_area == 'center'
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'
        assert len(box.get_state()['children']) == 4

        box = widgets.AppLayout(
            header=header,
            footer=footer,
            right_sidebar=right_sidebar,
            center=center
        )

        assert box.layout.grid_template_areas == ('"header header header"\n' +
                                                  '"center center right-sidebar"\n' +
                                                  '"footer footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.center.layout.grid_area == 'center'
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'
        assert len(box.get_state()['children']) == 4

        box = widgets.AppLayout(
            header=header,
            footer=footer,
            center=center
        )

        assert box.layout.grid_template_areas == ('"header header header"\n' +
                                                  '"center center center"\n' +
                                                  '"footer footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.center.layout.grid_area == 'center'
        assert len(box.get_state()['children']) == 3

        # if only center missing, remove it from view
        box = widgets.AppLayout(
            header=header,
            footer=footer,
            center=None,
            left_sidebar=left_sidebar,
            right_sidebar=right_sidebar
        )

        assert box.layout.grid_template_areas == ('"header header"\n' +
                                                  '"left-sidebar right-sidebar"\n' +
                                                  '"footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'
        assert box.center is None
        assert len(box.get_state()['children']) == 4

        # center and one sidebar missing -> 3 row arrangement
        box = widgets.AppLayout(
            header=header,
            footer=footer,
            center=None,
            left_sidebar=None,
            right_sidebar=right_sidebar
        )

        assert box.layout.grid_template_areas == ('"header header"\n' +
                                                  '"right-sidebar right-sidebar"\n' +
                                                  '"footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.left_sidebar is None
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'
        assert box.center is None
        assert len(box.get_state()['children']) == 3


        # remove middle row is both sidebars and center missing
        box = widgets.AppLayout(
            header=header,
            footer=footer,
            center=None,
            left_sidebar=None,
            right_sidebar=None
        )

        assert box.layout.grid_template_areas == ('"header"\n' +
                                                  '"footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.center is None
        assert box.left_sidebar is None
        assert box.right_sidebar is None
        assert len(box.get_state()['children']) == 2



        # do not merge if merge=False
        box = widgets.AppLayout(
            header=header,
            footer=footer,
            center=center,
            merge=False
        )

        assert box.layout.grid_template_areas == ('"header header header"\n' +
                                                  '"left-sidebar center right-sidebar"\n' +
                                                  '"footer footer footer"')
        assert box.footer.layout.grid_area == 'footer'
        assert box.header.layout.grid_area == 'header'
        assert box.center.layout.grid_area == 'center'
        assert box.left_sidebar is None
        assert box.right_sidebar is None
        assert len(box.get_state()['children']) == 3

        # merge header and footer simply removes it from view
        box = widgets.AppLayout(
            footer=footer,
            center=center,
            left_sidebar=left_sidebar,
            right_sidebar=right_sidebar
        )

        assert box.layout.grid_template_areas == ('"left-sidebar center right-sidebar"\n' +
                                                  '"footer footer footer"')
        assert box.center.layout.grid_area == 'center'
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'
        assert box.footer.layout.grid_area == 'footer'
        assert box.header is None
        assert len(box.get_state()['children']) == 4

        box = widgets.AppLayout(
            header=header,
            center=center,
            left_sidebar=left_sidebar,
            right_sidebar=right_sidebar
        )

        assert box.layout.grid_template_areas == ('"header header header"\n' +
                                                  '"left-sidebar center right-sidebar"')
        assert box.center.layout.grid_area == 'center'
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'
        assert box.header.layout.grid_area == 'header'
        assert box.footer is None
        assert len(box.get_state()['children']) == 4

        box = widgets.AppLayout(
            center=center,
            left_sidebar=left_sidebar,
            right_sidebar=right_sidebar
        )

        assert box.layout.grid_template_areas == '"left-sidebar center right-sidebar"'
        assert box.center.layout.grid_area == 'center'
        assert box.left_sidebar.layout.grid_area == 'left-sidebar'
        assert box.right_sidebar.layout.grid_area == 'right-sidebar'
        assert box.footer is None
        assert box.header is None
        assert len(box.get_state()['children']) == 3

        # merge all if only one widget
        box = widgets.AppLayout(
            center=center
        )

        assert box.layout.grid_template_areas == ('"center center center"\n' +
                                                  '"center center center"\n' +
                                                  '"center center center"')
        assert box.center.layout.grid_area == 'center'

        assert len(box.get_state()['children']) == 1

    def test_size_to_css(self):

        box = widgets.AppLayout()
        assert box._size_to_css("100px") == '100px'
        assert box._size_to_css("1fr") == '1fr'
        assert box._size_to_css("2.5fr") == '2.5fr'
        assert box._size_to_css('2.5') == '2.5fr'
        assert box._size_to_css('25%') == '25%'

        with pytest.raises(TypeError):
            box._size_to_css('this is not correct size')


    def test_set_pane_widths_heights(self):

        footer = widgets.Button()
        header = widgets.Button()
        center = widgets.Button()
        left_sidebar = widgets.Button()
        right_sidebar = widgets.Button()

        box = widgets.AppLayout(
            header=header,
            footer=footer,
            left_sidebar=left_sidebar,
            right_sidebar=left_sidebar,
            center=center
        )

        with pytest.raises(traitlets.TraitError):
            box.pane_widths = ['1fx', '1fx', '1fx', '1fx']
        with pytest.raises(traitlets.TraitError):
            box.pane_widths = ['1fx', '1fx']

        with pytest.raises(traitlets.TraitError):
            box.pane_heights = ['1fx', '1fx', '1fx', '1fx']
        with pytest.raises(traitlets.TraitError):
            box.pane_heights = ['1fx', '1fx']

        assert box.layout.grid_template_rows == "1fr 3fr 1fr"
        assert box.layout.grid_template_columns == "1fr 2fr 1fr"

        box.pane_heights = ['3fr', '100px', 20]
        assert box.layout.grid_template_rows == "3fr 100px 20fr"
        assert box.layout.grid_template_columns == "1fr 2fr 1fr"

        box.pane_widths = [3, 3, 1]
        assert box.layout.grid_template_rows == "3fr 100px 20fr"
        assert box.layout.grid_template_columns == "3fr 3fr 1fr"

class TestGridspecLayout(TestCase):
    "test GridspecLayout"

    def test_init(self):
        with pytest.raises(traitlets.TraitError):
            box = widgets.GridspecLayout()

        with pytest.raises(traitlets.TraitError):
            box = widgets.GridspecLayout(n_rows=-1, n_columns=1)

        box = widgets.GridspecLayout(n_rows=5, n_columns=3)
        assert box.n_rows == 5
        assert box.n_columns == 3
        assert len(box._grid_template_areas) == 5
        assert len(box._grid_template_areas[0]) == 3

        box = widgets.GridspecLayout(1, 2)
        assert box.n_rows == 1
        assert box.n_columns == 2

        with pytest.raises(traitlets.TraitError):
            box = widgets.GridspecLayout(0, 0)

    def test_setitem_index(self):

        box = widgets.GridspecLayout(2, 3)
        button1 = widgets.Button()
        button2 = widgets.Button()
        button3 = widgets.Button()
        button4 = widgets.Button()

        box[0, 0] = button1
        button1_label = button1.layout.grid_area
        assert button1 in box.children
        assert box.layout.grid_template_areas == '''"{} . ."\n". . ."'''.format(button1_label)

        box[-1, -1] = button2
        button2_label = button2.layout.grid_area
        assert button1_label != button2_label
        assert button2 in box.children
        assert box.layout.grid_template_areas == '''"{} . ."\n". . {}"'''.format(button1_label,
                                                                                 button2_label)

        box[1, 0] = button3
        button3_label = button3.layout.grid_area
        assert button1_label != button3_label
        assert button2_label != button3_label
        assert button3 in box.children
        assert box.layout.grid_template_areas == '''"{b1} . ."\n"{b3} . {b2}"'''.format(b1=button1_label,
                                                                                        b2=button2_label,
                                                                                        b3=button3_label)

        #replace widget
        box[1, 0] = button4
        button4_label = button4.layout.grid_area
        assert button1_label != button4_label
        assert button2_label != button4_label
        assert button4 in box.children
        assert button3 not in box.children

        assert box.layout.grid_template_areas == '''"{b1} . ."\n"{b4} . {b2}"'''.format(b1=button1_label,
                                                                                        b2=button2_label,
                                                                                        b4=button4_label)

    def test_setitem_slices(self):

        box = widgets.GridspecLayout(2, 3)
        button1 = widgets.Button()

        box[:2, 0] = button1
        assert len(box.children) == 1
        assert button1 in box.children
        button1_label = button1.layout.grid_area

        assert box.layout.grid_template_areas == '''"{b1} . ."\n"{b1} . ."'''.format(b1=button1_label)

        box = widgets.GridspecLayout(2, 3)
        button1 = widgets.Button()
        button2 = widgets.Button()

        box[:2, 1:] = button1
        assert len(box.children) == 1
        assert button1 in box.children
        button1_label = button1.layout.grid_area

        assert box.layout.grid_template_areas == '''". {b1} {b1}"\n". {b1} {b1}"'''.format(b1=button1_label)

        # replace button
        box[:2, 1:] = button2
        assert len(box.children) == 1
        assert button2 in box.children
        button2_label = button2.layout.grid_area

        assert box.layout.grid_template_areas == '''". {b1} {b1}"\n". {b1} {b1}"'''.format(b1=button2_label)

    def test_getitem_index(self):
        "test retrieving widget"

        box = widgets.GridspecLayout(2, 3)
        button1 = widgets.Button()

        box[0, 0] = button1

        assert box[0, 0] is button1

    def test_getitem_slices(self):
        "test retrieving widgets with slices"

        box = widgets.GridspecLayout(2, 3)
        button1 = widgets.Button()

        box[:2, 0] = button1
        assert box[:2, 0]  is button1

        box = widgets.GridspecLayout(2, 3)
        button1 = widgets.Button()
        button2 = widgets.Button()

        box[0, 0] = button1
        box[1, 0] = button2
        assert box[0, 0]  is button1
        assert box[1, 0]  is button2

        with pytest.raises(TypeError, match="The slice spans"):
            button = box[:2, 0]


class TestLayoutProperties(TestCase):
    """test mixin with layout properties"""

    class DummyTemplate(widgets.GridBox, LayoutProperties):
        location = traitlets.Instance(widgets.Widget, allow_none=True)

    def test_layout_updated_on_trait_change(self):
        "test whether respective layout traits are updated when traits change"

        template = self.DummyTemplate(width="100%")
        assert template.width == '100%'
        assert template.layout.width == '100%'

        template.width = 'auto'
        assert template.width == 'auto'
        assert template.layout.width == 'auto'

    def test_align_items_extra_options(self):

        template = self.DummyTemplate(align_items='top')
        assert template.align_items == 'top'
        assert template.layout.align_items == 'flex-start'

        template.align_items = 'bottom'
        assert template.align_items == 'bottom'
        assert template.layout.align_items == 'flex-end'

    def test_validate_properties(self):

        prop_obj = self.DummyTemplate()
        for prop in LayoutProperties.align_items.values:
            prop_obj.align_items = prop
            assert prop_obj.align_items == prop

        with pytest.raises(traitlets.TraitError):
            prop_obj.align_items = 'any default position'
