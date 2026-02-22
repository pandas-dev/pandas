# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from unittest import TestCase

from traitlets import TraitError

from ipywidgets.widgets import Accordion, Tab, Stack, HTML

class TestTab(TestCase):

    def setUp(self):
        self.children = [HTML('0'), HTML('1')]
        self.widget = Tab

    def test_selected_index_none(self):
        widget = self.widget(self.children, selected_index=None)
        state = widget.get_state()
        assert state['selected_index'] == 0

    def test_selected_index_default(self):
        widget = self.widget(self.children)
        state = widget.get_state()
        assert state['selected_index'] == 0

    def test_selected_index(self):
        widget = self.widget(self.children, selected_index=1)
        state = widget.get_state()
        assert state['selected_index'] == 1

    def test_selected_index_out_of_bounds(self):
        with self.assertRaises(TraitError):
            self.widget(self.children, selected_index=-1)

    def test_children_position_argument(self):
        self.widget(self.children)

    def test_titles(self):
        widget = self.widget(self.children, selected_index=None)
        assert widget.get_state()['titles'] == ('', '')
        assert widget.titles == ('', '')

        widget.set_title(1, 'Title 1')
        assert widget.get_state()['titles'] == ('', 'Title 1')
        assert widget.titles[1] == 'Title 1'
        assert widget.get_title(1) == 'Title 1'

        # Backwards compatible with 7.x api
        widget.set_title(1, None)
        assert widget.get_state()['titles'] == ('', '')
        assert widget.titles[1] == ''
        assert widget.get_title(1) == ''

        with self.assertRaises(IndexError):
            widget.set_title(2, 'out of bounds')
        with self.assertRaises(IndexError):
            widget.get_title(2)

        widget.children = tuple(widget.children[:1])
        assert len(widget.children) == 1
        assert widget.titles == ('',)

class TestAccordion(TestCase):

    def setUp(self):
        self.children = [HTML('0'), HTML('1')]
        self.widget = Accordion

    def test_selected_index_none(self):
        widget = self.widget(self.children, selected_index=None)
        state = widget.get_state()
        assert state['selected_index'] is None

    def test_selected_index_default(self):
        widget = self.widget(self.children)
        state = widget.get_state()
        assert state['selected_index'] is None

    def test_selected_index(self):
        widget = self.widget(self.children, selected_index=1)
        state = widget.get_state()
        assert state['selected_index'] == 1

    def test_selected_index_out_of_bounds(self):
        with self.assertRaises(TraitError):
            self.widget(self.children, selected_index=-1)

    def test_children_position_argument(self):
        self.widget(self.children)

    def test_titles(self):
        widget = self.widget(self.children, selected_index=None)
        assert widget.get_state()['titles'] == ('', '')
        assert widget.titles == ('', '')

        widget.set_title(1, 'Title 1')
        assert widget.get_state()['titles'] == ('', 'Title 1')
        assert widget.titles[1] == 'Title 1'
        assert widget.get_title(1) == 'Title 1'

        # Backwards compatible with 7.x api
        widget.set_title(1, None)
        assert widget.get_state()['titles'] == ('', '')
        assert widget.titles[1] == ''
        assert widget.get_title(1) == ''

        with self.assertRaises(IndexError):
            widget.set_title(2, 'out of bounds')
        with self.assertRaises(IndexError):
            widget.get_title(2)

        widget.children = tuple(widget.children[:1])
        assert len(widget.children) == 1
        assert widget.titles == ('',)

class TestStack(TestCase):

    def setUp(self):
        self.children = [HTML('0'), HTML('1')]
        self.widget = Stack

    def test_selected_index_none(self):
        widget = self.widget(self.children, selected_index=None)
        state = widget.get_state()
        assert state['selected_index'] is None

    def test_selected_index_default(self):
        widget = self.widget(self.children)
        state = widget.get_state()
        assert state['selected_index'] is None

    def test_selected_index(self):
        widget = self.widget(self.children, selected_index=1)
        state = widget.get_state()
        assert state['selected_index'] == 1

    def test_selected_index_out_of_bounds(self):
        with self.assertRaises(TraitError):
            self.widget(self.children, selected_index=-1)

    def test_children_position_argument(self):
        self.widget(self.children)

    def test_titles(self):
        widget = self.widget(self.children, selected_index=None)
        assert widget.get_state()['titles'] == ('', '')
        assert widget.titles == ('', '')

        widget.set_title(1, 'Title 1')
        assert widget.get_state()['titles'] == ('', 'Title 1')
        assert widget.titles[1] == 'Title 1'
        assert widget.get_title(1) == 'Title 1'

        # Backwards compatible with 7.x api
        widget.set_title(1, None)
        assert widget.get_state()['titles'] == ('', '')
        assert widget.titles[1] == ''
        assert widget.get_title(1) == ''

        with self.assertRaises(IndexError):
            widget.set_title(2, 'out of bounds')
        with self.assertRaises(IndexError):
            widget.get_title(2)

        widget.children = tuple(widget.children[:1])
        assert len(widget.children) == 1
        assert widget.titles == ('',)

