# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import inspect
from unittest import TestCase

from traitlets import TraitError

from ipywidgets import Dropdown, SelectionSlider, Select


class TestDropdown(TestCase):

    def test_construction(self):
        Dropdown()

    def test_dict_mapping_options(self):
        d = Dropdown(options={'One': 1, 'Two': 2, 'Three': 3})
        assert d.get_state('_options_labels') == {'_options_labels': ('One', 'Two', 'Three')}

    def test_setting_options_from_list(self):
        d = Dropdown()
        assert d.options == ()
        d.options = ['One', 'Two', 'Three']
        assert d.get_state('_options_labels') == {'_options_labels': ('One', 'Two', 'Three')}

    def test_setting_options_from_list_tuples(self):
        d = Dropdown()
        assert d.options == ()
        d.options = [('One', 1), ('Two', 2), ('Three', 3)]
        assert d.get_state('_options_labels') == {'_options_labels': ('One', 'Two', 'Three')}
        d.value = 2
        assert d.get_state('index') == {'index': 1}

    def test_setting_options_from_dict(self):
        d = Dropdown()
        assert d.options == ()
        d.options = {'One': 1, 'Two': 2, 'Three': 3}
        assert d.get_state('_options_labels') == {'_options_labels': ('One', 'Two', 'Three')}




class TestSelectionSlider(TestCase):

    def test_construction(self):
        SelectionSlider(options=['a', 'b', 'c'])

    def test_index_trigger(self):
        slider = SelectionSlider(options=['a', 'b', 'c'])
        observations = []
        def f(change):
            observations.append(change.new)
        slider.observe(f, 'index')
        assert slider.index == 0
        slider.options = [4, 5, 6]
        assert slider.index == 0
        assert slider.value == 4
        assert slider.label == '4'
        assert observations == [0]

class TestSelection(TestCase):

    def test_construction(self):
        select = Select(options=['a', 'b', 'c'])

    def test_index_trigger(self):
        select = Select(options=[1, 2, 3])
        observations = []
        def f(change):
            observations.append(change.new)
        select.observe(f, 'index')
        assert select.index == 0
        select.options = [4, 5, 6]
        assert select.index == 0
        assert select.value == 4
        assert select.label == '4'
        assert observations == [0]

    def test_duplicate(self):
        select = Select(options=['first', 1, 'dup', 'dup'])
        observations = []
        def f(change):
            observations.append(change.new)
        select.observe(f, 'index')
        select.index = 3
        assert select.index == 3
        assert select.value == 'dup'
        assert select.label == 'dup'
        assert observations == [3]
        select.index = 2
        assert select.index == 2
        assert select.value == 'dup'
        assert select.label == 'dup'
        assert observations == [3, 2]
        select.index = 0
        assert select.index == 0
        assert select.value == 'first'
        assert select.label == 'first'
        assert observations == [3, 2, 0]

        # picks the first matching value
        select.value = 'dup'
        assert select.index == 2
        assert select.value == 'dup'
        assert select.label == 'dup'
        assert observations == [3, 2, 0, 2]
