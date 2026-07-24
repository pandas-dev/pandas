# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from unittest import TestCase

from traitlets import TraitError

import ipywidgets as widgets


class TestBox(TestCase):

    def test_construction(self):
        box = widgets.Box()
        assert box.get_state()['children'] == []

    def test_construction_with_children(self):
        html = widgets.HTML('some html')
        slider = widgets.IntSlider()
        box = widgets.Box([html, slider])
        children_state = box.get_state()['children']
        assert children_state == [
            widgets.widget._widget_to_json(html, None),
            widgets.widget._widget_to_json(slider, None),
        ]

    def test_construction_style(self):
        box = widgets.Box(box_style='warning')
        assert box.get_state()['box_style'] == 'warning'

    def test_construction_invalid_style(self):
        with self.assertRaises(TraitError):
            widgets.Box(box_style='invalid')
