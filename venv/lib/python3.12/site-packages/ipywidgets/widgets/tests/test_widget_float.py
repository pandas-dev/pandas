# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from unittest import TestCase

from traitlets import TraitError

from ipywidgets import FloatSlider


class TestFloatSlider(TestCase):

    def test_construction(self):
        FloatSlider()

    def test_construction_readout_format(self):
        slider = FloatSlider(readout_format='$.1f')
        assert slider.get_state()['readout_format'] == '$.1f'

    def test_construction_invalid_readout_format(self):
        with self.assertRaises(TraitError):
            FloatSlider(readout_format='broken')
