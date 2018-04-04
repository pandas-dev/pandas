import pytest

from .base import BaseExtensionTests


class BaseSetitemTests(BaseExtensionTests):
    """Tests for ExtensionArray.__setitem__"""

    def test_set_scalar(self, data):
        expected = data.take([1, 1])
        subset = data[:2].copy()

        subset[0] = data[1]
        self.assert_extension_array_equal(subset, expected)

    def test_set_mask_scalar(self, data):
        expected = data.take([1, 1, 2, 1])
        subset = data[:4].copy()

        subset[[True, True, False, True]] = data[1]
        self.assert_extension_array_equal(subset, expected)

    @pytest.mark.parametrize('key', [
        [False, True, True, True],
        [1, 2, 3],
    ], ids=['mask', 'fancy'])
    def test_set_array(self, key, data):
        expected = data.take([0, 2, 2, 1])
        value = data.take([2, 2, 1])
        subset = data[:4].copy()

        subset[key] = value
        self.assert_extension_array_equal(subset, expected)

    def test_bad_mask_bad_length_raise(self, data):
        value = data[0]
        with pytest.raises(IndexError):
            data[[True, False]] = value
