"""
Tests for CategoricalIndex.__repr__ and related methods.
"""

import pytest

from pandas._config import using_string_dtype
import pandas._config.config as cf

from pandas import CategoricalIndex


class TestCategoricalIndexReprStringCategories:
    def test_string_categorical_index_repr(self, using_infer_string):
        # short
        idx = CategoricalIndex(["a", "bb", "ccc"])
        expected = """CategoricalIndex(['a', 'bb', 'ccc'], categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['a', 'bb', 'ccc']",
                "categories=[a, bb, ccc]",
            )
        assert repr(idx) == expected

    @pytest.mark.xfail(using_string_dtype(), reason="Different padding on multi-line")
    def test_categorical_index_repr_multiline(self, using_infer_string):
        # multiple lines
        idx = CategoricalIndex(["a", "bb", "ccc"] * 10)
        expected = """CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb',
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category')"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['a', 'bb', 'ccc']",
                "categories=[a, bb, ccc]",
            )
        assert repr(idx) == expected

    @pytest.mark.xfail(using_string_dtype(), reason="Different padding on multi-line")
    def test_categorical_index_repr_truncated(self, using_infer_string):
        # truncated
        idx = CategoricalIndex(["a", "bb", "ccc"] * 100)
        expected = """CategoricalIndex(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
                  ...
                  'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
                 categories=['a', 'bb', 'ccc'], ordered=False, dtype='category', length=300)"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['a', 'bb', 'ccc']",
                "categories=[a, bb, ccc]",
            )
        assert repr(idx) == expected

    def test_categorical_index_repr_many_categories(self, using_infer_string):
        # larger categories
        idx = CategoricalIndex(list("abcdefghijklmmo"))
        expected = """CategoricalIndex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                  'm', 'm', 'o'],
                 categories=['a', 'b', 'c', 'd', ..., 'k', 'l', 'm', 'o'], ordered=False, dtype='category')"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['a', 'b', 'c', 'd', ..., 'k', 'l', 'm', 'o']",
                "categories=[a, b, c, d, ..., k, l, m, o]",
            )
        assert repr(idx) == expected

    def test_categorical_index_repr_unicode(self, using_infer_string):
        # short
        idx = CategoricalIndex(["あ", "いい", "ううう"])
        expected = """CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['あ', 'いい', 'ううう']",
                "categories=[あ, いい, ううう]",
            )
        assert repr(idx) == expected

    @pytest.mark.xfail(using_string_dtype(), reason="Different padding on multi-line")
    def test_categorical_index_repr_unicode_multiline(self, using_infer_string):
        # multiple lines
        idx = CategoricalIndex(["あ", "いい", "ううう"] * 10)
        expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['あ', 'いい', 'ううう']",
                "categories=[あ, いい, ううう]",
            )
        assert repr(idx) == expected

    @pytest.mark.xfail(using_string_dtype(), reason="Different padding on multi-line")
    def test_categorical_index_repr_unicode_truncated(self, using_infer_string):
        # truncated
        idx = CategoricalIndex(["あ", "いい", "ううう"] * 100)
        expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['あ', 'いい', 'ううう']",
                "categories=[あ, いい, ううう]",
            )
        assert repr(idx) == expected

    def test_categorical_index_repr_unicode_many_categories(self, using_infer_string):
        # larger categories
        idx = CategoricalIndex(list("あいうえおかきくけこさしすせそ"))
        expected = """CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ', 'さ', 'し',
                  'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', ..., 'し', 'す', 'せ', 'そ'], ordered=False, dtype='category')"""  # noqa: E501
        if using_infer_string:
            expected = expected.replace(
                "categories=['あ', 'い', 'う', 'え', ..., 'し', 'す', 'せ', 'そ']",
                "categories=[あ, い, う, え, ..., し, す, せ, そ]",
            )
        assert repr(idx) == expected

    def test_categorical_index_repr_east_asian_width(self, using_infer_string):
        with cf.option_context("display.unicode.east_asian_width", True):
            # short
            idx = CategoricalIndex(["あ", "いい", "ううう"])
            expected = """CategoricalIndex(['あ', 'いい', 'ううう'], categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa: E501
            if using_infer_string:
                expected = expected.replace(
                    "categories=['あ', 'いい', 'ううう']",
                    "categories=[あ, いい, ううう]",
                )
            assert repr(idx) == expected

    @pytest.mark.xfail(using_string_dtype(), reason="Different padding on multi-line")
    def test_categorical_index_repr_east_asian_width_multiline(
        self, using_infer_string
    ):
        with cf.option_context("display.unicode.east_asian_width", True):
            # multiple lines
            idx = CategoricalIndex(["あ", "いい", "ううう"] * 10)
            expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category')"""  # noqa: E501

            if using_infer_string:
                expected = expected.replace(
                    "categories=['あ', 'いい', 'ううう']",
                    "categories=[あ, いい, ううう]",
                )
            assert repr(idx) == expected

    @pytest.mark.xfail(using_string_dtype(), reason="Different padding on multi-line")
    def test_categorical_index_repr_east_asian_width_truncated(
        self, using_infer_string
    ):
        with cf.option_context("display.unicode.east_asian_width", True):
            # truncated
            idx = CategoricalIndex(["あ", "いい", "ううう"] * 100)
            expected = """CategoricalIndex(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
                  'ううう', 'あ',
                  ...
                  'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
                  'あ', 'いい', 'ううう'],
                 categories=['あ', 'いい', 'ううう'], ordered=False, dtype='category', length=300)"""  # noqa: E501

            if using_infer_string:
                expected = expected.replace(
                    "categories=['あ', 'いい', 'ううう']",
                    "categories=[あ, いい, ううう]",
                )
            assert repr(idx) == expected

    def test_categorical_index_repr_east_asian_width_many_categories(
        self, using_infer_string
    ):
        with cf.option_context("display.unicode.east_asian_width", True):
            idx = CategoricalIndex(list("あいうえおかきくけこさしすせそ"))
            expected = """CategoricalIndex(['あ', 'い', 'う', 'え', 'お', 'か', 'き', 'く', 'け', 'こ',
                  'さ', 'し', 'す', 'せ', 'そ'],
                 categories=['あ', 'い', 'う', 'え', ..., 'し', 'す', 'せ', 'そ'], ordered=False, dtype='category')"""  # noqa: E501

            if using_infer_string:
                expected = expected.replace(
                    "categories=['あ', 'い', 'う', 'え', ..., 'し', 'す', 'せ', 'そ']",
                    "categories=[あ, い, う, え, ..., し, す, せ, そ]",
                )
            assert repr(idx) == expected
