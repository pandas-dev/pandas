# -*- coding: utf-8 -*-
from datetime import datetime

import pytest

from pandas import Timestamp
from pandas import compat

from pandas.tseries.offsets import (BMonthBegin, BMonthEnd,
                                    MonthBegin, MonthEnd)

from .test_offsets import Base
from .common import assert_offset_equal, assert_onOffset


class TestMonthBegin(Base):
    _offset = MonthBegin

    offset_cases = []
    # NOTE: I'm not entirely happy with the logic here for Begin -ss
    # see thread 'offset conventions' on the ML
    offset_cases.append((MonthBegin(), {
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2008, 2, 1): datetime(2008, 3, 1),
        datetime(2006, 12, 31): datetime(2007, 1, 1),
        datetime(2006, 12, 1): datetime(2007, 1, 1),
        datetime(2007, 1, 31): datetime(2007, 2, 1)}))

    offset_cases.append((MonthBegin(0), {
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2006, 12, 3): datetime(2007, 1, 1),
        datetime(2007, 1, 31): datetime(2007, 2, 1)}))

    offset_cases.append((MonthBegin(2), {
        datetime(2008, 2, 29): datetime(2008, 4, 1),
        datetime(2008, 1, 31): datetime(2008, 3, 1),
        datetime(2006, 12, 31): datetime(2007, 2, 1),
        datetime(2007, 12, 28): datetime(2008, 2, 1),
        datetime(2007, 1, 1): datetime(2007, 3, 1),
        datetime(2006, 11, 1): datetime(2007, 1, 1)}))

    offset_cases.append((MonthBegin(-1), {
        datetime(2007, 1, 1): datetime(2006, 12, 1),
        datetime(2008, 5, 31): datetime(2008, 5, 1),
        datetime(2008, 12, 31): datetime(2008, 12, 1),
        datetime(2006, 12, 29): datetime(2006, 12, 1),
        datetime(2006, 1, 2): datetime(2006, 1, 1)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)


class TestMonthEnd(Base):
    _offset = MonthEnd

    offset_cases = []
    offset_cases.append((MonthEnd(), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 2, 29),
        datetime(2006, 12, 29): datetime(2006, 12, 31),
        datetime(2006, 12, 31): datetime(2007, 1, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 31),
        datetime(2006, 12, 1): datetime(2006, 12, 31)}))

    offset_cases.append((MonthEnd(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 1, 31),
        datetime(2006, 12, 29): datetime(2006, 12, 31),
        datetime(2006, 12, 31): datetime(2006, 12, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 31)}))

    offset_cases.append((MonthEnd(2), {
        datetime(2008, 1, 1): datetime(2008, 2, 29),
        datetime(2008, 1, 31): datetime(2008, 3, 31),
        datetime(2006, 12, 29): datetime(2007, 1, 31),
        datetime(2006, 12, 31): datetime(2007, 2, 28),
        datetime(2007, 1, 1): datetime(2007, 2, 28),
        datetime(2006, 11, 1): datetime(2006, 12, 31)}))

    offset_cases.append((MonthEnd(-1), {
        datetime(2007, 1, 1): datetime(2006, 12, 31),
        datetime(2008, 6, 30): datetime(2008, 5, 31),
        datetime(2008, 12, 31): datetime(2008, 11, 30),
        datetime(2006, 12, 29): datetime(2006, 11, 30),
        datetime(2006, 12, 30): datetime(2006, 11, 30),
        datetime(2007, 1, 1): datetime(2006, 12, 31)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_day_of_month(self):
        dt = datetime(2007, 1, 1)
        offset = MonthEnd()

        result = dt + offset
        assert result == Timestamp(2007, 1, 31)

        result = result + offset
        assert result == Timestamp(2007, 2, 28)

    def test_normalize(self):
        dt = datetime(2007, 1, 1, 3)

        result = dt + MonthEnd(normalize=True)
        expected = dt.replace(hour=0) + MonthEnd()
        assert result == expected

    on_offset_cases = [(MonthEnd(), datetime(2007, 12, 31), True),
                       (MonthEnd(), datetime(2008, 1, 1), False)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, dt, expected = case
        assert_onOffset(offset, dt, expected)


class TestBMonthBegin(Base):
    _offset = BMonthBegin

    offset_cases = []
    offset_cases.append((BMonthBegin(), {
        datetime(2008, 1, 1): datetime(2008, 2, 1),
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2006, 12, 29): datetime(2007, 1, 1),
        datetime(2006, 12, 31): datetime(2007, 1, 1),
        datetime(2006, 9, 1): datetime(2006, 10, 2),
        datetime(2007, 1, 1): datetime(2007, 2, 1),
        datetime(2006, 12, 1): datetime(2007, 1, 1)}))

    offset_cases.append((BMonthBegin(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 1),
        datetime(2006, 10, 2): datetime(2006, 10, 2),
        datetime(2008, 1, 31): datetime(2008, 2, 1),
        datetime(2006, 12, 29): datetime(2007, 1, 1),
        datetime(2006, 12, 31): datetime(2007, 1, 1),
        datetime(2006, 9, 15): datetime(2006, 10, 2)}))

    offset_cases.append((BMonthBegin(2), {
        datetime(2008, 1, 1): datetime(2008, 3, 3),
        datetime(2008, 1, 15): datetime(2008, 3, 3),
        datetime(2006, 12, 29): datetime(2007, 2, 1),
        datetime(2006, 12, 31): datetime(2007, 2, 1),
        datetime(2007, 1, 1): datetime(2007, 3, 1),
        datetime(2006, 11, 1): datetime(2007, 1, 1)}))

    offset_cases.append((BMonthBegin(-1), {
        datetime(2007, 1, 1): datetime(2006, 12, 1),
        datetime(2008, 6, 30): datetime(2008, 6, 2),
        datetime(2008, 6, 1): datetime(2008, 5, 1),
        datetime(2008, 3, 10): datetime(2008, 3, 3),
        datetime(2008, 12, 31): datetime(2008, 12, 1),
        datetime(2006, 12, 29): datetime(2006, 12, 1),
        datetime(2006, 12, 30): datetime(2006, 12, 1),
        datetime(2007, 1, 1): datetime(2006, 12, 1)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    on_offset_cases = [(BMonthBegin(), datetime(2007, 12, 31), False),
                       (BMonthBegin(), datetime(2008, 1, 1), True),
                       (BMonthBegin(), datetime(2001, 4, 2), True),
                       (BMonthBegin(), datetime(2008, 3, 3), True)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, dt, expected = case
        assert_onOffset(offset, dt, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = BMonthBegin()
        offset2 = BMonthBegin()
        assert not offset1 != offset2


class TestBMonthEnd(Base):
    _offset = BMonthEnd

    offset_cases = []
    offset_cases.append((BMonthEnd(), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 2, 29),
        datetime(2006, 12, 29): datetime(2007, 1, 31),
        datetime(2006, 12, 31): datetime(2007, 1, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 31),
        datetime(2006, 12, 1): datetime(2006, 12, 29)}))

    offset_cases.append((BMonthEnd(0), {
        datetime(2008, 1, 1): datetime(2008, 1, 31),
        datetime(2008, 1, 31): datetime(2008, 1, 31),
        datetime(2006, 12, 29): datetime(2006, 12, 29),
        datetime(2006, 12, 31): datetime(2007, 1, 31),
        datetime(2007, 1, 1): datetime(2007, 1, 31)}))

    offset_cases.append((BMonthEnd(2), {
        datetime(2008, 1, 1): datetime(2008, 2, 29),
        datetime(2008, 1, 31): datetime(2008, 3, 31),
        datetime(2006, 12, 29): datetime(2007, 2, 28),
        datetime(2006, 12, 31): datetime(2007, 2, 28),
        datetime(2007, 1, 1): datetime(2007, 2, 28),
        datetime(2006, 11, 1): datetime(2006, 12, 29)}))

    offset_cases.append((BMonthEnd(-1), {
        datetime(2007, 1, 1): datetime(2006, 12, 29),
        datetime(2008, 6, 30): datetime(2008, 5, 30),
        datetime(2008, 12, 31): datetime(2008, 11, 28),
        datetime(2006, 12, 29): datetime(2006, 11, 30),
        datetime(2006, 12, 30): datetime(2006, 12, 29),
        datetime(2007, 1, 1): datetime(2006, 12, 29)}))

    @pytest.mark.parametrize('case', offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in compat.iteritems(cases):
            assert_offset_equal(offset, base, expected)

    def test_normalize(self):
        dt = datetime(2007, 1, 1, 3)

        result = dt + BMonthEnd(normalize=True)
        expected = dt.replace(hour=0) + BMonthEnd()
        assert result == expected

    on_offset_cases = [(BMonthEnd(), datetime(2007, 12, 31), True),
                       (BMonthEnd(), datetime(2008, 1, 1), False)]

    @pytest.mark.parametrize('case', on_offset_cases)
    def test_onOffset(self, case):
        offset, dt, expected = case
        assert_onOffset(offset, dt, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = BMonthEnd()
        offset2 = BMonthEnd()
        assert not offset1 != offset2
