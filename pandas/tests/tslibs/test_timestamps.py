import numpy as np
import pytest

from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds

from pandas import Timestamp


def test():
    Timestamp(year=2020)