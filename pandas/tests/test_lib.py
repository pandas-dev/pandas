import unittest

import numpy as np

from pandas.core.api import TimeStamp

class TestTimestamp(unittest.TestCase):


  def setUp(self):
    self.timestamp = Timestamp(datetime.datetime.utcnow())

  def assert_ns_timedelta(self, modified_timestamp, expected_value):
    value = self.timestamp.value
    modified_value = modified_timestamp.value

    self.assertEquals(modified_value - value, expected_value)

  def test_timedelta_ns_arithmetic(self):
    self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'ns'), -123)

  def test_timedelta_ns_based_arithmetic(self):
    self.assert_ns_timedelta(self.timestamp + np.timedelta64(1234567898, 'ns'), 1234567898)

  def test_timedelta_us_arithmetic(self):
    self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'us'), -123000)

  def test_timedelta_ns_arithmetic(self):
    self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'ms'), -123000000) 
