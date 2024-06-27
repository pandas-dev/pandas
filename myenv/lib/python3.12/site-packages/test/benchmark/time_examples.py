# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

if sys.version_info[0] == 3:
    xrange = range

import warnings


class TimeSuite:
    sample_time = 0.1

    def setup(self):
        self.n = 100

    def time_example_benchmark_1(self):
        s = ''
        for i in xrange(self.n):
            s = s + 'x'

    def time_example_benchmark_2(self):
        s = []
        for i in xrange(self.n):
            s.append('x')
        ''.join(s)


class TimeSuiteSub(TimeSuite):
    pass


def time_with_warnings():
    print('hi')
    warnings.warn('before')
    1 / 0
    warnings.warn('after')


time_with_warnings.sample_time = 0.1


def time_with_timeout():
    while True:
        pass


time_with_timeout.timeout = 0.1


class TimeWithRepeat:
    # Check that setup is re-run on each repeat
    called = None
    number = 1
    repeat = 10
    count = 0
    warmup_time = 0

    def setup(self):
        assert self.called is None
        self.called = False

    def teardown(self):
        assert self.called is True
        self.called = None
        print("<%d>" % (self.count,))

    def time_it(self):
        assert self.called is False
        self.called = True
        self.count += 1


class TimeWithRepeatCalibrate:
    # Check that setup is re-run on each repeat, apart from
    # autodetection of suitable `number`
    repeat = 1
    number = 0
    sample_time = 0.1

    def setup(self):
        print("setup")

    def time_it(self):
        pass


class TimeWithBadTimer:
    # Check that calibration of number is robust against bad timers
    repeat = 1
    number = 0
    sample_time = 0.1
    timeout = 5

    def timer(self):
        return 0.0

    def time_it(self):
        pass


def time_auto_repeat():
    pass


time_auto_repeat.number = 1
time_auto_repeat.processes = 1
time_auto_repeat.repeat = (2, 4, 10.0)
