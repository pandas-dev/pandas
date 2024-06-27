# Licensed under a 3-clause BSD style license - see LICENSE.rst

class ClassOne:
    pass


class ClassTwo:
    pass


def track_param(n):
    return 42


track_param.params = [ClassOne, ClassTwo]


def mem_param(n, m):
    return [[0] * m] * n


mem_param.params = ([10, 20], [2, 3])
mem_param.param_names = ['number', 'depth']


class ParamSuite:
    params = ['a', 'b', 'c']

    def setup(self, p):
        values = {'a': 1, 'b': 2, 'c': 3}
        self.count = 0
        self.value = values[p]

    def track_value(self, p):
        return self.value + self.count

    def teardown(self, p):
        self.count += 1
        del self.value


class FunctionParamSuite:
    params = [track_param,
              lambda x: x]
    param_names = ['func']

    def time_func(self, func):
        return func('foo')


class TuningTest:
    params = [1, 2]
    counter = [0, 0]
    number = 10
    repeat = 10
    processes = 1
    warmup_time = 0

    def setup(self, n):
        self.number = 1
        self.repeat = n
        self.counter[1] = 0

    def time_it(self, n):
        self.counter[0] += 1
        self.counter[1] += 1

    def teardown(self, n):
        # The time benchmark may call it one additional time
        if not (self.counter[0] <= n + 1 and self.counter[1] == 1):
            raise RuntimeError("Number and repeat didn't have effect: {} {}".format(
                self.counter, n))


def setup_skip(n):
    if n > 2000:
        raise NotImplementedError()


def time_skip(n):
    list(range(n))


time_skip.params = [1000, 2000, 3000]
time_skip.setup = setup_skip
time_skip.sample_time = 0.01


def track_find_test(n):
    import asv_test_repo

    return asv_test_repo.dummy_value[n - 1]


track_find_test.params = [1, 2]


def time_find_test_timeout():
    import time

    import asv_test_repo
    if asv_test_repo.dummy_value[1] < 0:
        time.sleep(100)


time_find_test_timeout.timeout = 1.0
time_find_test_timeout.repeat = 1
time_find_test_timeout.number = 1
time_find_test_timeout.warmup_time = 0


def track_param_selection(a, b):
    return a + b


track_param_selection.param_names = ['a', 'b']
track_param_selection.params = [[1, 2], [3, 5]]


def track_bytes():
    return 1000000


track_bytes.unit = "bytes"


def track_wrong_number_of_args(a, b):
    return 0


track_wrong_number_of_args.params = [[1, 2]]
