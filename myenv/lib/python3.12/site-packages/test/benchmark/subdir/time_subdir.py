x = None


def time_foo():
    if x != 42:
        raise RuntimeError()
    for y in range(1000):
        pass


def setup_foo():
    global x
    x = 42


time_foo.setup = setup_foo
time_foo.sample_time = 0.1
