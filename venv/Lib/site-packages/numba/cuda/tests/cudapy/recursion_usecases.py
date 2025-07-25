"""
Usecases of recursive functions in the CUDA target, many derived from
numba/tests/recursion_usecases.py.

Some functions are compiled at import time, hence a separate module.
"""

from numba import cuda


@cuda.jit("i8(i8)", device=True)
def fib1(n):
    if n < 2:
        return n
    # Note the second call does not use a named argument, unlike the CPU target
    # usecase
    return fib1(n - 1) + fib1(n - 2)


def make_fib2():
    @cuda.jit("i8(i8)", device=True)
    def fib2(n):
        if n < 2:
            return n
        return fib2(n - 1) + fib2(n - 2)

    return fib2


fib2 = make_fib2()


@cuda.jit
def type_change_self(x, y):
    if x > 1 and y > 0:
        return x + type_change_self(x - y, y)
    else:
        return y


# Implicit signature
@cuda.jit(device=True)
def fib3(n):
    if n < 2:
        return n

    return fib3(n - 1) + fib3(n - 2)


# Run-away self recursion
@cuda.jit(device=True)
def runaway_self(x):
    return runaway_self(x)


@cuda.jit(device=True)
def raise_self(x):
    if x == 1:
        raise ValueError("raise_self")
    elif x > 0:
        return raise_self(x - 1)
    else:
        return 1


@cuda.jit(debug=True, opt=False)
def raise_self_kernel(x):
    raise_self(x)


def make_optional_return_case(jit=lambda x: x):
    @jit
    def foo(x):
        if x > 5:
            return x - 1
        else:
            return

    @jit
    def bar(x):
        out = foo(x)
        if out is None:
            return out
        elif out < 8:
            return out
        else:
            return x * bar(out)

    return bar


def make_growing_tuple_case(jit=lambda x: x):
    # From issue #4387
    @jit
    def make_list(n):
        if n <= 0:
            return None

        return (n, make_list(n - 1))
    return make_list
