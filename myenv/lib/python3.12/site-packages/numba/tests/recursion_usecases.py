"""
Usecases of recursive functions.

Some functions are compiled at import time, hence a separate module.
"""

from numba import jit


@jit("i8(i8)", nopython=True)
def fib1(n):
    if n < 2:
        return n
    # Note the second call uses a named argument
    return fib1(n - 1) + fib1(n=n - 2)


def make_fib2():
    @jit("i8(i8)", nopython=True)
    def fib2(n):
        if n < 2:
            return n
        return fib2(n - 1) + fib2(n=n - 2)

    return fib2

fib2 = make_fib2()


def make_type_change_self(jit=lambda x: x):
    @jit
    def type_change_self(x, y):
        if x > 1 and y > 0:
            return x + type_change_self(x - y, y)
        else:
            return y
    return type_change_self


# Implicit signature
@jit(nopython=True)
def fib3(n):
    if n < 2:
        return n
    return fib3(n - 1) + fib3(n - 2)


# Run-away self recursion
@jit(nopython=True)
def runaway_self(x):
    return runaway_self(x)


@jit(nopython=True)
def raise_self(x):
    if x == 1:
        raise ValueError("raise_self")
    elif x > 0:
        return raise_self(x - 1)
    else:
        return 1


# Mutual recursion
@jit(nopython=True)
def outer_fac(n):
    if n < 1:
        return 1
    return n * inner_fac(n - 1)


@jit(nopython=True)
def inner_fac(n):
    if n < 1:
        return 1
    return n * outer_fac(n - 1)


# Mutual recursion with different arg names
def make_mutual2(jit=lambda x: x):
    @jit
    def foo(x):
        if x > 0:
            return 2 * bar(z=1, y=x)
        return 1 + x

    @jit
    def bar(y, z):
        return foo(x=y - z)

    return foo, bar


# Mutual runaway recursion

@jit(nopython=True)
def runaway_mutual(x):
    return runaway_mutual_inner(x)


@jit(nopython=True)
def runaway_mutual_inner(x):
    return runaway_mutual(x)


# Mutual type changing recursion

def make_type_change_mutual(jit=lambda x: x):
    @jit
    def foo(x, y):
        if x > 1 and y > 0:
            # call bar first to exercise partial type inference.
            # typeinferer suspended at the call to bar() and haven't determined
            # the potential return type from the else-branch
            return x + bar(x - y, y)
        else:
            return y

    @jit
    def bar(x, y):
        if x > 1 and y > 0:
            return x + foo(x - y, y)
        else:
            return y

    return foo


# Indirect mutual recursion
def make_four_level(jit=lambda x: x):
    @jit
    def first(x):
        # The recursing call must have a path that is non-recursing.
        if x > 0:
            return second(x) * 2
        else:
            return 1

    @jit
    def second(x):
        return third(x) * 3

    @jit
    def third(x):
        return fourth(x) * 4

    @jit
    def fourth(x):
        return first(x / 2 - 1)

    return first


def make_inner_error(jit=lambda x: x):
    @jit
    def outer(x):
        if x > 0:
            return inner(x)

        else:
            return 1

    @jit
    def inner(x):
        if x > 0:
            return outer(x - 1)
        else:
            # this branch is actually never executed
            return error_fun(x)

    @jit
    def error_fun(x):
        # to trigger an untyped attribute error
        return x.ndim

    return outer


def make_raise_mutual(jit=lambda x: x):
    @jit
    def outer(x):
        if x > 0:
            return inner(x)
        else:
            return 1

    @jit
    def inner(x):
        if x == 1:
            raise ValueError('raise_mutual')
        elif x > 0:
            return outer(x - 1)
        else:
            return 1

    return outer


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
