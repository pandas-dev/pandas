# Python 3 syntax only use cases, used in test_extending.py

# arg name is different, and there's no arg name to match before *args


def impl4(z, *args, kw=None):
    if z > 10:
        return 1
    else:
        return -1

# arg name is different but at a detectable location, with *args


def impl5(z, b, *args, kw=None):
    if z > 10:
        return 1
    else:
        return -1


def var_positional_impl(a, *star_args_token, kw=None, kw1=12):
    def impl(a, b, f, kw=None, kw1=12):
        if a > 10:
            return 1
        else:
            return -1
    return impl
