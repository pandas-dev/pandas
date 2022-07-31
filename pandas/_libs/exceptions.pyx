import inspect
import os


cpdef int find_stack_level():
    """
    Find the first place in the stack that is not inside pandas
    (tests notwithstanding).
    """
    cdef:
        int n
        str pkg_dir, test_dir, fname

    import pandas as pd


    pkg_dir = os.path.dirname(pd.__file__)
    test_dir = os.path.join(pkg_dir, "tests")

    # https://stackoverflow.com/questions/17407119/python-inspect-stack-is-slow
    frame = inspect.currentframe()
    n = 1
    while frame:
        fname = inspect.getfile(frame)
        if fname.startswith(pkg_dir) and not fname.startswith(test_dir):
            frame = frame.f_back
            n += 1
        else:
            break
    return n
