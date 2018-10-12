from pandas.core.common import *


def foo():
    return (1 + 2
            + 3
            + 4)



def bar():
    """
    DEPRECATED: linting should detect this and complain because it's not
        a sphinx directive
    """
    foobar = 1+2+3+4
    return (1 + 2 +
            3 +
            4)
