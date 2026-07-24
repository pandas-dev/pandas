"""Check that functions can handle scalar input"""

from numpy.testing import assert_array_almost_equal
import bottleneck as bn
import pytest


@pytest.mark.parametrize(
    "func",
    bn.get_functions("reduce") + bn.get_functions("nonreduce_axis"),  # noqa: W504
    ids=lambda x: x.__name__,
)
def test_scalar_input(func, args=tuple()):
    """Test that bn.xxx gives the same output as bn.slow.xxx for scalar input."""
    if func.__name__ in ("partition", "argpartition", "push"):
        return
    func0 = eval("bn.slow.%s" % func.__name__)
    msg = "\nfunc %s | input %s\n"
    a = -9
    argsi = [a] + list(args)
    actual = func(*argsi)
    desired = func0(*argsi)
    err_msg = msg % (func.__name__, a)
    assert_array_almost_equal(actual, desired, err_msg=err_msg)
