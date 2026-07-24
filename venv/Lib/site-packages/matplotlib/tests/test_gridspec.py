import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pytest


def test_equal():
    gs = gridspec.GridSpec(2, 1)
    assert gs[0, 0] == gs[0, 0]
    assert gs[:, 0] == gs[:, 0]


def test_update():
    gs = gridspec.GridSpec(2, 1)

    gs.update(left=.1)
    assert gs.left == .1


def test_width_ratios():
    """
    Addresses issue #5835.
    See at https://github.com/matplotlib/matplotlib/issues/5835.
    """
    with pytest.raises(ValueError):
        gridspec.GridSpec(1, 1, width_ratios=[2, 1, 3])


def test_height_ratios():
    """
    Addresses issue #5835.
    See at https://github.com/matplotlib/matplotlib/issues/5835.
    """
    with pytest.raises(ValueError):
        gridspec.GridSpec(1, 1, height_ratios=[2, 1, 3])


def test_SubplotParams():
    s = gridspec.SubplotParams(.1, .1, .9, .9)
    assert s.left == 0.1

    s.reset()
    assert s.left == matplotlib.rcParams['figure.subplot.left']

    with pytest.raises(ValueError, match='left cannot be >= right'):
        s.update(left=s.right + .01)

    with pytest.raises(ValueError, match='bottom cannot be >= top'):
        s.update(bottom=s.top + .01)

    with pytest.raises(ValueError, match='left cannot be >= right'):
        gridspec.SubplotParams(.1, .1, .09, .9)


def test_repr():
    ss = gridspec.GridSpec(3, 3)[2, 1:3]
    assert repr(ss) == "GridSpec(3, 3)[2:3, 1:3]"

    ss = gridspec.GridSpec(2, 2,
                           height_ratios=(3, 1),
                           width_ratios=(1, 3))
    assert repr(ss) == \
        "GridSpec(2, 2, height_ratios=(3, 1), width_ratios=(1, 3))"


def test_subplotspec_args():
    fig, axs = plt.subplots(1, 2)
    # should work:
    gs = gridspec.GridSpecFromSubplotSpec(2, 1,
                                          subplot_spec=axs[0].get_subplotspec())
    assert gs.get_topmost_subplotspec() == axs[0].get_subplotspec()
    with pytest.raises(TypeError, match="subplot_spec must be type SubplotSpec"):
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs[0])
    with pytest.raises(TypeError, match="subplot_spec must be type SubplotSpec"):
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axs)
