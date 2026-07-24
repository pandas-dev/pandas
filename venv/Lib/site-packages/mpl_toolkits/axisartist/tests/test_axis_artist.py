import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
from matplotlib.testing.decorators import image_comparison
from matplotlib.transforms import Affine2D

from mpl_toolkits.axisartist import (AxisArtistHelperRectlinear, GridHelperCurveLinear,
                                     HostAxes)
from mpl_toolkits.axisartist.axis_artist import (AxisArtist, AxisLabel,
                                                 LabelBase, Ticks, TickLabels)


@image_comparison(['axis_artist_ticks.png'], style='default')
def test_ticks():
    fig, ax = plt.subplots()

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    locs_angles = [((i / 10, 0.0), i * 30) for i in range(-1, 12)]

    ticks_in = Ticks(ticksize=10, axis=ax.xaxis)
    ticks_in.set_locs_angles(locs_angles)
    ax.add_artist(ticks_in)

    ticks_out = Ticks(ticksize=10, tick_direction="out", color='C3', axis=ax.xaxis)
    ticks_out.set_locs_angles(locs_angles)
    ax.add_artist(ticks_out)


@image_comparison(['axis_artist_labelbase.png'], style='default')
def test_labelbase():
    fig, ax = plt.subplots()

    ax.plot([0.5], [0.5], "o")

    label = LabelBase(0.5, 0.5, "Test")
    label._ref_angle = -90
    label._offset_radius = 50
    label.set_rotation(-90)
    label.set(ha="center", va="top")
    ax.add_artist(label)


@image_comparison(['axis_artist_ticklabels.png'], style='default')
def test_ticklabels():
    fig, ax = plt.subplots()

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    ax.plot([0.2, 0.4], [0.5, 0.5], "o")

    ticks = Ticks(ticksize=10, axis=ax.xaxis)
    ax.add_artist(ticks)
    locs_angles_labels = [((0.2, 0.5), -90, "0.2"),
                          ((0.4, 0.5), -120, "0.4")]
    tick_locs_angles = [(xy, a + 180) for xy, a, l in locs_angles_labels]
    ticks.set_locs_angles(tick_locs_angles)

    ticklabels = TickLabels(axis_direction="left")
    ticklabels._locs_angles_labels = locs_angles_labels
    ticklabels.set_pad(10)
    ax.add_artist(ticklabels)

    ax.plot([0.5], [0.5], "s")
    axislabel = AxisLabel(0.5, 0.5, "Test")
    axislabel._offset_radius = 20
    axislabel._ref_angle = 0
    axislabel.set_axis_direction("bottom")
    ax.add_artist(axislabel)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


@image_comparison(['axis_artist.png'], style='default')
def test_axis_artist():
    fig, ax = plt.subplots()

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    for loc in ('left', 'right', 'bottom'):
        helper = AxisArtistHelperRectlinear.Fixed(ax, loc=loc)
        axisline = AxisArtist(ax, helper, offset=None, axis_direction=loc)
        axisline.major_ticks.set_tick_direction({
            "left": "in", "right": "out", "bottom": "inout",
        }[loc])
        ax.add_artist(axisline)

    # Settings for bottom AxisArtist.
    axisline.set_label("TTT")
    axisline.label.set_pad(5)

    ax.set_ylabel("Test")


@mpl.style.context('default')
def test_axisartist_tightbbox():
    fig = plt.figure()
    tr = Affine2D().scale(np.pi / 180., 1.) + PolarAxes.PolarTransform()
    grid_helper = GridHelperCurveLinear(tr)
    ax = fig.add_subplot(axes_class=HostAxes, grid_helper=grid_helper)
    ax.axis["lon"] = ax.new_floating_axis(1, 9)

    ax.set_xlim(-5, 12)
    ax.set_ylim(-5, 10)

    ax.axis['lon'].major_ticklabels.set_visible(False)

    # Since the labels are invisible and the lines are clipped to the axes,
    # the axis's tight bbox should be contained in the axes box.
    renderer = fig._get_renderer()
    tight_points = ax.axis['lon'].get_tightbbox(renderer).get_points()
    for point in tight_points:
        assert ax.bbox.contains(*point)
