from datetime import datetime
import gc
import inspect
import io
import warnings

import numpy as np
from numpy.testing import assert_almost_equal
from packaging.version import parse as parse_version
import pyparsing
import pytest

import matplotlib as mpl
from matplotlib.backend_bases import MouseEvent
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties, fontManager, get_font
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as mtransforms
from matplotlib.testing.decorators import check_figures_equal, image_comparison
from matplotlib.testing._markers import needs_usetex
from matplotlib.text import Text, Annotation, OffsetFrom

pyparsing_version = parse_version(pyparsing.__version__)


@image_comparison(['font_styles'], style='mpl20')
def test_font_styles():

    def find_matplotlib_font(**kw):
        prop = FontProperties(**kw)
        path = findfont(prop, directory=mpl.get_data_path())
        return FontProperties(fname=path)

    from matplotlib.font_manager import FontProperties, findfont
    warnings.filterwarnings(
        'ignore',
        r"findfont: Font family \[u?'Foo'\] not found. Falling back to .",
        UserWarning,
        module='matplotlib.font_manager')

    fig, ax = plt.subplots()

    normal_font = find_matplotlib_font(
        family="sans-serif",
        style="normal",
        variant="normal",
        size=14)
    a = ax.annotate(
        "Normal Font",
        (0.1, 0.1),
        xycoords='axes fraction',
        fontproperties=normal_font)
    assert a.get_fontname() == 'DejaVu Sans'
    assert a.get_fontstyle() == 'normal'
    assert a.get_fontvariant() == 'normal'
    assert a.get_weight() == 'normal'
    assert a.get_stretch() == 'normal'

    bold_font = find_matplotlib_font(
        family="Foo",
        style="normal",
        variant="normal",
        weight="bold",
        stretch=500,
        size=14)
    ax.annotate(
        "Bold Font",
        (0.1, 0.2),
        xycoords='axes fraction',
        fontproperties=bold_font)

    bold_italic_font = find_matplotlib_font(
        family="sans serif",
        style="italic",
        variant="normal",
        weight=750,
        stretch=500,
        size=14)
    ax.annotate(
        "Bold Italic Font",
        (0.1, 0.3),
        xycoords='axes fraction',
        fontproperties=bold_italic_font)

    light_font = find_matplotlib_font(
        family="sans-serif",
        style="normal",
        variant="normal",
        weight=200,
        stretch=500,
        size=14)
    ax.annotate(
        "Light Font",
        (0.1, 0.4),
        xycoords='axes fraction',
        fontproperties=light_font)

    condensed_font = find_matplotlib_font(
        family="sans-serif",
        style="normal",
        variant="normal",
        weight=500,
        stretch=100,
        size=14)
    ax.annotate(
        "Condensed Font",
        (0.1, 0.5),
        xycoords='axes fraction',
        fontproperties=condensed_font)

    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(['multiline'], style='mpl20')
def test_multiline():
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.set_title("multiline\ntext alignment")

    plt.text(
        0.2, 0.5, "TpTpTp\n$M$\nTpTpTp", size=20, ha="center", va="top")

    plt.text(
        0.5, 0.5, "TpTpTp\n$M^{M^{M^{M}}}$\nTpTpTp", size=20,
        ha="center", va="top")

    plt.text(
        0.8, 0.5, "TpTpTp\n$M_{q_{q_{q}}}$\nTpTpTp", size=20,
        ha="center", va="top")

    plt.xlim(0, 1)
    plt.ylim(0, 0.8)

    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(['multiline2'], style='mpl20')
def test_multiline2():
    fig, ax = plt.subplots()

    ax.set_xlim(0, 1.4)
    ax.set_ylim(0, 2)
    ax.axhline(0.5, color='C2', linewidth=0.3)
    sts = ['Line', '2 Lineg\n 2 Lg', '$\\sum_i x $', 'hi $\\sum_i x $\ntest',
           'test\n $\\sum_i x $', '$\\sum_i x $\n $\\sum_i x $']
    renderer = fig.canvas.get_renderer()

    def draw_box(ax, tt):
        r = mpatches.Rectangle((0, 0), 1, 1, clip_on=False,
                               transform=ax.transAxes)
        r.set_bounds(
            tt.get_window_extent(renderer)
            .transformed(ax.transAxes.inverted())
            .bounds)
        ax.add_patch(r)

    horal = 'left'
    for nn, st in enumerate(sts):
        tt = ax.text(0.2 * nn + 0.1, 0.5, st, horizontalalignment=horal,
                     verticalalignment='bottom')
        draw_box(ax, tt)
    ax.text(1.2, 0.5, 'Bottom align', color='C2')

    ax.axhline(1.3, color='C2', linewidth=0.3)
    for nn, st in enumerate(sts):
        tt = ax.text(0.2 * nn + 0.1, 1.3, st, horizontalalignment=horal,
                     verticalalignment='top')
        draw_box(ax, tt)
    ax.text(1.2, 1.3, 'Top align', color='C2')

    ax.axhline(1.8, color='C2', linewidth=0.3)
    for nn, st in enumerate(sts):
        tt = ax.text(0.2 * nn + 0.1, 1.8, st, horizontalalignment=horal,
                     verticalalignment='baseline')
        draw_box(ax, tt)
    ax.text(1.2, 1.8, 'Baseline align', color='C2')

    ax.axhline(0.1, color='C2', linewidth=0.3)
    for nn, st in enumerate(sts):
        tt = ax.text(0.2 * nn + 0.1, 0.1, st, horizontalalignment=horal,
                     verticalalignment='bottom', rotation=20)
        draw_box(ax, tt)
    ax.text(1.2, 0.1, 'Bot align, rot20', color='C2')


@image_comparison(['antialiased.png'], style='mpl20')
def test_antialiasing():
    mpl.rcParams['text.antialiased'] = False  # Passed arguments should override.

    fig = plt.figure(figsize=(5.25, 0.75))
    fig.text(0.3, 0.75, "antialiased", horizontalalignment='center',
             verticalalignment='center', antialiased=True)
    fig.text(0.3, 0.25, r"$\sqrt{x}$", horizontalalignment='center',
             verticalalignment='center', antialiased=True)

    mpl.rcParams['text.antialiased'] = True  # Passed arguments should override.
    fig.text(0.7, 0.75, "not antialiased", horizontalalignment='center',
             verticalalignment='center', antialiased=False)
    fig.text(0.7, 0.25, r"$\sqrt{x}$", horizontalalignment='center',
             verticalalignment='center', antialiased=False)

    mpl.rcParams['text.antialiased'] = False  # Should not affect existing text.


@image_comparison(['text_contains.png'], style='mpl20')
def test_contains():
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    mevent = MouseEvent('button_press_event', fig.canvas, 0.5, 0.5, 1, None)

    xs = np.linspace(0.25, 0.75, 30)
    ys = np.linspace(0.25, 0.75, 30)
    xs, ys = np.meshgrid(xs, ys)

    txt = plt.text(
        0.5, 0.4, 'hello world', ha='center', fontsize=30, rotation=30)
    # uncomment to draw the text's bounding box
    # txt.set_bbox(dict(edgecolor='black', facecolor='none'))

    # draw the text. This is important, as the contains method can only work
    # when a renderer exists.
    fig.canvas.draw()

    for x, y in zip(xs.flat, ys.flat):
        mevent.x, mevent.y = plt.gca().transAxes.transform([x, y])
        contains, _ = txt.contains(mevent)
        color = 'yellow' if contains else 'red'

        # capture the viewLim, plot a point, and reset the viewLim
        vl = ax.viewLim.frozen()
        ax.plot(x, y, 'o', color=color)
        ax.viewLim.set(vl)


def test_annotation_contains():
    # Check that Annotation.contains looks at the bboxes of the text and the
    # arrow separately, not at the joint bbox.
    fig, ax = plt.subplots()
    ann = ax.annotate(
        "hello", xy=(.4, .4), xytext=(.6, .6), arrowprops={"arrowstyle": "->"})
    fig.canvas.draw()  # Needed for the same reason as in test_contains.
    event = MouseEvent(
        "button_press_event", fig.canvas, *ax.transData.transform((.5, .6)))
    assert ann.contains(event) == (False, {})


@pytest.mark.parametrize('err, xycoords, match', (
    (TypeError, print, "xycoords callable must return a BboxBase or Transform, not a"),
    (TypeError, [0, 0], r"'xycoords' must be an instance of str, tuple"),
    (ValueError, "foo", "'foo' is not a valid coordinate"),
    (ValueError, "foo bar", "'foo bar' is not a valid coordinate"),
    (ValueError, "offset foo", "xycoords cannot be an offset coordinate"),
    (ValueError, "axes foo", "'foo' is not a recognized unit"),
))
def test_annotate_errors(err, xycoords, match):
    fig, ax = plt.subplots()
    with pytest.raises(err, match=match):
        ax.annotate('xy', (0, 0), xytext=(0.5, 0.5), xycoords=xycoords)
        fig.canvas.draw()


@image_comparison(['titles'], style='mpl20')
def test_titles():
    # left and right side titles
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.set_title("left title", loc="left")
    ax.set_title("right title", loc="right")
    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(['text_alignment'], style='mpl20')
def test_alignment():
    plt.figure()
    ax = plt.subplot(1, 1, 1)

    x = 0.1
    for rotation in (0, 30):
        for alignment in ('top', 'bottom', 'baseline', 'center'):
            ax.text(
                x, 0.5, alignment + " Tj", va=alignment, rotation=rotation,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            ax.text(
                x, 1.0, r'$\sum_{i=0}^{j}$', va=alignment, rotation=rotation)
            x += 0.1

    ax.plot([0, 1], [0.5, 0.5])
    ax.plot([0, 1], [1.0, 1.0])

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.5)
    ax.set_xticks([])
    ax.set_yticks([])


@image_comparison(baseline_images=['rotation_anchor.png'], style='mpl20',
                  remove_text=True)
def test_rotation_mode_anchor():
    fig, ax = plt.subplots()

    ax.plot([0, 1], lw=0)
    ax.axvline(.5, linewidth=.5, color='.5')
    ax.axhline(.5, linewidth=.5, color='.5')

    N = 4
    for r in range(N):
        ax.text(.5, .5, 'pP', color=f'C{r}', size=100,
                rotation=r/N*360, rotation_mode='anchor',
                verticalalignment='center_baseline')


@image_comparison(['axes_titles.png'], style='mpl20')
def test_axes_titles():
    # Related to issue #3327
    plt.figure()
    ax = plt.subplot(1, 1, 1)
    ax.set_title('center', loc='center', fontsize=20, fontweight=700)
    ax.set_title('left', loc='left', fontsize=12, fontweight=400)
    ax.set_title('right', loc='right', fontsize=12, fontweight=400)


def test_set_position():
    fig, ax = plt.subplots()

    # test set_position
    ann = ax.annotate(
        'test', (0, 0), xytext=(0, 0), textcoords='figure pixels')
    fig.canvas.draw()

    init_pos = ann.get_window_extent(fig.canvas.renderer)
    shift_val = 15
    ann.set_position((shift_val, shift_val))
    fig.canvas.draw()
    post_pos = ann.get_window_extent(fig.canvas.renderer)

    for a, b in zip(init_pos.min, post_pos.min):
        assert a + shift_val == b

    # test xyann
    ann = ax.annotate(
        'test', (0, 0), xytext=(0, 0), textcoords='figure pixels')
    fig.canvas.draw()

    init_pos = ann.get_window_extent(fig.canvas.renderer)
    shift_val = 15
    ann.xyann = (shift_val, shift_val)
    fig.canvas.draw()
    post_pos = ann.get_window_extent(fig.canvas.renderer)

    for a, b in zip(init_pos.min, post_pos.min):
        assert a + shift_val == b


def test_char_index_at():
    fig = plt.figure()
    text = fig.text(0.1, 0.9, "")

    text.set_text("i")
    bbox = text.get_window_extent()
    size_i = bbox.x1 - bbox.x0

    text.set_text("m")
    bbox = text.get_window_extent()
    size_m = bbox.x1 - bbox.x0

    text.set_text("iiiimmmm")
    bbox = text.get_window_extent()
    origin = bbox.x0

    assert text._char_index_at(origin - size_i) == 0  # left of first char
    assert text._char_index_at(origin) == 0
    assert text._char_index_at(origin + 0.499*size_i) == 0
    assert text._char_index_at(origin + 0.501*size_i) == 1
    assert text._char_index_at(origin + size_i*3) == 3
    assert text._char_index_at(origin + size_i*4 + size_m*3) == 7
    assert text._char_index_at(origin + size_i*4 + size_m*4) == 8
    assert text._char_index_at(origin + size_i*4 + size_m*10) == 8


@pytest.mark.parametrize('text', ['', 'O'], ids=['empty', 'non-empty'])
def test_non_default_dpi(text):
    fig, ax = plt.subplots()

    t1 = ax.text(0.5, 0.5, text, ha='left', va='bottom')
    fig.canvas.draw()
    dpi = fig.dpi

    bbox1 = t1.get_window_extent()
    bbox2 = t1.get_window_extent(dpi=dpi * 10)
    np.testing.assert_allclose(bbox2.get_points(), bbox1.get_points() * 10,
                               rtol=5e-2)
    # Text.get_window_extent should not permanently change dpi.
    assert fig.dpi == dpi


def test_get_rotation_string():
    assert Text(rotation='horizontal').get_rotation() == 0.
    assert Text(rotation='vertical').get_rotation() == 90.


def test_get_rotation_float():
    for i in [15., 16.70, 77.4]:
        assert Text(rotation=i).get_rotation() == i


def test_get_rotation_int():
    for i in [67, 16, 41]:
        assert Text(rotation=i).get_rotation() == float(i)


def test_get_rotation_raises():
    with pytest.raises(ValueError):
        Text(rotation='hozirontal')


def test_get_rotation_none():
    assert Text(rotation=None).get_rotation() == 0.0


def test_get_rotation_mod360():
    for i, j in zip([360., 377., 720+177.2], [0., 17., 177.2]):
        assert_almost_equal(Text(rotation=i).get_rotation(), j)


@pytest.mark.parametrize("ha", ["center", "right", "left"])
@pytest.mark.parametrize("va", ["center", "top", "bottom",
                                "baseline", "center_baseline"])
def test_null_rotation_with_rotation_mode(ha, va):
    fig, ax = plt.subplots()
    kw = dict(rotation=0, va=va, ha=ha)
    t0 = ax.text(.5, .5, 'test', rotation_mode='anchor', **kw)
    t1 = ax.text(.5, .5, 'test', rotation_mode='default', **kw)
    fig.canvas.draw()
    assert_almost_equal(t0.get_window_extent(fig.canvas.renderer).get_points(),
                        t1.get_window_extent(fig.canvas.renderer).get_points())


@image_comparison(['text_bboxclip'], style='mpl20')
def test_bbox_clipping():
    plt.text(0.9, 0.2, 'Is bbox clipped?', backgroundcolor='r', clip_on=True)
    t = plt.text(0.9, 0.5, 'Is fancy bbox clipped?', clip_on=True)
    t.set_bbox({"boxstyle": "round, pad=0.1"})


@image_comparison(['annotation_negative_ax_coords.png'], style='mpl20')
def test_annotation_negative_ax_coords():
    fig, ax = plt.subplots()

    ax.annotate('+ pts',
                xytext=[30, 20], textcoords='axes points',
                xy=[30, 20], xycoords='axes points', fontsize=32)
    ax.annotate('- pts',
                xytext=[30, -20], textcoords='axes points',
                xy=[30, -20], xycoords='axes points', fontsize=32,
                va='top')
    ax.annotate('+ frac',
                xytext=[0.75, 0.05], textcoords='axes fraction',
                xy=[0.75, 0.05], xycoords='axes fraction', fontsize=32)
    ax.annotate('- frac',
                xytext=[0.75, -0.05], textcoords='axes fraction',
                xy=[0.75, -0.05], xycoords='axes fraction', fontsize=32,
                va='top')

    ax.annotate('+ pixels',
                xytext=[160, 25], textcoords='axes pixels',
                xy=[160, 25], xycoords='axes pixels', fontsize=32)
    ax.annotate('- pixels',
                xytext=[160, -25], textcoords='axes pixels',
                xy=[160, -25], xycoords='axes pixels', fontsize=32,
                va='top')


@image_comparison(['annotation_negative_fig_coords.png'], style='mpl20')
def test_annotation_negative_fig_coords():
    fig, ax = plt.subplots()

    ax.annotate('+ pts',
                xytext=[10, 250], textcoords='figure points',
                xy=[10, 250], xycoords='figure points', fontsize=32)
    ax.annotate('- pts',
                xytext=[-10, 310], textcoords='figure points',
                xy=[-10, 310], xycoords='figure points', fontsize=32,
                va='top')
    ax.annotate('+ frac',
                xytext=[0.05, 0.5], textcoords='figure fraction',
                xy=[0.05, 0.5], xycoords='figure fraction', fontsize=32)
    ax.annotate('- frac',
                xytext=[-0.05, 0.45], textcoords='figure fraction',
                xy=[-0.05, 0.45], xycoords='figure fraction', fontsize=32,
                va='top')

    ax.annotate('+ pixels',
                xytext=[50, 50], textcoords='figure pixels',
                xy=[50, 50], xycoords='figure pixels', fontsize=32)
    ax.annotate('- pixels',
                xytext=[-50, 150], textcoords='figure pixels',
                xy=[-50, 150], xycoords='figure pixels', fontsize=32,
                va='top')


def test_text_stale():
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.draw_all()
    assert not ax1.stale
    assert not ax2.stale
    assert not fig.stale

    txt1 = ax1.text(.5, .5, 'aardvark')
    assert ax1.stale
    assert txt1.stale
    assert fig.stale

    ann1 = ax2.annotate('aardvark', xy=[.5, .5])
    assert ax2.stale
    assert ann1.stale
    assert fig.stale

    plt.draw_all()
    assert not ax1.stale
    assert not ax2.stale
    assert not fig.stale


@image_comparison(['agg_text_clip.png'], style='mpl20')
def test_agg_text_clip():
    np.random.seed(1)
    fig, (ax1, ax2) = plt.subplots(2)
    for x, y in np.random.rand(10, 2):
        ax1.text(x, y, "foo", clip_on=True)
        ax2.text(x, y, "foo")


def test_text_size_binding():
    mpl.rcParams['font.size'] = 10
    fp = mpl.font_manager.FontProperties(size='large')
    sz1 = fp.get_size_in_points()
    mpl.rcParams['font.size'] = 100

    assert sz1 == fp.get_size_in_points()


@image_comparison(['font_scaling.pdf'], style='mpl20')
def test_font_scaling():
    mpl.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(6.4, 12.4))
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.set_ylim(-10, 600)

    for i, fs in enumerate(range(4, 43, 2)):
        ax.text(0.1, i*30, f"{fs} pt font size", fontsize=fs)


@pytest.mark.parametrize('spacing1, spacing2', [(0.4, 2), (2, 0.4), (2, 2)])
def test_two_2line_texts(spacing1, spacing2):
    text_string = 'line1\nline2'
    fig = plt.figure()
    renderer = fig.canvas.get_renderer()

    text1 = fig.text(0.25, 0.5, text_string, linespacing=spacing1)
    text2 = fig.text(0.25, 0.5, text_string, linespacing=spacing2)
    fig.canvas.draw()

    box1 = text1.get_window_extent(renderer=renderer)
    box2 = text2.get_window_extent(renderer=renderer)

    # line spacing only affects height
    assert box1.width == box2.width
    if spacing1 == spacing2:
        assert box1.height == box2.height
    else:
        assert box1.height != box2.height


def test_validate_linespacing():
    with pytest.raises(TypeError):
        plt.text(.25, .5, "foo", linespacing="abc")


def test_nonfinite_pos():
    fig, ax = plt.subplots()
    ax.text(0, np.nan, 'nan')
    ax.text(np.inf, 0, 'inf')
    fig.canvas.draw()


@needs_usetex
def test_usetex_is_copied():
    # Indirectly tests that update_from (which is used to copy tick label
    # properties) copies usetex state.
    fig = plt.figure()
    plt.rcParams["text.usetex"] = False
    ax1 = fig.add_subplot(121)
    plt.rcParams["text.usetex"] = True
    ax2 = fig.add_subplot(122)
    fig.canvas.draw()
    for ax, usetex in [(ax1, False), (ax2, True)]:
        for t in ax.xaxis.majorTicks:
            assert t.label1.get_usetex() == usetex


@needs_usetex
def test_single_artist_usetex():
    # Check that a single artist marked with usetex does not get passed through
    # the mathtext parser at all (for the Agg backend) (the mathtext parser
    # currently fails to parse \frac12, requiring \frac{1}{2} instead).
    fig = plt.figure()
    fig.text(.5, .5, r"$\frac12$", usetex=True)
    fig.canvas.draw()


@pytest.mark.parametrize("fmt", ["png", "pdf", "svg"])
def test_single_artist_usenotex(fmt):
    # Check that a single artist can be marked as not-usetex even though the
    # rcParam is on ("2_2_2" fails if passed to TeX).  This currently skips
    # postscript output as the ps renderer doesn't support mixing usetex and
    # non-usetex.
    plt.rcParams["text.usetex"] = True
    fig = plt.figure()
    fig.text(.5, .5, "2_2_2", usetex=False)
    fig.savefig(io.BytesIO(), format=fmt)


@image_comparison(['text_as_path_opacity.svg'], style='_classic_test')
def test_text_as_path_opacity():
    plt.figure()
    plt.gca().set_axis_off()
    plt.text(0.25, 0.25, 'c', color=(0, 0, 0, 0.5))
    plt.text(0.25, 0.5, 'a', alpha=0.5)
    plt.text(0.25, 0.75, 'x', alpha=0.5, color=(0, 0, 0, 1))


@image_comparison(['text_as_text_opacity.svg'], style='_classic_test')
def test_text_as_text_opacity():
    mpl.rcParams['svg.fonttype'] = 'none'
    plt.figure()
    plt.gca().set_axis_off()
    plt.text(0.25, 0.25, '50% using `color`', color=(0, 0, 0, 0.5))
    plt.text(0.25, 0.5, '50% using `alpha`', alpha=0.5)
    plt.text(0.25, 0.75, '50% using `alpha` and 100% `color`', alpha=0.5,
             color=(0, 0, 0, 1))


def test_text_repr():
    # smoketest to make sure text repr doesn't error for category
    plt.plot(['A', 'B'], [1, 2])
    repr(plt.text(['A'], 0.5, 'Boo'))


def test_annotation_update():
    fig, ax = plt.subplots(1, 1)
    an = ax.annotate('annotation', xy=(0.5, 0.5))
    extent1 = an.get_window_extent(fig.canvas.get_renderer())
    fig.tight_layout()
    extent2 = an.get_window_extent(fig.canvas.get_renderer())

    assert not np.allclose(extent1.get_points(), extent2.get_points(),
                           rtol=1e-6)


@check_figures_equal()
def test_annotation_units(fig_test, fig_ref):
    ax = fig_test.add_subplot()
    ax.plot(datetime.now(), 1, "o")  # Implicitly set axes extents.
    ax.annotate("x", (datetime.now(), 0.5), xycoords=("data", "axes fraction"),
                # This used to crash before.
                xytext=(0, 0), textcoords="offset points")
    ax = fig_ref.add_subplot()
    ax.plot(datetime.now(), 1, "o")
    ax.annotate("x", (datetime.now(), 0.5), xycoords=("data", "axes fraction"))


@image_comparison(['large_subscript_title.png'], style='mpl20')
def test_large_subscript_title():
    plt.rcParams['axes.titley'] = None

    fig, axs = plt.subplots(1, 2, figsize=(9, 2.5), constrained_layout=True)
    ax = axs[0]
    ax.set_title(r'$\sum_{i} x_i$')
    ax.set_title('New way', loc='left')
    ax.set_xticklabels([])

    ax = axs[1]
    ax.set_title(r'$\sum_{i} x_i$', y=1.01)
    ax.set_title('Old Way', loc='left')
    ax.set_xticklabels([])


@pytest.mark.parametrize(
    "x, rotation, halign",
    [(0.7, 0, 'left'),
     (0.5, 95, 'left'),
     (0.3, 0, 'right'),
     (0.3, 185, 'left')])
def test_wrap(x, rotation, halign):
    fig = plt.figure(figsize=(18, 18))
    gs = GridSpec(nrows=3, ncols=3, figure=fig)
    subfig = fig.add_subfigure(gs[1, 1])
    # we only use the central subfigure, which does not align with any
    # figure boundary, to ensure only subfigure boundaries are relevant
    s = 'This is a very long text that should be wrapped multiple times.'
    text = subfig.text(x, 0.7, s, wrap=True, rotation=rotation, ha=halign)
    fig.canvas.draw()
    assert text._get_wrapped_text() == ('This is a very long\n'
                                        'text that should be\n'
                                        'wrapped multiple\n'
                                        'times.')


def test_mathwrap():
    fig = plt.figure(figsize=(5, 4))
    s = r'This is a very $\overline{\mathrm{long}}$ line of Mathtext.'
    text = fig.text(0, 0.5, s, size=40, wrap=True)
    fig.canvas.draw()
    assert text._get_wrapped_text() == ('This is a very $\\overline{\\mathrm{long}}$\n'
                                        'line of Mathtext.')


def test_get_window_extent_wrapped():
    # Test that a long title that wraps to two lines has the same vertical
    # extent as an explicit two line title.

    fig1 = plt.figure(figsize=(3, 3))
    fig1.suptitle("suptitle that is clearly too long in this case", wrap=True)
    window_extent_test = fig1._suptitle.get_window_extent()

    fig2 = plt.figure(figsize=(3, 3))
    fig2.suptitle("suptitle that is clearly\ntoo long in this case")
    window_extent_ref = fig2._suptitle.get_window_extent()

    assert window_extent_test.y0 == window_extent_ref.y0
    assert window_extent_test.y1 == window_extent_ref.y1


def test_long_word_wrap():
    fig = plt.figure(figsize=(6, 4))
    text = fig.text(9.5, 8, 'Alonglineoftexttowrap', wrap=True)
    fig.canvas.draw()
    assert text._get_wrapped_text() == 'Alonglineoftexttowrap'


def test_wrap_no_wrap():
    fig = plt.figure(figsize=(6, 4))
    text = fig.text(0, 0, 'non wrapped text', wrap=True)
    fig.canvas.draw()
    assert text._get_wrapped_text() == 'non wrapped text'


@check_figures_equal()
def test_buffer_size(fig_test, fig_ref):
    # On old versions of the Agg renderer, large non-ascii single-character
    # strings (here, "€") would be rendered clipped because the rendering
    # buffer would be set by the physical size of the smaller "a" character.
    ax = fig_test.add_subplot()
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["€", "a"])
    ax.yaxis.majorTicks[1].label1.set_color("w")
    ax = fig_ref.add_subplot()
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["€", ""])


def test_fontproperties_kwarg_precedence():
    """Test that kwargs take precedence over fontproperties defaults."""
    plt.figure()
    text1 = plt.xlabel("value", fontproperties='Times New Roman', size=40.0)
    text2 = plt.ylabel("counts", size=40.0, fontproperties='Times New Roman')
    assert text1.get_size() == 40.0
    assert text2.get_size() == 40.0


def test_transform_rotates_text():
    ax = plt.gca()
    transform = mtransforms.Affine2D().rotate_deg(30)
    text = ax.text(0, 0, 'test', transform=transform,
                   transform_rotates_text=True)
    result = text.get_rotation()
    assert_almost_equal(result, 30)


def test_update_mutate_input():
    inp = dict(fontproperties=FontProperties(weight="bold"),
               bbox=None)
    cache = dict(inp)
    t = Text()
    t.update(inp)
    assert inp['fontproperties'] == cache['fontproperties']
    assert inp['bbox'] == cache['bbox']


@pytest.mark.parametrize('rotation', ['invalid string', [90]])
def test_invalid_rotation_values(rotation):
    with pytest.raises(
            ValueError,
            match=("rotation must be 'vertical', 'horizontal' or a number")):
        Text(0, 0, 'foo', rotation=rotation)


def test_invalid_color():
    with pytest.raises(ValueError):
        plt.figtext(.5, .5, "foo", c="foobar")


# See gh-26152 for more information on this xfail
@pytest.mark.xfail(pyparsing_version.release == (3, 1, 0),
                   reason="Error messages are incorrect with pyparsing 3.1.0")
def test_parse_math():
    fig, ax = plt.subplots()
    ax.text(0, 0, r"$ \wrong{math} $", parse_math=False)
    fig.canvas.draw()

    ax.text(0, 0, r"$ \wrong{math} $", parse_math=True)
    with pytest.raises(ValueError, match='Unknown symbol'):
        fig.canvas.draw()


# See gh-26152 for more information on this xfail
@pytest.mark.xfail(pyparsing_version.release == (3, 1, 0),
                   reason="Error messages are incorrect with pyparsing 3.1.0")
def test_parse_math_rcparams():
    # Default is True
    fig, ax = plt.subplots()
    ax.text(0, 0, r"$ \wrong{math} $")
    with pytest.raises(ValueError, match='Unknown symbol'):
        fig.canvas.draw()

    # Setting rcParams to False
    with mpl.rc_context({'text.parse_math': False}):
        fig, ax = plt.subplots()
        ax.text(0, 0, r"$ \wrong{math} $")
        fig.canvas.draw()


@needs_usetex
def test_metrics_cache():
    # dig into the signature to get the mutable default used as a cache
    renderer_cache = inspect.signature(
        mpl.text._get_text_metrics_function
    ).parameters['_cache'].default

    renderer_cache.clear()

    fig = plt.figure()
    fig.text(.3, .5, "foo\nbar")
    fig.text(.3, .5, "foo\nbar", usetex=True)
    fig.text(.5, .5, "foo\nbar", usetex=True)
    fig.canvas.draw()
    renderer = fig._get_renderer()
    assert renderer in renderer_cache

    ys = {}  # mapping of strings to where they were drawn in y with draw_tex.

    def call(*args, **kwargs):
        renderer, x, y, s, *_ = args
        ys.setdefault(s, set()).add(y)

    renderer.draw_tex = call
    fig.canvas.draw()
    assert [*ys] == ["foo", "bar"]
    # Check that both TeX strings were drawn with the same y-position for both
    # single-line substrings.  Previously, there used to be an incorrect cache
    # collision with the non-TeX string (drawn first here) whose metrics would
    # get incorrectly reused by the first TeX string.
    assert len(ys["foo"]) == len(ys["bar"]) == 1

    info = renderer_cache[renderer].cache_info()
    # Every string gets a miss for the first layouting (extents), then a hit
    # when drawing, but "foo\nbar" gets two hits as it's drawn twice.
    assert info.hits > info.misses


def test_metrics_cache2():
    # dig into the signature to get the mutable default used as a cache
    renderer_cache = inspect.signature(
        mpl.text._get_text_metrics_function
    ).parameters['_cache'].default
    gc.collect()
    renderer_cache.clear()

    def helper():
        fig, ax = plt.subplots()
        fig.draw_without_rendering()
        # show we hit the outer cache
        assert len(renderer_cache) == 1
        func = renderer_cache[fig.canvas.get_renderer()]
        cache_info = func.cache_info()
        # show we hit the inner cache
        assert cache_info.currsize > 0
        assert cache_info.currsize == cache_info.misses
        assert cache_info.hits > cache_info.misses
        plt.close(fig)

    helper()
    gc.collect()
    # show the outer cache has a lifetime tied to the renderer (via the figure)
    assert len(renderer_cache) == 0


def test_annotate_offset_fontsize():
    # Test that offset_fontsize parameter works and uses accurate values
    fig, ax = plt.subplots()
    text_coords = ['offset points', 'offset fontsize']
    # 10 points should be equal to 1 fontsize unit at fontsize=10
    xy_text = [(10, 10), (1, 1)]
    anns = [ax.annotate('test', xy=(0.5, 0.5),
                        xytext=xy_text[i],
                        fontsize='10',
                        xycoords='data',
                        textcoords=text_coords[i]) for i in range(2)]
    points_coords, fontsize_coords = (ann.get_window_extent() for ann in anns)
    fig.canvas.draw()
    assert str(points_coords) == str(fontsize_coords)


def test_get_set_antialiased():
    txt = Text(.5, .5, "foo\nbar")
    assert txt._antialiased == mpl.rcParams['text.antialiased']
    assert txt.get_antialiased() == mpl.rcParams['text.antialiased']

    txt.set_antialiased(True)
    assert txt._antialiased is True
    assert txt.get_antialiased() == txt._antialiased

    txt.set_antialiased(False)
    assert txt._antialiased is False
    assert txt.get_antialiased() == txt._antialiased


def test_annotation_antialiased():
    annot = Annotation("foo\nbar", (.5, .5), antialiased=True)
    assert annot._antialiased is True
    assert annot.get_antialiased() == annot._antialiased

    annot2 = Annotation("foo\nbar", (.5, .5), antialiased=False)
    assert annot2._antialiased is False
    assert annot2.get_antialiased() == annot2._antialiased

    annot3 = Annotation("foo\nbar", (.5, .5), antialiased=False)
    annot3.set_antialiased(True)
    assert annot3.get_antialiased() is True
    assert annot3._antialiased is True

    annot4 = Annotation("foo\nbar", (.5, .5))
    assert annot4._antialiased == mpl.rcParams['text.antialiased']


@check_figures_equal()
def test_annotate_and_offsetfrom_copy_input(fig_test, fig_ref):
    # Both approaches place the text (10, 0) pixels away from the center of the line.
    ax = fig_test.add_subplot()
    l, = ax.plot([0, 2], [0, 2])
    of_xy = np.array([.5, .5])
    ax.annotate("foo", textcoords=OffsetFrom(l, of_xy), xytext=(10, 0),
                xy=(0, 0))  # xy is unused.
    of_xy[:] = 1
    ax = fig_ref.add_subplot()
    l, = ax.plot([0, 2], [0, 2])
    an_xy = np.array([.5, .5])
    ax.annotate("foo", xy=an_xy, xycoords=l, xytext=(10, 0), textcoords="offset points")
    an_xy[:] = 2


@check_figures_equal(extensions=['png', 'pdf', 'svg'])
def test_text_antialiased_off_default_vs_manual(fig_test, fig_ref):
    fig_test.text(0.5, 0.5, '6 inches x 2 inches',
                             antialiased=False)

    mpl.rcParams['text.antialiased'] = False
    fig_ref.text(0.5, 0.5, '6 inches x 2 inches')


@check_figures_equal(extensions=['png', 'pdf', 'svg'])
def test_text_antialiased_on_default_vs_manual(fig_test, fig_ref):
    fig_test.text(0.5, 0.5, '6 inches x 2 inches', antialiased=True)

    mpl.rcParams['text.antialiased'] = True
    fig_ref.text(0.5, 0.5, '6 inches x 2 inches')


def test_text_annotation_get_window_extent():
    figure = Figure(dpi=100)
    renderer = RendererAgg(200, 200, 100)

    # Only text annotation
    annotation = Annotation('test', xy=(0, 0), xycoords='figure pixels')
    annotation.set_figure(figure)

    text = Text(text='test', x=0, y=0)
    text.set_figure(figure)

    bbox = annotation.get_window_extent(renderer=renderer)

    text_bbox = text.get_window_extent(renderer=renderer)
    assert bbox.width == text_bbox.width
    assert bbox.height == text_bbox.height

    _, _, d = renderer.get_text_width_height_descent(
        'text', annotation._fontproperties, ismath=False)
    font = get_font(fontManager._find_fonts_by_props(annotation._fontproperties))
    for name, key in [('OS/2', 'sTypoDescender'), ('hhea', 'descent')]:
        if (table := font.get_sfnt_table(name)) is not None:
            units_per_em = font.get_sfnt_table('head')['unitsPerEm']
            fontsize = annotation._fontproperties.get_size_in_points()
            lp_d = -table[key] / units_per_em * fontsize * figure.dpi / 72
            break
    else:
        _, _, lp_d = renderer.get_text_width_height_descent(
            'lp', annotation._fontproperties, ismath=False)
    below_line = max(d, lp_d)

    # These numbers are specific to the current implementation of Text
    points = bbox.get_points()
    assert points[0, 0] == 0.0
    assert points[1, 0] == text_bbox.width
    assert points[0, 1] == -below_line
    assert points[1, 1] == text_bbox.height - below_line


def test_text_with_arrow_annotation_get_window_extent():
    headwidth = 21
    fig, ax = plt.subplots(dpi=100)
    txt = ax.text(s='test', x=0, y=0)
    ann = ax.annotate(
        'test',
        xy=(0.0, 50.0),
        xytext=(50.0, 50.0), xycoords='figure pixels',
        arrowprops={
            'facecolor': 'black', 'width': 2,
            'headwidth': headwidth, 'shrink': 0.0})

    plt.draw()
    renderer = fig.canvas.renderer
    # bounding box of text
    text_bbox = txt.get_window_extent(renderer=renderer)
    # bounding box of annotation (text + arrow)
    bbox = ann.get_window_extent(renderer=renderer)
    # bounding box of arrow
    arrow_bbox = ann.arrow_patch.get_window_extent(renderer)
    # bounding box of annotation text
    ann_txt_bbox = Text.get_window_extent(ann)

    # make sure annotation width is 50 px wider than
    # just the text
    assert bbox.width == text_bbox.width + 50.0
    # make sure the annotation text bounding box is same size
    # as the bounding box of the same string as a Text object
    assert_almost_equal(ann_txt_bbox.height, text_bbox.height)
    assert ann_txt_bbox.width == text_bbox.width
    # compute the expected bounding box of arrow + text
    expected_bbox = mtransforms.Bbox.union([ann_txt_bbox, arrow_bbox])
    assert_almost_equal(bbox.height, expected_bbox.height)


def test_arrow_annotation_get_window_extent():
    dpi = 100
    dots_per_point = dpi / 72
    figure = Figure(dpi=dpi)
    figure.set_figwidth(2.0)
    figure.set_figheight(2.0)
    renderer = RendererAgg(200, 200, 100)

    # Text annotation with arrow; arrow dimensions are in points
    annotation = Annotation(
        '', xy=(0.0, 50.0), xytext=(50.0, 50.0), xycoords='figure pixels',
        arrowprops={
            'facecolor': 'black', 'width': 8, 'headwidth': 10, 'shrink': 0.0})
    annotation.set_figure(figure)
    annotation.draw(renderer)

    bbox = annotation.get_window_extent()
    points = bbox.get_points()

    assert bbox.width == 50.0
    assert_almost_equal(bbox.height, 10.0 * dots_per_point)
    assert points[0, 0] == 0.0
    assert points[0, 1] == 50.0 - 5 * dots_per_point


def test_empty_annotation_get_window_extent():
    figure = Figure(dpi=100)
    figure.set_figwidth(2.0)
    figure.set_figheight(2.0)
    renderer = RendererAgg(200, 200, 100)

    # Text annotation with arrow
    annotation = Annotation(
        '', xy=(0.0, 50.0), xytext=(0.0, 50.0), xycoords='figure pixels')
    annotation.set_figure(figure)
    annotation.draw(renderer)

    bbox = annotation.get_window_extent()
    points = bbox.get_points()

    assert points[0, 0] == 0.0
    assert points[1, 0] == 0.0
    assert points[1, 1] == 50.0
    assert points[0, 1] == 50.0


@image_comparison(['basictext_wrap.png'], style='mpl20')
def test_basic_wrap():
    fig = plt.figure(figsize=(8, 6))
    plt.axis([0, 10, 0, 10])
    t = "This is a really long string that I'd rather have wrapped so that" \
        " it doesn't go outside of the figure, but if it's long enough it" \
        " will go off the top or bottom!"
    plt.text(4, 1, t, ha='left', rotation=15, wrap=True)
    plt.text(6, 5, t, ha='left', rotation=15, wrap=True)
    plt.text(5, 5, t, ha='right', rotation=-15, wrap=True)
    plt.text(5, 10, t, fontsize=18, style='oblique', ha='center',
             va='top', wrap=True)
    plt.text(3, 4, t, family='serif', style='italic', ha='right', wrap=True)
    plt.text(-1, 0, t, ha='left', rotation=-15, wrap=True)


@image_comparison(['fonttext_wrap.png'], style='mpl20')
def test_font_wrap():
    fig = plt.figure(figsize=(8, 6))
    plt.axis([0, 10, 0, 10])
    t = "This is a really long string that I'd rather have wrapped so that" \
        " it doesn't go outside of the figure, but if it's long enough it" \
        " will go off the top or bottom!"
    plt.text(4, -1, t, fontsize=18, family='serif', ha='left', rotation=15,
             wrap=True)
    plt.text(6, 5, t, family='sans serif', ha='left', rotation=15, wrap=True)
    plt.text(5, 10, t, weight='heavy', ha='center', va='top', wrap=True)
    plt.text(3, 4, t, family='monospace', ha='right', wrap=True)
    plt.text(-1, 0, t, fontsize=14, style='italic', ha='left', rotation=-15,
             wrap=True)


def test_ha_for_angle():
    text_instance = Text()
    angles = np.arange(0, 360.1, 0.1)
    for angle in angles:
        alignment = text_instance._ha_for_angle(angle)
        assert alignment in ['center', 'left', 'right']


def test_va_for_angle():
    text_instance = Text()
    angles = np.arange(0, 360.1, 0.1)
    for angle in angles:
        alignment = text_instance._va_for_angle(angle)
        assert alignment in ['center', 'top', 'baseline']


@image_comparison(['xtick_rotation_mode.png'], remove_text=False, style='mpl20')
def test_xtick_rotation_mode():
    fig, ax = plt.subplots(figsize=(12, 1))
    ax.set_yticks([])
    ax2 = ax.twiny()

    ax.set_xticks(range(37), ['foo'] * 37, rotation_mode="xtick")
    ax2.set_xticks(range(37), ['foo'] * 37, rotation_mode="xtick")

    angles = np.linspace(0, 360, 37)

    for tick, angle in zip(ax.get_xticklabels(), angles):
        tick.set_rotation(angle)
    for tick, angle in zip(ax2.get_xticklabels(), angles):
        tick.set_rotation(angle)

    plt.subplots_adjust(left=0.01, right=0.99, top=.6, bottom=.4)


@image_comparison(['ytick_rotation_mode.png'], remove_text=False, style='mpl20')
def test_ytick_rotation_mode():
    fig, ax = plt.subplots(figsize=(1, 12))
    ax.set_xticks([])
    ax2 = ax.twinx()

    ax.set_yticks(range(37), ['foo'] * 37, rotation_mode="ytick")
    ax2.set_yticks(range(37), ['foo'] * 37, rotation_mode='ytick')

    angles = np.linspace(0, 360, 37)
    for tick, angle in zip(ax.get_yticklabels(), angles):
        tick.set_rotation(angle)
    for tick, angle in zip(ax2.get_yticklabels(), angles):
        tick.set_rotation(angle)

    plt.subplots_adjust(left=0.4, right=0.6, top=.99, bottom=.01)


def test_text_tightbbox_outside_scale_domain():
    # Test that text at positions outside the valid domain of axes scales
    # (e.g., negative coordinates with log scale) returns a null bbox.
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_ylim(1, 100)

    invalid_text = ax.text(0, -5, 'invalid')
    invalid_bbox = invalid_text.get_tightbbox(fig.canvas.get_renderer())
    assert not np.isfinite(invalid_bbox.width)


def _test_complex_shaping(fig):
    # Raqm is Arabic for writing; note that because Arabic is RTL, the characters here
    # may seem to be in a different order than expected, but libraqm will order them
    # correctly for us.
    text = (
        'Arabic: \N{Arabic Letter REH}\N{Arabic FATHA}\N{Arabic Letter QAF}'
        '\N{Arabic SUKUN}\N{Arabic Letter MEEM}')
    math_signs = '\N{N-ary Product}\N{N-ary Coproduct}\N{N-ary summation}\N{Integral}'
    text = math_signs + text + math_signs
    fig.text(0.5, 0.75, text, size=32, ha='center', va='center')
    # Also check fallback behaviour:
    # - English should use cmr10
    # - Math signs should use DejaVu Sans Display (and thus be larger than the rest)
    # - Arabic should use DejaVu Sans
    fig.text(0.5, 0.25, text, size=32, ha='center', va='center',
             family=['cmr10', 'DejaVu Sans Display', 'DejaVu Sans'])


@image_comparison(['complex'], extensions=['png', 'pdf', 'svg', 'eps'], style='mpl20')
def test_complex_shaping():
    fig = plt.figure(figsize=(6, 2))
    _test_complex_shaping(fig)


def _test_text_features(fig):
    t = fig.text(1, 0.7, 'Default: fi ffi fl st',
                 fontsize=32, horizontalalignment='right')
    assert t.get_fontfeatures() is None
    t = fig.text(1, 0.4, 'Disabled: fi ffi fl st',
                 fontsize=32, horizontalalignment='right',
                 fontfeatures=['-liga'])
    assert t.get_fontfeatures() == ('-liga', )
    t = fig.text(1, 0.1, 'Discretionary: fi ffi fl st',
                 fontsize=32, horizontalalignment='right')
    t.set_fontfeatures(['dlig'])
    assert t.get_fontfeatures() == ('dlig', )


@image_comparison(['features'], remove_text=False, style='mpl20',
                  extensions=['png', 'pdf', 'svg', 'eps'])
def test_text_features():
    fig = plt.figure(figsize=(5, 1.5))
    _test_text_features(fig)


@pytest.mark.parametrize(
    'input, match',
    [
        ([1, 2, 3], 'must be list of tuple'),
        ([(1, 2)], 'must be list of tuple'),
        ([('en', 'foo', 2)], 'start location must be int'),
        ([('en', 1, 'foo')], 'end location must be int'),
    ],
)
def test_text_language_invalid(input, match):
    with pytest.raises(TypeError, match=match):
        Text(0, 0, 'foo', language=input)


def _test_text_language(fig):
    t = fig.text(0, 0.8, 'Default', fontsize=32)
    assert t.get_language() is None
    t = fig.text(0, 0.55, 'Lang A', fontsize=32)
    assert t.get_language() is None
    t = fig.text(0, 0.3, 'Lang B', fontsize=32)
    assert t.get_language() is None
    t = fig.text(0, 0.05, 'Mixed', fontsize=32)
    assert t.get_language() is None

    # DejaVu Sans supports language-specific glyphs in the Serbian and Macedonian
    # languages in the Cyrillic alphabet.
    cyrillic = '\U00000431'
    t = fig.text(0.4, 0.8, cyrillic, fontsize=32)
    assert t.get_language() is None
    t = fig.text(0.4, 0.55, cyrillic, fontsize=32, language='sr')
    assert t.get_language() == 'sr'
    t = fig.text(0.4, 0.3, cyrillic, fontsize=32)
    t.set_language('ru')
    assert t.get_language() == 'ru'
    t = fig.text(0.4, 0.05, cyrillic * 4, fontsize=32,
                 language=[('ru', 0, 1), ('sr', 1, 2), ('ru', 2, 3), ('sr', 3, 4)])
    assert t.get_language() == (('ru', 0, 1), ('sr', 1, 2), ('ru', 2, 3), ('sr', 3, 4))

    # Or the Sámi family of languages in the Latin alphabet.
    latin = '\U0000014a'
    t = fig.text(0.7, 0.8, latin, fontsize=32)
    assert t.get_language() is None
    with plt.rc_context({'text.language': 'en'}):
        t = fig.text(0.7, 0.55, latin, fontsize=32)
    assert t.get_language() == 'en'
    t = fig.text(0.7, 0.3, latin, fontsize=32, language='smn')
    assert t.get_language() == 'smn'
    # Tuples are not documented, but we'll allow it.
    t = fig.text(0.7, 0.05, latin * 4, fontsize=32)
    t.set_language((('en', 0, 1), ('smn', 1, 2), ('en', 2, 3), ('smn', 3, 4)))
    assert t.get_language() == (
        ('en', 0, 1), ('smn', 1, 2), ('en', 2, 3), ('smn', 3, 4))


@image_comparison(['language'], remove_text=False, style='mpl20',
                  extensions=['png', 'pdf', 'svg', 'eps'])
def test_text_language():
    fig = plt.figure(figsize=(5, 3))
    _test_text_language(fig)


@image_comparison(['draw_text_fallback.png'], style='mpl20')
def test_draw_text_as_path_fallback(monkeypatch):
    # Delete RendererAgg.draw_text so that we use the RendererBase.draw_text fallback.
    monkeypatch.delattr('matplotlib.backends.backend_agg.RendererAgg.draw_text')
    heights = [2, 1.5, 3]
    fig = plt.figure(figsize=(6, sum(heights)))
    subfig = fig.subfigures(3, 1, height_ratios=heights)
    _test_complex_shaping(subfig[0])
    _test_text_features(subfig[1])
    _test_text_language(subfig[2])
