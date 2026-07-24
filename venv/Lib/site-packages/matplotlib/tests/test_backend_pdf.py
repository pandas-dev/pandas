import datetime
import decimal
import io
from pathlib import Path
import string

import numpy as np
from packaging.version import parse as parse_version
import pytest

import matplotlib as mpl
from matplotlib import (
    pyplot as plt, rcParams, font_manager as fm
)
from matplotlib.cbook import _get_data_path
from matplotlib.ft2font import FT2Font
from matplotlib.backends._backend_pdf_ps import get_glyphs_subset, font_as_file
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.testing import _gen_multi_font_text, _has_tex_package
from matplotlib.testing.decorators import check_figures_equal, image_comparison
from matplotlib.testing._markers import needs_usetex


@image_comparison(['pdf_use14corefonts.pdf'], style='_classic_test')
def test_use14corefonts():
    rcParams['pdf.use14corefonts'] = True
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.size'] = 8
    rcParams['font.sans-serif'] = ['Helvetica']
    rcParams['pdf.compression'] = 0

    text = '''A three-line text positioned just above a blue line
and containing some French characters and the euro symbol:
"Merci pépé pour les 10 €"'''

    fig, ax = plt.subplots()
    ax.set_title('Test PDF backend with option use14corefonts=True')
    ax.text(0.5, 0.5, text, horizontalalignment='center',
            verticalalignment='bottom',
            fontsize=14)
    ax.axhline(0.5, linewidth=0.5)


def test_multipage_pagecount():
    with PdfPages(io.BytesIO()) as pdf:
        assert pdf.get_pagecount() == 0
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3])
        fig.savefig(pdf, format="pdf")
        assert pdf.get_pagecount() == 1
        pdf.savefig()
        assert pdf.get_pagecount() == 2


def test_multipage_properfinalize():
    pdfio = io.BytesIO()
    with PdfPages(pdfio) as pdf:
        for i in range(10):
            fig, ax = plt.subplots()
            ax.set_title('This is a long title')
            fig.savefig(pdf, format="pdf")
    s = pdfio.getvalue()
    assert s.count(b'startxref') == 1
    assert len(s) < 40000


def test_multipage_keep_empty(tmp_path):
    # An empty pdf deletes itself afterwards.
    fn = tmp_path / "a.pdf"
    with PdfPages(fn) as pdf:
        pass
    assert not fn.exists()

    # Test pdf files with content, they should never be deleted.
    fn = tmp_path / "b.pdf"
    with PdfPages(fn) as pdf:
        pdf.savefig(plt.figure())
    assert fn.exists()


def test_composite_image():
    # Test that figures can be saved with and without combining multiple images
    # (on a single set of axes) into a single composite image.
    X, Y = np.meshgrid(np.arange(-5, 5, 1), np.arange(-5, 5, 1))
    Z = np.sin(Y ** 2)
    fig, ax = plt.subplots()
    ax.set_xlim(0, 3)
    ax.imshow(Z, extent=[0, 1, 0, 1])
    ax.imshow(Z[::-1], extent=[2, 3, 0, 1])
    plt.rcParams['image.composite_image'] = True
    with PdfPages(io.BytesIO()) as pdf:
        fig.savefig(pdf, format="pdf")
        assert len(pdf._file._images) == 1
    plt.rcParams['image.composite_image'] = False
    with PdfPages(io.BytesIO()) as pdf:
        fig.savefig(pdf, format="pdf")
        assert len(pdf._file._images) == 2


def test_indexed_image():
    # An image with low color count should compress to a palette-indexed format.
    pikepdf = pytest.importorskip('pikepdf')

    data = np.zeros((256, 1, 3), dtype=np.uint8)
    data[:, 0, 0] = np.arange(256)  # Maximum unique colours for an indexed image.

    rcParams['pdf.compression'] = True
    fig = plt.figure()
    fig.figimage(data, resize=True)
    buf = io.BytesIO()
    fig.savefig(buf, format='pdf', dpi='figure')

    with pikepdf.Pdf.open(buf) as pdf:
        page, = pdf.pages
        if parse_version(pikepdf.__version__) < parse_version('10.9.0'):
            image, = page.images.values()
        else:
            image, = page.get_images().values()
        pdf_image = pikepdf.PdfImage(image)
        assert pdf_image.indexed
        pil_image = pdf_image.as_pil_image()
        rgb = np.asarray(pil_image.convert('RGB'))

    np.testing.assert_array_equal(data, rgb)


def test_savefig_metadata(monkeypatch):
    pikepdf = pytest.importorskip('pikepdf')
    monkeypatch.setenv('SOURCE_DATE_EPOCH', '0')

    fig, ax = plt.subplots()
    ax.plot(range(5))

    md = {
        'Author': 'me',
        'Title': 'Multipage PDF',
        'Subject': 'Test page',
        'Keywords': 'test,pdf,multipage',
        'ModDate': datetime.datetime(
            1968, 8, 1, tzinfo=datetime.timezone(datetime.timedelta(0))),
        'Trapped': 'True'
    }
    buf = io.BytesIO()
    fig.savefig(buf, metadata=md, format='pdf')

    with pikepdf.Pdf.open(buf) as pdf:
        info = {k: str(v) for k, v in pdf.docinfo.items()}

    assert info == {
        '/Author': 'me',
        '/CreationDate': 'D:19700101000000Z',
        '/Creator': f'Matplotlib v{mpl.__version__}, https://matplotlib.org',
        '/Keywords': 'test,pdf,multipage',
        '/ModDate': 'D:19680801000000Z',
        '/Producer': f'Matplotlib pdf backend v{mpl.__version__}',
        '/Subject': 'Test page',
        '/Title': 'Multipage PDF',
        '/Trapped': '/True',
    }


def test_invalid_metadata():
    fig, ax = plt.subplots()

    with pytest.warns(UserWarning,
                      match="Unknown infodict keyword: 'foobar'."):
        fig.savefig(io.BytesIO(), format='pdf', metadata={'foobar': 'invalid'})

    with pytest.warns(UserWarning,
                      match='not an instance of datetime.datetime.'):
        fig.savefig(io.BytesIO(), format='pdf',
                    metadata={'ModDate': '1968-08-01'})

    with pytest.warns(UserWarning,
                      match='not one of {"True", "False", "Unknown"}'):
        fig.savefig(io.BytesIO(), format='pdf', metadata={'Trapped': 'foo'})

    with pytest.warns(UserWarning, match='not an instance of str.'):
        fig.savefig(io.BytesIO(), format='pdf', metadata={'Title': 1234})


def test_multipage_metadata(monkeypatch):
    pikepdf = pytest.importorskip('pikepdf')
    monkeypatch.setenv('SOURCE_DATE_EPOCH', '0')

    fig, ax = plt.subplots()
    ax.plot(range(5))

    md = {
        'Author': 'me',
        'Title': 'Multipage PDF',
        'Subject': 'Test page',
        'Keywords': 'test,pdf,multipage',
        'ModDate': datetime.datetime(
            1968, 8, 1, tzinfo=datetime.timezone(datetime.timedelta(0))),
        'Trapped': 'True'
    }
    buf = io.BytesIO()
    with PdfPages(buf, metadata=md) as pdf:
        pdf.savefig(fig)
        pdf.savefig(fig)

    with pikepdf.Pdf.open(buf) as pdf:
        info = {k: str(v) for k, v in pdf.docinfo.items()}

    assert info == {
        '/Author': 'me',
        '/CreationDate': 'D:19700101000000Z',
        '/Creator': f'Matplotlib v{mpl.__version__}, https://matplotlib.org',
        '/Keywords': 'test,pdf,multipage',
        '/ModDate': 'D:19680801000000Z',
        '/Producer': f'Matplotlib pdf backend v{mpl.__version__}',
        '/Subject': 'Test page',
        '/Title': 'Multipage PDF',
        '/Trapped': '/True',
    }


def test_text_urls():
    pikepdf = pytest.importorskip('pikepdf')

    test_url = 'https://test_text_urls.matplotlib.org/'

    fig = plt.figure(figsize=(2, 1))
    fig.text(0.1, 0.1, 'test plain 123', url=f'{test_url}plain')
    fig.text(0.1, 0.4, 'test mathtext $123$', url=f'{test_url}mathtext')

    with io.BytesIO() as fd:
        fig.savefig(fd, format='pdf')

        with pikepdf.Pdf.open(fd) as pdf:
            annots = pdf.pages[0].Annots

            # Iteration over Annots must occur within the context manager,
            # otherwise it may fail depending on the pdf structure.
            for y, fragment in [('0.1', 'plain'), ('0.4', 'mathtext')]:
                annot = next(
                    (a for a in annots if a.A.URI == f'{test_url}{fragment}'),
                    None)
                assert annot is not None
                assert getattr(annot, 'QuadPoints', None) is None
                # Positions in points (72 per inch.)
                assert annot.Rect[1] == decimal.Decimal(y) * 72


def test_text_rotated_urls():
    pikepdf = pytest.importorskip('pikepdf')

    test_url = 'https://test_text_urls.matplotlib.org/'

    fig = plt.figure(figsize=(1, 1))
    fig.text(0.1, 0.1, 'N', rotation=45, url=f'{test_url}')

    with io.BytesIO() as fd:
        fig.savefig(fd, format='pdf')

        with pikepdf.Pdf.open(fd) as pdf:
            annots = pdf.pages[0].Annots

            # Iteration over Annots must occur within the context manager,
            # otherwise it may fail depending on the pdf structure.
            annot = next(
                (a for a in annots if a.A.URI == f'{test_url}'),
                None)
            assert annot is not None
            assert getattr(annot, 'QuadPoints', None) is not None
            # Positions in points (72 per inch)
            assert annot.Rect[0] == \
               annot.QuadPoints[6] - decimal.Decimal('0.00001')


@needs_usetex
def test_text_urls_tex():
    pikepdf = pytest.importorskip('pikepdf')

    test_url = 'https://test_text_urls.matplotlib.org/'

    fig = plt.figure(figsize=(2, 1))
    fig.text(0.1, 0.7, 'test tex $123$', usetex=True, url=f'{test_url}tex')

    with io.BytesIO() as fd:
        fig.savefig(fd, format='pdf')

        with pikepdf.Pdf.open(fd) as pdf:
            annots = pdf.pages[0].Annots

            # Iteration over Annots must occur within the context manager,
            # otherwise it may fail depending on the pdf structure.
            annot = next(
                (a for a in annots if a.A.URI == f'{test_url}tex'),
                None)
            assert annot is not None
            # Positions in points (72 per inch.)
            assert annot.Rect[1] == decimal.Decimal('0.7') * 72


def test_pdfpages_fspath(tmp_path):
    with PdfPages(tmp_path / 'unused.pdf') as pdf:
        pdf.savefig(plt.figure())


@image_comparison(['hatching_legend.pdf'], style='mpl20')
def test_hatching_legend(text_placeholders):
    """Test for correct hatching on patches in legend"""
    fig = plt.figure(figsize=(1, 2))

    a = Rectangle([0, 0], 0, 0, facecolor="green", hatch="XXXX")
    b = Rectangle([0, 0], 0, 0, facecolor="blue", hatch="XXXX")

    # Verify that hatches in PDFs work after empty labels. See
    # https://github.com/matplotlib/matplotlib/issues/4469
    fig.legend([a, b, a, b], ["", "", "", ""])


@image_comparison(['grayscale_alpha.pdf'], style='_classic_test')
def test_grayscale_alpha():
    """Masking images with NaN did not work for grayscale images"""
    x, y = np.ogrid[-2:2:.1, -2:2:.1]
    dd = np.exp(-(x**2 + y**2))
    dd[dd < .1] = np.nan
    fig, ax = plt.subplots()
    ax.imshow(dd, interpolation='none', cmap='gray_r')
    ax.set_xticks([])
    ax.set_yticks([])


@mpl.style.context('default')
@check_figures_equal(extensions=["pdf", "eps"])
def test_pdf_eps_savefig_when_color_is_none(fig_test, fig_ref):
    ax_test = fig_test.add_subplot()
    ax_test.set_axis_off()
    ax_test.plot(np.sin(np.linspace(-5, 5, 100)), "v", c="none")
    ax_ref = fig_ref.add_subplot()
    ax_ref.set_axis_off()


@needs_usetex
def test_failing_latex():
    """Test failing latex subprocess call"""
    plt.xlabel("$22_2_2$", usetex=True)  # This fails with "Double subscript"
    with pytest.raises(RuntimeError):
        plt.savefig(io.BytesIO(), format="pdf")


def test_empty_rasterized():
    # Check that empty figures that are rasterised save to pdf files fine
    fig, ax = plt.subplots()
    ax.plot([], [], rasterized=True)
    fig.savefig(io.BytesIO(), format="pdf")


@image_comparison(['kerning.pdf'], style='mpl20')
def test_kerning():
    fig = plt.figure()
    s = "AVAVAVAVAVAVAVAV€AAVV"
    fig.text(0, .25, s, size=5)
    fig.text(0, .75, s, size=20)


def test_glyphs_subset():
    fpath = str(_get_data_path("fonts/ttf/DejaVuSerif.ttf"))
    chars = "these should be subsetted! 1234567890"

    # non-subsetted FT2Font
    nosubfont = FT2Font(fpath)
    nosubfont.set_text(chars)
    nosubcmap = nosubfont.get_charmap()

    # subsetted FT2Font
    glyph_indices = {nosubcmap[ord(c)] for c in chars}
    with get_glyphs_subset(fm.FontPath(fpath, 0), glyph_indices) as subset:
        subfont = FT2Font(font_as_file(subset.font))
    subfont.set_text(chars)
    subcmap = subfont.get_charmap()

    # all unique chars must be available in subsetted font
    assert {*chars} == {chr(key) for key in subcmap}

    # subsetted font's charmap should have less entries
    assert len(subcmap) < len(nosubcmap)

    # since both objects are assigned same characters
    assert subfont.get_num_glyphs() == nosubfont.get_num_glyphs()


@image_comparison(["multi_font_type3.pdf"], style='mpl20')
def test_multi_font_type3():
    fonts, test_str = _gen_multi_font_text()
    plt.rc('font', family=fonts, size=16)
    plt.rc('pdf', fonttype=3)

    fig = plt.figure(figsize=(8, 6))
    fig.text(0.5, 0.5, test_str,
             horizontalalignment='center', verticalalignment='center')


@image_comparison(["multi_font_type42.pdf"], style='mpl20')
def test_multi_font_type42():
    fonts, test_str = _gen_multi_font_text()
    plt.rc('font', family=fonts, size=16)
    plt.rc('pdf', fonttype=42)

    fig = plt.figure(figsize=(8, 6))
    fig.text(0.5, 0.5, test_str,
             horizontalalignment='center', verticalalignment='center')


@image_comparison(['ttc_type3.pdf'], style='mpl20')
def test_ttc_type3():
    fp = fm.FontProperties(family=['WenQuanYi Zen Hei'])
    if Path(fm.findfont(fp)).name != 'wqy-zenhei.ttc':
        pytest.skip('Font wqy-zenhei.ttc may be missing')

    fonts = ['WenQuanYi Zen Hei', 'WenQuanYi Zen Hei Mono']
    plt.rc('font', size=16)
    plt.rc('pdf', fonttype=3)

    figs = plt.figure(figsize=(7, len(fonts) / 2)).subfigures(len(fonts))
    for font, fig in zip(fonts, figs):
        fig.text(0.5, 0.5, f'{font}: {string.ascii_uppercase}', font=font,
                 horizontalalignment='center', verticalalignment='center')


@image_comparison(['ttc_type42.pdf'], style='mpl20')
def test_ttc_type42():
    fp = fm.FontProperties(family=['WenQuanYi Zen Hei'])
    if Path(fm.findfont(fp)).name != 'wqy-zenhei.ttc':
        pytest.skip('Font wqy-zenhei.ttc may be missing')

    fonts = ['WenQuanYi Zen Hei', 'WenQuanYi Zen Hei Mono']
    plt.rc('font', size=16)
    plt.rc('pdf', fonttype=42)

    figs = plt.figure(figsize=(7, len(fonts) / 2)).subfigures(len(fonts))
    for font, fig in zip(fonts, figs):
        fig.text(0.5, 0.5, f'{font}: {string.ascii_uppercase}', font=font,
                 horizontalalignment='center', verticalalignment='center')


@pytest.mark.parametrize('family_name, file_name',
                         [("Noto Sans", "NotoSans-Regular.otf"),
                          ("FreeMono", "FreeMono.otf")])
def test_otf_font_smoke(family_name, file_name):
    # checks that there's no segfault
    fp = fm.FontProperties(family=[family_name])
    if Path(fm.findfont(fp)).name != file_name:
        pytest.skip(f"Font {family_name} may be missing")

    plt.rc('font', family=[family_name], size=27)

    fig = plt.figure()
    fig.text(0.15, 0.475, "Привет мир!")
    fig.savefig(io.BytesIO(), format="pdf")


@image_comparison(["truetype-conversion.pdf"], style='mpl20')
# mpltest.ttf does not have "l"/"p" glyphs so we get a warning when trying to
# get the font extents.
def test_truetype_conversion(recwarn):
    mpl.rcParams['pdf.fonttype'] = 3
    fig, ax = plt.subplots()
    ax.text(0, 0, "ABCDE",
            font=Path(__file__).parent / "data/mpltest.ttf", fontsize=72)
    ax.set_xticks([])
    ax.set_yticks([])


@pytest.mark.skipif(not _has_tex_package("heuristica"),
                    reason="LaTeX lacks heuristica package")
@image_comparison(["font-heuristica.pdf"], style='_classic_test')
def test_font_heuristica():
    # Heuristica uses the callothersubr operator for some glyphs
    mpl.rcParams['text.latex.preamble'] = '\n'.join((
        r'\usepackage{heuristica}',
        r'\usepackage[T1]{fontenc}',
        r'\usepackage[utf8]{inputenc}'
    ))
    fig, ax = plt.subplots()
    ax.text(0.1, 0.1, r"BHTem fi ffl 1234", usetex=True, fontsize=50)
    ax.set_xticks([])
    ax.set_yticks([])


@pytest.mark.skipif(not _has_tex_package("DejaVuSans"),
                    reason="LaTeX lacks DejaVuSans package")
@image_comparison(["font-dejavusans.pdf"], style='_classic_test')
def test_font_dejavusans():
    # DejaVuSans uses the seac operator to compose characters with diacritics
    mpl.rcParams['text.latex.preamble'] = '\n'.join((
        r'\usepackage{DejaVuSans}',
        r'\usepackage[T1]{fontenc}',
        r'\usepackage[utf8]{inputenc}'
    ))

    fig, ax = plt.subplots()
    ax.text(0.1, 0.1, r"\textsf{ñäö ABCDabcd}", usetex=True, fontsize=50)
    ax.text(0.1, 0.3, r"\textsf{fi ffl 1234}", usetex=True, fontsize=50)
    ax.set_xticks([])
    ax.set_yticks([])


@pytest.mark.skipif(not _has_tex_package("charter"),
                    reason="LaTeX lacks charter package")
@image_comparison(["font-bitstream-charter.pdf"], style='_classic_test')
def test_font_bitstream_charter():
    mpl.rcParams['text.latex.preamble'] = '\n'.join((
        r'\usepackage{charter}',
        r'\usepackage[T1]{fontenc}',
        r'\usepackage[utf8]{inputenc}'
    ))
    fig, ax = plt.subplots()
    ax.text(0.1, 0.1, r"åüš ABCDabcd", usetex=True, fontsize=50)
    ax.text(0.1, 0.3, r"fi ffl 1234", usetex=True, fontsize=50)
    ax.set_xticks([])
    ax.set_yticks([])


def test_scatter_offaxis_colored_pdf_size():
    """
    Test that off-axis scatter plots with per-point colors don't bloat PDFs.

    Regression test for issue #2488. When scatter points with per-point colors
    are completely outside the visible axes, the PDF backend should skip
    writing those markers to significantly reduce file size.
    """
    # Use John Hunter's birthday as random seed for reproducibility
    rng = np.random.default_rng(19680801)

    n_points = 1000
    x = rng.random(n_points) * 10
    y = rng.random(n_points) * 10
    c = rng.random(n_points)

    # Test 1: Scatter with per-point colors, all points OFF-AXIS
    fig1, ax1 = plt.subplots()
    ax1.scatter(x, y, c=c)
    ax1.set_xlim(20, 30)  # Move view completely away from data (x is 0-10)
    ax1.set_ylim(20, 30)  # Move view completely away from data (y is 0-10)

    buf1 = io.BytesIO()
    fig1.savefig(buf1, format='pdf')
    size_offaxis_colored = buf1.tell()
    plt.close(fig1)

    # Test 2: Empty scatter (baseline - accounts for scatter call overhead)
    fig2, ax2 = plt.subplots()
    ax2.scatter([], [])  # Empty scatter to match the axes structure
    ax2.set_xlim(20, 30)
    ax2.set_ylim(20, 30)

    buf2 = io.BytesIO()
    fig2.savefig(buf2, format='pdf')
    size_empty = buf2.tell()
    plt.close(fig2)

    # Test 3: Scatter with visible markers (should be much larger)
    fig3, ax3 = plt.subplots()
    ax3.scatter(x + 20, y + 20, c=c)  # Shift points to be visible
    ax3.set_xlim(20, 30)
    ax3.set_ylim(20, 30)

    buf3 = io.BytesIO()
    fig3.savefig(buf3, format='pdf')
    size_visible = buf3.tell()
    plt.close(fig3)

    # The off-axis colored scatter should be close to empty size.
    # Since the axes are identical, the difference should be minimal
    # (just the scatter collection setup, no actual marker data).
    # Use a tight tolerance since axes output is identical.
    assert size_offaxis_colored < size_empty + 5_000, (
        f"Off-axis colored scatter PDF ({size_offaxis_colored} bytes) is too large. "
        f"Expected close to empty scatter size ({size_empty} bytes). "
        f"Markers may not be properly skipped."
    )

    # The visible scatter should be significantly larger than both empty and
    # off-axis, demonstrating the optimization is working.
    assert size_visible > size_empty + 15_000, (
        f"Visible scatter PDF ({size_visible} bytes) should be much larger "
        f"than empty ({size_empty} bytes) to validate the test."
    )
    assert size_visible > size_offaxis_colored + 15_000, (
        f"Visible scatter PDF ({size_visible} bytes) should be much larger "
        f"than off-axis ({size_offaxis_colored} bytes) to validate optimization."
    )


@check_figures_equal(extensions=["pdf"])
def test_scatter_offaxis_colored_visual(fig_test, fig_ref):
    """
    Test that on-axis scatter with per-point colors still renders correctly.

    Ensures the optimization for off-axis markers doesn't break normal
    scatter rendering.
    """
    rng = np.random.default_rng(19680801)

    n_points = 100
    x = rng.random(n_points) * 5
    y = rng.random(n_points) * 5
    c = rng.random(n_points)

    # Test figure: scatter with clipping optimization
    ax_test = fig_test.subplots()
    ax_test.scatter(x, y, c=c, s=50)
    ax_test.set_xlim(0, 10)
    ax_test.set_ylim(0, 10)

    # Reference figure: should look identical
    ax_ref = fig_ref.subplots()
    ax_ref.scatter(x, y, c=c, s=50)
    ax_ref.set_xlim(0, 10)
    ax_ref.set_ylim(0, 10)


@check_figures_equal(extensions=["pdf"])
def test_scatter_mixed_onoff_axis(fig_test, fig_ref):
    """
    Test scatter with some points on-axis and some off-axis.

    Ensures the optimization correctly handles the common case where only
    some markers are outside the visible area.
    """
    rng = np.random.default_rng(19680801)

    # Create points: half on-axis (0-5), half off-axis (15-20)
    n_points = 50
    x_on = rng.random(n_points) * 5
    y_on = rng.random(n_points) * 5
    x_off = rng.random(n_points) * 5 + 15
    y_off = rng.random(n_points) * 5 + 15

    x = np.concatenate([x_on, x_off])
    y = np.concatenate([y_on, y_off])
    c = rng.random(2 * n_points)

    # Test figure: scatter with mixed points
    ax_test = fig_test.subplots()
    ax_test.scatter(x, y, c=c, s=50)
    ax_test.set_xlim(0, 10)
    ax_test.set_ylim(0, 10)

    # Reference figure: only the on-axis points should be visible
    ax_ref = fig_ref.subplots()
    ax_ref.scatter(x_on, y_on, c=c[:n_points], s=50)
    ax_ref.set_xlim(0, 10)
    ax_ref.set_ylim(0, 10)


@check_figures_equal(extensions=["pdf"])
def test_scatter_large_markers_partial_clip(fig_test, fig_ref):
    """
    Test that large markers are rendered when partially visible.

    Addresses reviewer concern: markers with centers outside the canvas but
    with edges extending into the visible area should still be rendered.
    """
    # Create markers just outside the visible area
    # Canvas is 0-10, markers at x=-0.5 and x=10.5
    x = np.array([-0.5, 10.5, 5])  # left edge, right edge, center
    y = np.array([5, 5, -0.5])  # center, center, bottom edge
    c = np.array([0.2, 0.5, 0.8])

    # Test figure: large markers (s=500 ≈ 11 points radius)
    # Centers are outside, but marker edges extend into visible area
    ax_test = fig_test.subplots()
    ax_test.scatter(x, y, c=c, s=500)
    ax_test.set_xlim(0, 10)
    ax_test.set_ylim(0, 10)

    # Reference figure: same plot (should render identically)
    ax_ref = fig_ref.subplots()
    ax_ref.scatter(x, y, c=c, s=500)
    ax_ref.set_xlim(0, 10)
    ax_ref.set_ylim(0, 10)


@check_figures_equal(extensions=["pdf"])
def test_scatter_logscale(fig_test, fig_ref):
    """
    Test scatter optimization with logarithmic scales.

    Ensures bounds checking works correctly in log-transformed coordinates.
    """
    rng = np.random.default_rng(19680801)

    # Create points across several orders of magnitude
    n_points = 50
    x = 10 ** (rng.random(n_points) * 4)  # 1 to 10000
    y = 10 ** (rng.random(n_points) * 4)
    c = rng.random(n_points)

    # Test figure: log scale with points mostly outside view
    ax_test = fig_test.subplots()
    ax_test.scatter(x, y, c=c, s=50)
    ax_test.set_xscale('log')
    ax_test.set_yscale('log')
    ax_test.set_xlim(100, 1000)  # Only show middle range
    ax_test.set_ylim(100, 1000)

    # Reference figure: should render identically
    ax_ref = fig_ref.subplots()
    ax_ref.scatter(x, y, c=c, s=50)
    ax_ref.set_xscale('log')
    ax_ref.set_yscale('log')
    ax_ref.set_xlim(100, 1000)
    ax_ref.set_ylim(100, 1000)


@check_figures_equal(extensions=["pdf"])
def test_scatter_polar(fig_test, fig_ref):
    """
    Test scatter optimization with polar coordinates.

    Ensures bounds checking works correctly in polar projections.
    """
    rng = np.random.default_rng(19680801)

    n_points = 50
    theta = rng.random(n_points) * 2 * np.pi
    r = rng.random(n_points) * 3
    c = rng.random(n_points)

    # Test figure: polar projection
    ax_test = fig_test.subplots(subplot_kw={'projection': 'polar'})
    ax_test.scatter(theta, r, c=c, s=50)
    ax_test.set_ylim(0, 2)  # Limit radial range

    # Reference figure: should render identically
    ax_ref = fig_ref.subplots(subplot_kw={'projection': 'polar'})
    ax_ref.scatter(theta, r, c=c, s=50)
    ax_ref.set_ylim(0, 2)
