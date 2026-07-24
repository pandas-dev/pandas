from io import BytesIO
import gc
import multiprocessing
import os
from pathlib import Path
from PIL import Image
import shutil
import sys
import warnings

import pytest

from unittest.mock import MagicMock, patch

import matplotlib as mpl
import matplotlib.font_manager as fm_mod
from matplotlib.font_manager import (
    findfont, findSystemFonts, FontEntry, FontPath, FontProperties, fontManager,
    json_dump, json_load, get_font, is_opentype_cff_font,
    MSUserFontDirectories, ttfFontProperty, _get_font_alt_names,
    _get_fontconfig_fonts, _normalize_weight)
from matplotlib import cbook, ft2font, pyplot as plt, rc_context, figure as mfigure
from matplotlib.testing import subprocess_run_helper, subprocess_run_for_testing


has_fclist = sys.platform != 'emscripten' and shutil.which('fc-list') is not None


def test_font_path():
    fp = FontPath('foo', 123)
    fp2 = FontPath('foo', 321)
    assert str(fp) == 'foo'
    assert repr(fp) == "FontPath('foo', 123)"
    assert fp.path == 'foo'
    assert fp.face_index == 123
    # Should be immutable.
    with pytest.raises(AttributeError, match='has no setter'):
        fp.path = 'bar'
    with pytest.raises(AttributeError, match='has no setter'):
        fp.face_index = 321
    # Should be comparable with str and itself.
    assert fp == 'foo'
    assert fp == FontPath('foo', 123)
    assert fp <= fp
    assert fp >= fp
    assert fp != fp2
    assert fp < fp2
    assert fp <= fp2
    assert fp2 > fp
    assert fp2 >= fp
    # Should be hashable, but not the same as str.
    d = {fp: 1, 'bar': 2}
    assert fp in d
    assert d[fp] == 1
    assert d[FontPath('foo', 123)] == 1
    assert fp2 not in d
    assert 'foo' not in d
    assert FontPath('bar', 0) not in d


def test_font_priority():
    with rc_context(rc={
            'font.sans-serif':
            ['cmmi10', 'Bitstream Vera Sans']}):
        fontfile = findfont(FontProperties(family=["sans-serif"]))
    assert Path(fontfile).name == 'cmmi10.ttf'

    # Smoketest get_charmap, which isn't used internally anymore
    font = get_font(fontfile)
    cmap = font.get_charmap()
    assert len(cmap) == 131
    assert cmap[8729] == 30


def test_score_weight():
    assert 0 == fontManager.score_weight("regular", "regular")
    assert 0 == fontManager.score_weight("bold", "bold")
    assert (0 < fontManager.score_weight(400, 400) <
            fontManager.score_weight("normal", "bold"))
    assert (0 < fontManager.score_weight("normal", "regular") <
            fontManager.score_weight("normal", "bold"))
    assert (fontManager.score_weight("normal", "regular") ==
            fontManager.score_weight(400, 400))


def test_json_serialization(tmp_path):
    # Can't open a NamedTemporaryFile twice on Windows, so use a temporary
    # directory instead.
    json_dump(fontManager, tmp_path / "fontlist.json")
    copy = json_load(tmp_path / "fontlist.json")
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', 'findfont: Font family.*not found')
        for prop in ({'family': 'STIXGeneral'},
                     {'family': 'Bitstream Vera Sans', 'weight': 700},
                     {'family': 'no such font family'}):
            fp = FontProperties(**prop)
            assert (fontManager.findfont(fp, rebuild_if_missing=False) ==
                    copy.findfont(fp, rebuild_if_missing=False))


def test_otf():
    fname = '/usr/share/fonts/opentype/freefont/FreeMono.otf'
    if Path(fname).exists():
        with pytest.warns(mpl.MatplotlibDeprecationWarning):
            assert is_opentype_cff_font(fname)
    for f in fontManager.ttflist:
        if 'otf' in f.fname:
            with open(f.fname, 'rb') as fd:
                res = fd.read(4) == b'OTTO'
            with pytest.warns(mpl.MatplotlibDeprecationWarning):
                assert res == is_opentype_cff_font(f.fname)


@pytest.mark.skipif(sys.platform == "win32" or not has_fclist,
                    reason='no fontconfig installed')
def test_get_fontconfig_fonts():
    assert len(_get_fontconfig_fonts()) > 1


def test_utf16m_sfnt():
    try:
        # seguisbi = Microsoft Segoe UI Semibold
        entry = next(entry for entry in fontManager.ttflist
                     if Path(entry.fname).name == "seguisbi.ttf")
    except StopIteration:
        pytest.skip("Couldn't find seguisbi.ttf font to test against.")
    else:
        # Check that we successfully read "semibold" from the font's sfnt table
        # and set its weight accordingly.
        assert entry.weight == 600


def test_find_ttc():
    fp = FontProperties(family=["WenQuanYi Zen Hei"])
    fontpath = findfont(fp)
    if Path(fontpath).name != "wqy-zenhei.ttc":
        pytest.skip("Font wqy-zenhei.ttc may be missing")
    # All fonts from this collection should have loaded as well.
    for name in ["WenQuanYi Zen Hei Mono", "WenQuanYi Zen Hei Sharp"]:
        subfontpath = findfont(FontProperties(family=[name]), fallback_to_default=False)
        assert subfontpath.path == fontpath.path
        assert subfontpath.face_index != fontpath.face_index
        subfont = get_font(subfontpath)
        assert subfont.fname == subfontpath.path
        assert subfont.face_index == subfontpath.face_index
    fig, ax = plt.subplots()
    ax.text(.5, .5, "\N{KANGXI RADICAL DRAGON}", fontproperties=fp)
    for fmt in ["raw", "svg", "pdf", "ps"]:
        fig.savefig(BytesIO(), format=fmt)


def test_find_noto():
    fp = FontProperties(family=["Noto Sans CJK SC", "Noto Sans CJK JP"])
    name = Path(findfont(fp)).name
    if name not in ("NotoSansCJKsc-Regular.otf", "NotoSansCJK-Regular.ttc"):
        pytest.skip(f"Noto Sans CJK SC font may be missing (found {name})")

    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, 'Hello, 你好', fontproperties=fp)
    for fmt in ["raw", "svg", "pdf", "ps"]:
        fig.savefig(BytesIO(), format=fmt)


def test_find_valid():
    class PathLikeClass:
        def __init__(self, filename):
            self.filename = filename

        def __fspath__(self):
            return self.filename

    file_str = findfont('DejaVu Sans')
    file_bytes = os.fsencode(file_str)

    font = get_font(file_str)
    assert font.fname == file_str
    font = get_font(file_bytes)
    assert font.fname == file_bytes
    font = get_font(PathLikeClass(file_str))
    assert font.fname == file_str
    font = get_font(PathLikeClass(file_bytes))
    assert font.fname == file_bytes
    font = get_font(FontPath(file_str, 0))
    assert font.fname == file_str

    # Note, fallbacks are not currently accessible.
    font = get_font([file_str, file_bytes,
                     PathLikeClass(file_str), PathLikeClass(file_bytes)])
    assert font.fname == file_str


def test_find_invalid(tmp_path):

    with pytest.raises(FileNotFoundError):
        get_font(tmp_path / 'non-existent-font-name.ttf')

    with pytest.raises(FileNotFoundError):
        get_font(str(tmp_path / 'non-existent-font-name.ttf'))

    with pytest.raises(FileNotFoundError):
        get_font(bytes(tmp_path / 'non-existent-font-name.ttf'))


@pytest.mark.skipif(sys.platform != 'linux' or not has_fclist,
                    reason='only Linux with fontconfig installed')
def test_user_fonts_linux(tmpdir, monkeypatch):
    font_test_file = 'mpltest.ttf'

    # Precondition: the test font should not be available
    fonts = findSystemFonts()
    if any(font_test_file in font for font in fonts):
        pytest.skip(f'{font_test_file} already exists in system fonts')

    # Prepare a temporary user font directory
    user_fonts_dir = tmpdir.join('fonts')
    user_fonts_dir.ensure(dir=True)
    shutil.copyfile(Path(__file__).parent / 'data' / font_test_file,
                    user_fonts_dir.join(font_test_file))

    with monkeypatch.context() as m:
        m.setenv('XDG_DATA_HOME', str(tmpdir))
        _get_fontconfig_fonts.cache_clear()
        # Now, the font should be available
        fonts = findSystemFonts()
        assert any(font_test_file in font for font in fonts)

    # Make sure the temporary directory is no longer cached.
    _get_fontconfig_fonts.cache_clear()


def test_addfont_as_path(monkeypatch):
    """Smoke test that addfont() accepts pathlib.Path."""
    font_test_file = 'mpltest.ttf'
    path = Path(__file__).parent / 'data' / font_test_file
    try:
        fontManager.addfont(path)
        assert fontManager.findfont('mpltest:weight=500') == FontPath(path, 0)
        with monkeypatch.context() as m, pytest.raises(ValueError):
            m.setenv('MPL_IGNORE_SYSTEM_FONTS', 'true')  # Can only find internal fonts.
            fontManager.findfont('mpltest:weight=500', fallback_to_default=False)
        added, = (font for font in fontManager.ttflist
                  if font.fname.endswith(font_test_file))
        fontManager.ttflist.remove(added)
    finally:
        to_remove = [font for font in fontManager.ttflist
                     if font.fname.endswith(font_test_file)]
        for font in to_remove:
            fontManager.ttflist.remove(font)


@pytest.mark.skipif(sys.platform != 'win32', reason='Windows only')
def test_user_fonts_win32():
    if not (os.environ.get('APPVEYOR') or os.environ.get('TF_BUILD')):
        pytest.xfail("This test should only run on CI (appveyor or azure) "
                     "as the developer's font directory should remain "
                     "unchanged.")
    pytest.xfail("We need to update the registry for this test to work")
    font_test_file = 'mpltest.ttf'

    # Precondition: the test font should not be available
    fonts = findSystemFonts()
    if any(font_test_file in font for font in fonts):
        pytest.skip(f'{font_test_file} already exists in system fonts')

    user_fonts_dir = MSUserFontDirectories[0]

    # Make sure that the user font directory exists (this is probably not the
    # case on Windows versions < 1809)
    os.makedirs(user_fonts_dir)

    # Copy the test font to the user font directory
    shutil.copy(Path(__file__).parent / 'data' / font_test_file, user_fonts_dir)

    # Now, the font should be available
    fonts = findSystemFonts()
    assert any(font_test_file in font for font in fonts)


def _model_handler(_):
    fig, ax = plt.subplots()
    fig.savefig(BytesIO(), format="pdf")
    plt.close()


@pytest.mark.skipif(sys.platform == 'emscripten',
                    reason='emscripten does not support subprocesses')
@pytest.mark.skipif(not hasattr(os, "register_at_fork"),
                    reason="Cannot register at_fork handlers")
# Python 3.15+ raises DeprecationWarning for fork in multi-threaded process
@pytest.mark.filterwarnings("ignore:.*multi-threaded.*fork.*:DeprecationWarning")
@pytest.mark.filterwarnings("ignore:.*multi-threaded.*fork.*:RuntimeWarning")
def test_fork():
    _model_handler(0)  # Make sure the font cache is filled.
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(processes=2) as pool:
        pool.map(_model_handler, range(2))


def test_missing_family(caplog):
    plt.rcParams["font.sans-serif"] = ["this-font-does-not-exist"]
    with caplog.at_level("WARNING"):
        findfont("sans")
    assert [rec.getMessage() for rec in caplog.records] == [
        "findfont: Font family ['sans'] not found. "
        "Falling back to DejaVu Sans.",
        "findfont: Generic family 'sans' not found because none of the "
        "following families were found: this-font-does-not-exist",
    ]


def _test_threading():
    import threading
    from matplotlib.ft2font import LoadFlags
    import matplotlib.font_manager as fm

    def loud_excepthook(args):
        raise RuntimeError("error in thread!")

    threading.excepthook = loud_excepthook

    N = 10
    b = threading.Barrier(N)

    def bad_idea(n):
        b.wait(timeout=5)
        for j in range(100):
            font = fm.get_font(fm.findfont("DejaVu Sans"))
            font.set_text(str(n), 0.0, flags=LoadFlags.NO_HINTING)

    threads = [
        threading.Thread(target=bad_idea, name=f"bad_thread_{j}", args=(j,))
        for j in range(N)
    ]

    for t in threads:
        t.start()

    for t in threads:
        t.join(timeout=9)
        if t.is_alive():
            raise RuntimeError("thread failed to join")


def test_fontcache_thread_safe():
    pytest.importorskip('threading')

    subprocess_run_helper(_test_threading, timeout=10)


def test_lockfilefailure(tmp_path):
    # The logic here:
    # 1. get a temp directory from pytest
    # 2. import matplotlib which makes sure it exists
    # 3. get the cache dir (where we check it is writable)
    # 4. make it not writable
    # 5. try to write into it via font manager
    proc = subprocess_run_for_testing(
        [
            sys.executable,
            "-c",
            "import matplotlib;"
            "import os;"
            "p = matplotlib.get_cachedir();"
            "os.chmod(p, 0o555);"
            "import matplotlib.font_manager;"
        ],
        env={**os.environ, 'MPLCONFIGDIR': str(tmp_path)},
        check=True
    )


def test_fontentry_dataclass():
    fontent = FontEntry(name='font-name')

    png = fontent._repr_png_()
    img = Image.open(BytesIO(png))
    assert img.width > 0
    assert img.height > 0

    html = fontent._repr_html_()
    assert html.startswith("<img src=\"data:image/png;base64")


def test_fontentry_dataclass_invalid_path():
    with pytest.raises(FileNotFoundError):
        fontent = FontEntry(fname='/random', name='font-name')
        fontent._repr_html_()


@pytest.mark.skipif(sys.platform == 'win32', reason='Linux or OS only')
def test_get_font_names():
    paths_mpl = [cbook._get_data_path('fonts', subdir) for subdir in ['ttf']]
    fonts_mpl = findSystemFonts(paths_mpl, fontext='ttf')
    fonts_system = findSystemFonts(fontext='ttf')
    ttf_fonts = set()
    for path in fonts_mpl + fonts_system:
        try:
            font = ft2font.FT2Font(path)
            prop = ttfFontProperty(font)
            ttf_fonts.add(prop.name)
            for face_index in range(1, font.num_faces):
                font = ft2font.FT2Font(path, face_index=face_index)
                prop = ttfFontProperty(font)
                ttf_fonts.add(prop.name)
        except Exception:
            pass
    # fontManager may contain additional entries for alternative family names
    # (e.g. typographic family, platform-specific Name ID 1) registered by
    # addfont(), so primary names must be a subset of the manager's names.
    assert ttf_fonts <= set(fontManager.get_font_names())


def test_addfont_alternative_names(tmp_path):
    """
    Fonts that advertise different family names across platforms or name IDs
    should be registered under all of those names so users can address the font
    by any of them.

    Two real-world patterns are covered:

    - **MS platform ID 1 differs from Mac platform ID 1** (e.g. Ubuntu Light):
      FreeType returns the Mac ID 1 value as ``family_name``; the MS ID 1
      value ("Ubuntu Light") is an equally valid name that users expect to work.
    - **Name ID 16 (Typographic Family) differs from ID 1** (older fonts):
      some fonts store a broader family name in ID 16.
    """
    mac_key = (1, 0, 0)
    ms_key = (3, 1, 0x0409)

    # Case 1: MS ID1 differs from Mac ID1 (Ubuntu Light pattern)
    # Mac ID1="Test Family" → FreeType family_name (primary)
    # MS  ID1="Test Family Light" → alternate name users expect to work
    ubuntu_style_sfnt = {
        (*mac_key, 1): "Test Family".encode("latin-1"),
        (*ms_key,  1): "Test Family Light".encode("utf-16-be"),
        (*mac_key, 2): "Light".encode("latin-1"),
        (*ms_key,  2): "Regular".encode("utf-16-be"),
    }
    fake_font = MagicMock()
    fake_font.get_sfnt.return_value = ubuntu_style_sfnt

    assert _get_font_alt_names(fake_font, "Test Family") == [("Test Family Light", 400)]
    assert _get_font_alt_names(fake_font, "Test Family Light") == [
        ("Test Family", 300)]

    # Case 2: ID 16 differs from ID 1 (older typographic-family pattern)
    # ID 17 (typographic subfamily) is absent → defaults to weight 400
    id16_sfnt = {
        (*mac_key, 1):  "Test Family".encode("latin-1"),
        (*ms_key,  1):  "Test Family".encode("utf-16-be"),
        (*ms_key,  16): "Test Family Light".encode("utf-16-be"),
    }
    fake_font_id16 = MagicMock()
    fake_font_id16.get_sfnt.return_value = id16_sfnt

    assert _get_font_alt_names(
        fake_font_id16, "Test Family"
    ) == [("Test Family Light", 400)]

    # Case 3: all entries agree → no alternates
    same_sfnt = {
        (*mac_key, 1): "Test Family".encode("latin-1"),
        (*ms_key,  1): "Test Family".encode("utf-16-be"),
    }
    fake_font_same = MagicMock()
    fake_font_same.get_sfnt.return_value = same_sfnt
    assert _get_font_alt_names(fake_font_same, "Test Family") == []

    # Case 4: get_sfnt() raises ValueError (e.g. non-SFNT font) → empty list
    fake_font_no_sfnt = MagicMock()
    fake_font_no_sfnt.get_sfnt.side_effect = ValueError
    assert _get_font_alt_names(fake_font_no_sfnt, "Test Family") == []

    fake_path = str(tmp_path / "fake.ttf")
    primary_entry = FontEntry(fname=fake_path, name="Test Family",
                              style="normal", variant="normal",
                              weight=300, stretch="normal", size="scalable")

    with patch("matplotlib.font_manager.ft2font.FT2Font",
               return_value=fake_font), \
         patch("matplotlib.font_manager.ttfFontProperty",
               return_value=primary_entry):
        fm_instance = fm_mod.FontManager.__new__(fm_mod.FontManager)
        fm_instance.ttflist = []
        fm_instance.afmlist = []
        fm_instance._findfont_cached = MagicMock()
        fm_instance._findfont_cached.cache_clear = MagicMock()
        fm_instance.addfont(fake_path)

    names = [e.name for e in fm_instance.ttflist]
    assert names == ["Test Family", "Test Family Light"]
    alt_entry = fm_instance.ttflist[1]
    assert alt_entry.weight == 400
    assert alt_entry.style == primary_entry.style
    assert alt_entry.fname == primary_entry.fname


@pytest.mark.parametrize("subfam,expected", [
    ("Thin",        100),
    ("ExtraLight",  200),
    ("UltraLight",  200),
    ("DemiLight",   350),
    ("SemiLight",   350),
    ("Light",       300),
    ("Book",        380),
    ("Regular",     400),
    ("Normal",      400),
    ("Medium",      500),
    ("DemiBold",    600),
    ("Demi",        600),
    ("SemiBold",    600),
    ("ExtraBold",   800),
    ("SuperBold",   800),
    ("UltraBold",   800),
    ("Bold",        700),
    ("UltraBlack", 1000),
    ("SuperBlack", 1000),
    ("ExtraBlack", 1000),
    ("Ultra",      1000),
    ("Black",       900),
    ("Heavy",       900),
    ("",            400),  # fallback: unrecognised → regular
])
def test_alt_name_weight_from_subfamily(subfam, expected):
    """_get_font_alt_names derives weight from the paired subfamily string."""
    ms_key = (3, 1, 0x0409)
    fake_font = MagicMock()
    fake_font.get_sfnt.return_value = {
        (*ms_key, 1): "Family Alt".encode("utf-16-be"),
        (*ms_key, 2): subfam.encode("utf-16-be"),
    }
    result = _get_font_alt_names(fake_font, "Family")
    assert result == [("Family Alt", expected)]


def test_donot_cache_tracebacks():

    class SomeObject:
        pass

    def inner():
        x = SomeObject()
        fig = mfigure.Figure()
        ax = fig.subplots()
        fig.text(.5, .5, 'aardvark', family='doesnotexist')
        with BytesIO() as out:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                fig.savefig(out, format='raw')

    inner()

    for obj in gc.get_objects():
        if isinstance(obj, SomeObject):
            pytest.fail("object from inner stack still alive")


def test_fontproperties_init_deprecation():
    """
    Test the deprecated API of FontProperties.__init__.

    The deprecation does not change behavior, it only adds a deprecation warning
    via a decorator. Therefore, the purpose of this test is limited to check
    which calls do and do not issue deprecation warnings. Behavior is still
    tested via the existing regular tests.
    """
    with pytest.warns(mpl.MatplotlibDeprecationWarning):
        # multiple positional arguments
        FontProperties("Times", "italic")

    with pytest.warns(mpl.MatplotlibDeprecationWarning):
        # Mixed positional and keyword arguments
        FontProperties("Times", size=10)

    with pytest.warns(mpl.MatplotlibDeprecationWarning):
        # passing a family list positionally
        FontProperties(["Times"])

    # still accepted:
    FontProperties(family="Times", style="italic")
    FontProperties(family="Times")
    FontProperties("Times")  # works as pattern and family
    FontProperties("serif-24:style=oblique:weight=bold")  # pattern

    # also still accepted:
    # passing as pattern via family kwarg was not covered by the docs but
    # historically worked. This is left unchanged for now.
    # AFAICT, we cannot detect this: We can determine whether a string
    # works as pattern, but that doesn't help, because there are strings
    # that are both pattern and family. We would need to identify, whether
    # a string is *not* a valid family.
    # Since this case is not covered by docs, I've refrained from jumping
    # extra hoops to detect this possible API misuse.
    FontProperties(family="serif-24:style=oblique:weight=bold")


def test_normalize_weights():
    assert _normalize_weight(300) == 300  # passthrough
    assert _normalize_weight('ultralight') == 100
    assert _normalize_weight('light') == 200
    assert _normalize_weight('normal') == 400
    assert _normalize_weight('regular') == 400
    assert _normalize_weight('book') == 400
    assert _normalize_weight('medium') == 500
    assert _normalize_weight('roman') == 500
    assert _normalize_weight('semibold') == 600
    assert _normalize_weight('demibold') == 600
    assert _normalize_weight('demi') == 600
    assert _normalize_weight('bold') == 700
    assert _normalize_weight('heavy') == 800
    assert _normalize_weight('extra bold') == 800
    assert _normalize_weight('black') == 900
    with pytest.raises(KeyError):
        _normalize_weight('invalid')


def test_font_match_warning(caplog):
    findfont(FontProperties(family=["DejaVu Sans"], weight=750))
    logs = [rec.message for rec in caplog.records]
    assert 'findfont: Failed to find font weight 750, now using 700.' in logs


def test_mutable_fontproperty_cache_invalidation():
    fp = FontProperties()
    assert findfont(fp).endswith("DejaVuSans.ttf")
    fp.set_weight("bold")
    assert findfont(fp).endswith("DejaVuSans-Bold.ttf")


def test_fontproperty_default_cache_invalidation():
    mpl.rcParams["font.weight"] = "normal"
    assert findfont("DejaVu Sans").endswith("DejaVuSans.ttf")
    mpl.rcParams["font.weight"] = "bold"
    assert findfont("DejaVu Sans").endswith("DejaVuSans-Bold.ttf")
