from pathlib import Path
import io

import pytest

from matplotlib import ft2font
from matplotlib.testing.decorators import check_figures_equal
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt


def test_fallback_errors():
    file_name = fm.findfont('DejaVu Sans')

    with pytest.raises(TypeError, match="Fallback list must be a list"):
        # failing to be a list will fail before the 0
        ft2font.FT2Font(file_name, _fallback_list=(0,))  # type: ignore[arg-type]

    with pytest.raises(
            TypeError, match="Fallback fonts must be FT2Font objects."
    ):
        ft2font.FT2Font(file_name, _fallback_list=[0])  # type: ignore[list-item]


def test_ft2font_positive_hinting_factor():
    file_name = fm.findfont('DejaVu Sans')
    with pytest.raises(
            ValueError, match="hinting_factor must be greater than 0"
    ):
        ft2font.FT2Font(file_name, 0)


@pytest.mark.parametrize('family_name, file_name',
                          [("WenQuanYi Zen Hei",  "wqy-zenhei.ttc"),
                           ("Noto Sans CJK JP", "NotoSansCJK.ttc"),
                           ("Noto Sans TC", "NotoSansTC-Regular.otf")]
                         )
def test_fallback_smoke(family_name, file_name):
    fp = fm.FontProperties(family=[family_name])
    if Path(fm.findfont(fp)).name != file_name:
        pytest.skip(f"Font {family_name} ({file_name}) is missing")
    plt.rcParams['font.size'] = 20
    fig = plt.figure(figsize=(4.75, 1.85))
    fig.text(0.05, 0.45, "There are 几个汉字 in between!",
             family=['DejaVu Sans', family_name])
    fig.text(0.05, 0.85, "There are 几个汉字 in between!",
             family=[family_name])

    # TODO enable fallback for other backends!
    for fmt in ['png', 'raw']:  # ["svg", "pdf", "ps"]:
        fig.savefig(io.BytesIO(), format=fmt)


@pytest.mark.parametrize('family_name, file_name',
                         [("WenQuanYi Zen Hei",  "wqy-zenhei"),
                          ("Noto Sans CJK JP", "NotoSansCJK"),
                          ("Noto Sans TC", "NotoSansTC-Regular.otf")]
                         )
@check_figures_equal(extensions=["png", "pdf", "eps", "svg"])
def test_font_fallback_chinese(fig_test, fig_ref, family_name, file_name):
    fp = fm.FontProperties(family=[family_name])
    if file_name not in Path(fm.findfont(fp)).name:
        pytest.skip(f"Font {family_name} ({file_name}) is missing")

    text = ["There are", "几个汉字", "in between!"]

    plt.rcParams["font.size"] = 20
    test_fonts = [["DejaVu Sans", family_name]] * 3
    ref_fonts = [["DejaVu Sans"], [family_name], ["DejaVu Sans"]]

    for j, (txt, test_font, ref_font) in enumerate(
            zip(text, test_fonts, ref_fonts)
    ):
        fig_ref.text(0.05, .85 - 0.15*j, txt, family=ref_font)
        fig_test.text(0.05, .85 - 0.15*j, txt, family=test_font)


@pytest.mark.parametrize("font_list",
                          [['DejaVu Serif', 'DejaVu Sans'],
                           ['DejaVu Sans Mono']],
                         ids=["two fonts", "one font"])
def test_fallback_missing(recwarn, font_list):
    fig = plt.figure()
    fig.text(.5, .5, "Hello 🙃 World!", family=font_list)
    fig.canvas.draw()
    assert all(isinstance(warn.message, UserWarning) for warn in recwarn)
    # not sure order is guaranteed on the font listing so
    assert recwarn[0].message.args[0].startswith(
           "Glyph 128579 (\\N{UPSIDE-DOWN FACE}) missing from font(s)")
    assert all([font in recwarn[0].message.args[0] for font in font_list])


@pytest.mark.parametrize(
    "family_name, file_name",
    [
        ("WenQuanYi Zen Hei", "wqy-zenhei"),
        ("Noto Sans CJK JP", "NotoSansCJK"),
        ("Noto Sans TC", "NotoSansTC-Regular.otf")
    ],
)
def test__get_fontmap(family_name, file_name):
    fp = fm.FontProperties(family=[family_name])
    found_file_name = Path(fm.findfont(fp)).name
    if file_name not in found_file_name:
        pytest.skip(f"Font {family_name} ({file_name}) is missing")

    text = "There are 几个汉字 in between!"
    ft = fm.get_font(
        fm.fontManager._find_fonts_by_props(
            fm.FontProperties(family=["DejaVu Sans", family_name])
        )
    )
    fontmap = ft._get_fontmap(text)
    for char, font in fontmap.items():
        if ord(char) > 127:
            assert Path(font.fname).name == found_file_name
        else:
            assert Path(font.fname).name == "DejaVuSans.ttf"
