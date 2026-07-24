from __future__ import annotations

__lazy_modules__ = {f"{__spec__.parent}.termsize"}

import sys
import typing

from plumbum import cli, colors

from .termsize import get_terminal_size

if typing.TYPE_CHECKING:
    import PIL.Image


class Image:
    __slots__ = ["char_ratio", "size"]

    def __init__(self, size: tuple[int, int] | None = None, char_ratio: float = 2.45):
        self.size = size
        self.char_ratio = char_ratio

    def best_aspect(
        self, orig: tuple[int, int], term: tuple[int, int]
    ) -> tuple[int, int]:
        """Select a best possible size matching the original aspect ratio.
        Size is width, height.
        The char_ratio option gives the height of each char with respect
        to its width, zero for no effect."""

        if not self.char_ratio:  # Don't use if char ratio is 0
            return term

        # Cells are char_ratio times taller than they are wide, so a row count
        # spans char_ratio "width-units" of height. To preserve the image's
        # width/height aspect, the width in cells should equal
        # height_in_cells * char_ratio * (orig_width / orig_height).
        cell_ratio = orig[0] / orig[1] * self.char_ratio

        # First try filling the full terminal height; if that overflows the
        # available width, fall back to filling the full width instead.
        width = int(term[1] * cell_ratio)
        if width <= term[0]:
            return width, term[1]

        return term[0], int(term[0] / cell_ratio)

    def show(self, filename: str, double: bool = False) -> None:
        """Display an image on the command line. Can select a size or show in double resolution."""

        import PIL.Image

        return (
            self.show_pil_double(PIL.Image.open(filename))
            if double
            else self.show_pil(PIL.Image.open(filename))
        )

    def _init_size(self, im: PIL.Image.Image) -> tuple[int, int]:
        """Return the expected image size"""
        if self.size is None:
            term_size = get_terminal_size()
            return self.best_aspect(im.size, term_size)

        return self.size

    def show_pil(self, im: PIL.Image.Image) -> None:
        "Standard show routine"
        size = self._init_size(im)
        new_im = im.resize(size).convert("RGB")

        for y in range(size[1]):
            for x in range(size[0] - 1):
                pix = new_im.getpixel((x, y))
                sys.stdout.write(
                    colors.bg.rgb(*pix) + " "  # type: ignore[misc]
                )  # '\u2588'
            sys.stdout.write(colors.reset + " \n")
        sys.stdout.write(colors.reset + "\n")
        sys.stdout.flush()

    def show_pil_double(self, im: PIL.Image.Image) -> None:
        "Show double resolution on some fonts"

        size = self._init_size(im)
        size = (size[0], size[1] * 2)
        new_im = im.resize(size).convert("RGB")

        for y in range(size[1] // 2):
            for x in range(size[0] - 1):
                pix = new_im.getpixel((x, y * 2))
                pixl = new_im.getpixel((x, y * 2 + 1))
                sys.stdout.write(
                    (colors.bg.rgb(*pixl) & colors.fg.rgb(*pix)) + "\u2580"  # type: ignore[misc]
                )
            sys.stdout.write(colors.reset + " \n")
        sys.stdout.write(colors.reset + "\n")
        sys.stdout.flush()


class ShowImageApp(cli.Application):
    "Display an image on the terminal"

    double = cli.Flag(
        ["-d", "--double"], help="Double resolution (looks good only with some fonts)"
    )

    @cli.switch(["-c", "--colors"], cli.Range(1, 4), help="Level of color, 1-4")
    def colors_set(self, n: int) -> None:  # pylint: disable=no-self-use
        colors.use_color = n

    size = cli.SwitchAttr(["-s", "--size"], help="Size, should be in the form 100x150")

    ratio = cli.SwitchAttr(
        ["--ratio"], float, default=2.45, help="Aspect ratio of the font"
    )

    @cli.positional(cli.ExistingFile)
    def main(self, filename: str) -> None:
        size = tuple(map(int, self.size.split("x"))) if self.size else None

        Image(size, self.ratio).show(filename, self.double)  # type: ignore[arg-type]


__all__ = [
    "Image",
    "ShowImageApp",
]


def __dir__() -> list[str]:
    return list(__all__)


if __name__ == "__main__":
    ShowImageApp.run()
