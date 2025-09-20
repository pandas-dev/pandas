from __future__ import annotations

import sys

from plumbum import colors

from .. import cli
from .termsize import get_terminal_size


class Image:
    __slots__ = "size char_ratio".split()

    def __init__(self, size=None, char_ratio=2.45):
        self.size = size
        self.char_ratio = char_ratio

    def best_aspect(self, orig, term):
        """Select a best possible size matching the original aspect ratio.
        Size is width, height.
        The char_ratio option gives the height of each char with respect
        to its width, zero for no effect."""

        if not self.char_ratio:  # Don't use if char ratio is 0
            return term

        orig_ratio = orig[0] / orig[1] / self.char_ratio

        if int(term[1] / orig_ratio) <= term[0]:
            return int(term[1] / orig_ratio), term[1]

        return term[0], int(term[0] * orig_ratio)

    def show(self, filename, double=False):
        """Display an image on the command line. Can select a size or show in double resolution."""

        import PIL.Image

        return (
            self.show_pil_double(PIL.Image.open(filename))
            if double
            else self.show_pil(PIL.Image.open(filename))
        )

    def _init_size(self, im):
        """Return the expected image size"""
        if self.size is None:
            term_size = get_terminal_size()
            return self.best_aspect(im.size, term_size)

        return self.size

    def show_pil(self, im):
        "Standard show routine"
        size = self._init_size(im)
        new_im = im.resize(size).convert("RGB")

        for y in range(size[1]):
            for x in range(size[0] - 1):
                pix = new_im.getpixel((x, y))
                sys.stdout.write(colors.bg.rgb(*pix) + " ")  # '\u2588'
            sys.stdout.write(colors.reset + " \n")
        sys.stdout.write(colors.reset + "\n")
        sys.stdout.flush()

    def show_pil_double(self, im):
        "Show double resolution on some fonts"

        size = self._init_size(im)
        size = (size[0], size[1] * 2)
        new_im = im.resize(size).convert("RGB")

        for y in range(size[1] // 2):
            for x in range(size[0] - 1):
                pix = new_im.getpixel((x, y * 2))
                pixl = new_im.getpixel((x, y * 2 + 1))
                sys.stdout.write(
                    (colors.bg.rgb(*pixl) & colors.fg.rgb(*pix)) + "\u2580"
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
    def colors_set(self, n):  # pylint: disable=no-self-use
        colors.use_color = n

    size = cli.SwitchAttr(["-s", "--size"], help="Size, should be in the form 100x150")

    ratio = cli.SwitchAttr(
        ["--ratio"], float, default=2.45, help="Aspect ratio of the font"
    )

    @cli.positional(cli.ExistingFile)
    def main(self, filename):
        size = None
        if self.size:
            size = map(int, self.size.split("x"))

        Image(size, self.ratio).show(filename, self.double)


if __name__ == "__main__":
    ShowImageApp.run()
