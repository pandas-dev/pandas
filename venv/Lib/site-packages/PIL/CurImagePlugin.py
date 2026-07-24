#
# The Python Imaging Library.
# $Id$
#
# Windows Cursor support for PIL
#
# notes:
#       uses BmpImagePlugin.py to read the bitmap data.
#
# history:
#       96-05-27 fl     Created
#
# Copyright (c) Secret Labs AB 1997.
# Copyright (c) Fredrik Lundh 1996.
#
# See the README file for information on usage and redistribution.
#
from __future__ import annotations

from . import BmpImagePlugin, Image, ImageOps
from ._binary import i16le as i16
from ._binary import i32le as i32

#
# --------------------------------------------------------------------


def _accept(prefix: bytes) -> bool:
    return prefix.startswith(b"\0\0\2\0")


##
# Image plugin for Windows Cursor files.


class CurImageFile(BmpImagePlugin.BmpImageFile):
    format = "CUR"
    format_description = "Windows Cursor"

    def _open(self) -> None:
        assert self.fp is not None
        offset = self.fp.tell()

        # check magic
        s = self.fp.read(6)
        if not _accept(s):
            msg = "not a CUR file"
            raise SyntaxError(msg)

        # pick the largest cursor in the file
        m = b""
        for i in range(i16(s, 4)):
            s = self.fp.read(16)
            if not m:
                m = s
            elif s[0] > m[0] and s[1] > m[1]:
                m = s
        if not m:
            msg = "No cursors were found"
            raise TypeError(msg)

        # load as bitmap
        self._bitmap(i32(m, 12) + offset)
        self._masked = self.mode in ("1", "L")
        if self._masked:
            self._rawmode = self.mode
            self._mode = "LA"

        # patch up the bitmap height
        self._size = self.width, self.height // 2
        if not self._masked:
            self.tile = [self.tile[0]._replace(extents=(0, 0) + self.size)]

    def load_prepare(self) -> None:
        if self._masked:
            self._mode = self._rawmode
            self._size = self.width, self.height * 2
        super().load_prepare()

    def load_end(self) -> None:
        if not self._masked:
            return
        self._mode = "LA"
        new_height = self.height // 2

        and_mask = self.im.crop((0, 0, self.width, new_height))
        xor_mask = self.im.crop((0, new_height, self.width, self.height))

        self._size = self.width, new_height
        self._im = Image.core.fill(self.mode, self.size)
        self._im.paste(
            xor_mask.convert(self.mode),
            (0, 0) + self.size,
            ImageOps.invert(Image.Image()._new(and_mask)).im,
        )


#
# --------------------------------------------------------------------

Image.register_open(CurImageFile.format, CurImageFile, _accept)

Image.register_extension(CurImageFile.format, ".cur")
