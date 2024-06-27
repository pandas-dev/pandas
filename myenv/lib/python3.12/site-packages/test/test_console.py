# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import sys
import locale
import itertools

import pytest
from asv_runner.console import _write_with_fallback, color_print

from asv.console import log


def test_write_with_fallback(tmpdir, capfd):
    tmpdir = str(tmpdir)

    def check_write(value, expected, stream_encoding, preferred_encoding):
        old_getpreferredencoding = locale.getpreferredencoding
        try:
            locale.getpreferredencoding = lambda: preferred_encoding

            # Check writing to io.StringIO
            stream = io.StringIO()
            _write_with_fallback(value, stream)
            assert stream.getvalue() == value

            # Check writing to a text stream
            buf = io.BytesIO()
            stream = io.TextIOWrapper(buf, encoding=stream_encoding)
            _write_with_fallback(value, stream)
            stream.flush()
            got = buf.getvalue()
            assert got == expected

            # Writing to io.BytesIO
            stream = io.BytesIO()
            with pytest.raises(TypeError):
                _write_with_fallback(value, stream)

            # Check writing to a file
            fn = os.path.join(tmpdir, 'tmp.txt')
            with io.open(fn, 'w', encoding=stream_encoding) as stream:
                _write_with_fallback(value, stream)
            with open(fn, 'rb') as stream:
                got = stream.read()
                assert got == expected
        finally:
            locale.getpreferredencoding = old_getpreferredencoding

    # What is printed should follow the following rules:
    #
    # - Try printing in stream encoding.
    # - Try printing in locale preferred encoding.
    # - Otherwise, map characters produced by asv to ascii equivalents, and
    #   - Try to print in latin1
    #   - Try to print in ascii, replacing all non-ascii characters
    encodings = ['utf-8', 'latin1', 'ascii', 'euc-jp']
    strings = ["helloμ", "hello·", "hello難", "helloä", "hello±"]
    repmap = {"helloμ": "hellou", "hello·": "hello-", "hello±": "hello~"}

    for pref_enc, stream_enc, s in itertools.product(encodings, encodings, strings):
        expected = None
        encodings = [stream_enc, pref_enc]
        for enc in encodings:
            try:
                expected = s.encode(enc)
                break
            except UnicodeError:
                pass
        else:
            s2 = repmap.get(s, s)
            expected = s2.encode(pref_enc, errors='replace')

        check_write(s, expected, stream_enc, pref_enc)

    # Should bail out on bytes input
    with pytest.raises(ValueError):
        _write_with_fallback(b"a", sys.stdout)


def test_color_print_nofail(capfd):
    # Try out color print

    color_print("hello", "red")
    color_print("indeed難", "blue")
    with pytest.raises(ValueError):
        color_print(b"really\xfe", "green", "not really")

    out, err = capfd.readouterr()
    assert 'hello' in out
    assert 'indeed' in out


def test_log_indent(capsys):
    log.set_nitems(0)
    log.info("First\nSecond")

    out, err = capsys.readouterr()
    lines = out.lstrip().splitlines()
    assert lines[0].index('First') == lines[1].index('Second')

    log.set_nitems(1)
    log.info("First\nSecond")

    out, err = capsys.readouterr()
    lines = out.lstrip().splitlines()
    assert lines[0].index('First') == lines[1].index('Second')
