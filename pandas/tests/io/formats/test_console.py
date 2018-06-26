import pytest
from functools import partial

from pandas.io.formats.console import detect_console_encoding, locale


def mock_raises_exception(error=Exception):
    raise error


@pytest.mark.parametrize('empty,filled', [
    ['stdin', 'stdout'],
    ['stdout', 'stdin']
])
def test_detect_console_encoding_from_stdout_stdin(monkeypatch, empty, filled):
    # Ensures that when sys.stdout.encoding or sys.stdin.encoding is used when
    # they have values filled.
    class MockEncoding(object):
        def __init__(self, encoding):
            super(MockEncoding, self).__init__()
            self.encoding = encoding

    monkeypatch.setattr('sys.{}'.format(empty), MockEncoding(''))
    monkeypatch.setattr('sys.{}'.format(filled), MockEncoding(filled))
    assert detect_console_encoding() == filled


@pytest.mark.parametrize('stdEncoding', [
    partial(mock_raises_exception, AttributeError),
    partial(mock_raises_exception, IOError),
    lambda: 'ascii'
])
def test_detect_console_encoding_fallback_to_locale(monkeypatch, stdEncoding):
    monkeypatch.setattr('locale.getpreferredencoding', lambda: 'foo')
    monkeypatch.setattr('sys.stdout', stdEncoding)
    assert detect_console_encoding() == 'foo'


@pytest.mark.parametrize('std,locale', [
    ['ascii', 'ascii'],
    ['ascii', mock_raises_exception],
    [mock_raises_exception, 'ascii'],
    [mock_raises_exception, mock_raises_exception]
])
def test_detect_console_encoding_fallback_to_default(monkeypatch, std, locale):
    # When both the stdout/stdin encoding and locale preferred encoding checks
    # fail (or return 'ascii', we should default to the sys default encoding.
    monkeypatch.setattr('sys.stdout', std)
    monkeypatch.setattr('locale.getpreferredencoding', locale)
    monkeypatch.setattr('sys.getdefaultencoding', lambda: 'sysDefaultEncoding')
    assert detect_console_encoding() == 'sysDefaultEncoding'
