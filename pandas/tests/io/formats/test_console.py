import pytest

from pandas.io.formats.console import detect_console_encoding


@pytest.mark.parametrize('empty,filled', [
    ['stdin', 'stdout'],
    ['stdout', 'stdin']
])
def test_detect_console_encoding_from_stdout_stdin(monkeypatch, empty, filled):
    # Ensures that when sys.stdout.encoding or sys.stdin.encoding is used when
    # they have values filled.
    # GH 21552
    class MockEncoding(object):
        def __init__(self, encoding):
            super(MockEncoding, self).__init__()
            self.encoding = encoding

    with monkeypatch.context() as context:
        context.setattr('sys.{}'.format(empty), MockEncoding(''))
        context.setattr('sys.{}'.format(filled), MockEncoding(filled))
        assert detect_console_encoding() == filled


@pytest.mark.parametrize('stdEncoding', [
    pytest.raises(AttributeError),
    pytest.raises(IOError),
    lambda: 'ascii'
])
def test_detect_console_encoding_fallback_to_locale(monkeypatch, stdEncoding):
    # GH 21552
    with monkeypatch.context() as context:
        context.setattr('locale.getpreferredencoding', lambda: 'foo')
        context.setattr('sys.stdout', stdEncoding)
        assert detect_console_encoding() == 'foo'


@pytest.mark.parametrize('std,locale', [
    ['ascii', 'ascii'],
    ['ascii', pytest.raises(Exception)],
    [pytest.raises(Exception), 'ascii'],
    [pytest.raises(Exception), pytest.raises(Exception)]
])
def test_detect_console_encoding_fallback_to_default(monkeypatch, std, locale):
    # When both the stdout/stdin encoding and locale preferred encoding checks
    # fail (or return 'ascii', we should default to the sys default encoding.
    # GH 21552
    with monkeypatch.context() as context:
        context.setattr('locale.getpreferredencoding', locale)
        context.setattr('sys.stdout', std)
        context.setattr('sys.getdefaultencoding', lambda: 'sysDefaultEncoding')
        assert detect_console_encoding() == 'sysDefaultEncoding'
