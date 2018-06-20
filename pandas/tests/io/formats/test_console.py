import pytest
from functools import partial

from pandas.io.formats.console import detect_console_encoding, locale


def mock_raises_exception(error=Exception):
    raise error


def test_detect_console_encoding_stdout(monkeypatch):
    monkeypatch.setattr('sys.stdin.encoding', '')
    monkeypatch.setattr('sys.stdout.encoding', 'foo')
    assert detect_console_encoding() == 'foo'


def test_detect_console_encoding_stdin(monkeypatch):
    monkeypatch.setattr('sys.stdout.encoding', '')
    monkeypatch.setattr('sys.stdin.encoding', 'foo')
    assert detect_console_encoding() == 'foo'


@pytest.mark.parametrize('error', [AttributeError, IOError])
def test_detect_console_encoding_stdout_error_uses_locale(monkeypatch, error):
    monkeypatch.setattr('locale.getpreferredencoding', lambda: 'foo')
    monkeypatch.setattr('sys.stdout', partial(mock_raises_exception, error))
    assert detect_console_encoding() == 'foo'


def test_detect_console_encoding_sys_default_encoding(monkeypatch):
    monkeypatch.setattr('locale.getpreferredencoding', mock_raises_exception)
    monkeypatch.setattr('sys.stdout', mock_raises_exception)
    monkeypatch.setattr('sys.getdefaultencoding', lambda: 'sysDefaultEncoding')
    assert detect_console_encoding() == 'sysDefaultEncoding'
