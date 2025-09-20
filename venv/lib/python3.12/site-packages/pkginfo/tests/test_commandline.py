import collections
import io
import json as json_module
import sys
from unittest import mock

import pytest

def test__parse_options_empty():
    from pkginfo.commandline import __doc__ as usage
    from pkginfo.commandline import _parse_options

    firstline = usage.splitlines()[0]
    buf = io.StringIO()

    with mock.patch.object(sys, "stderr", buf):
        with pytest.raises(SystemExit):
            _parse_options([])

    assert(firstline in buf.getvalue())

def test__parse_options_nonempty():
    from pkginfo.commandline import _parse_options

    options, args = _parse_options(['foo'])
    assert(args == ['foo'])

def _make_base(options):
    from pkginfo.commandline import Base

    return Base(options)

def test_base_ctor_defaults():
    base = _make_base(_Options(fields=()))
    assert(base._fields is None)

def test_base_ctor_w_fields():
    fields = object()
    base = _make_base(_Options(fields=fields))
    assert(base._fields is fields)

def _capture_output(func, *args, **kw):
    buf = io.StringIO()

    with mock.patch.object(sys, "stdout", buf):
        func(*args, **kw)

    return buf.getvalue()

def _no_output(simple, meta):
    with mock.patch.object(sys, "stdout", object()):  # raise if write
        simple(meta)

def _make_simple(options):
    from pkginfo.commandline import Simple
    return Simple(options)

def test_simple___init__():
    simple = _make_simple(_Options(fields=None, skip=True))
    assert(simple._skip)

def test_simple__call___w_empty_fields():
    simple = _make_simple(_Options(fields=(), skip=False))
    meta = _Meta()
    _no_output(simple, meta)

def test_simple__call___w_skip_and_value_None_no_fields():
    simple = _make_simple(_Options(fields=(), skip=True))
    meta = _Meta(foo=None)
    _no_output(simple, meta)

def test_simple__call___w_skip_and_value_empty_tuple_explicit_fields():
    simple = _make_simple(_Options(fields=('foo',), skip=True))
    meta = _Meta(foo=(), bar='Bar')
    _no_output(simple, meta)

def test_simple__call___w_skip_but_values_explicit_fields():
    simple = _make_simple(_Options(fields=('foo',), skip=True))
    meta = _Meta(foo='Foo')
    output = _capture_output(simple, meta)
    assert(output == 'foo: Foo\n')

def _make_single_line(options):
    from pkginfo.commandline import SingleLine

    return SingleLine(options)

def test_singleline__init___():
    single = _make_single_line(
        _Options(fields=None, item_delim='I', sequence_delim='S'))
    assert(single._item_delim == 'I')
    assert(single._sequence_delim == 'S')

def test_singleline__call__wo_fields_wo_list():
    single = _make_single_line(
        _Options(fields=(), item_delim='|',
                    sequence_delim=object()))  # raise if used
    meta = _Meta(foo='Foo', bar='Bar')
    output = _capture_output(single, meta)
    assert(output == 'Bar|Foo\n')

def test_singleline__call__w_fields_w_list():
    single = _make_single_line(
        _Options(fields=('foo', 'bar'), item_delim='|',
                    sequence_delim='*'))
    meta = _Meta(foo='Foo', bar=['Bar1', 'Bar2'], baz='Baz')
    output = _capture_output(single, meta)
    assert(output == 'Foo|Bar1*Bar2\n')

def _make_csv(options):
    from pkginfo.commandline import CSV

    return CSV(options)

def test_csv__init___():
    csv = _make_csv(
        _Options(fields=None, sequence_delim='S'))
    assert(csv._sequence_delim == 'S')

def test_csv__call__wo_fields_wo_list():
    meta = _Meta(foo='Foo', bar='Bar')
    csv = _make_csv(
        _Options(fields=None,
                    sequence_delim=object()))  # raise if used
    output = _capture_output(csv, meta)
    assert(output == 'bar,foo\r\nBar,Foo\r\n')

def test_csv__call__w_fields_w_list():
    meta = _Meta(foo='Foo', bar=['Bar1', 'Bar2'], baz='Baz')
    csv = _make_csv(
        _Options(fields=('foo', 'bar'), item_delim='|',
                    sequence_delim='*'))
    output = _capture_output(csv, meta)
    assert(output == 'foo,bar\r\nFoo,Bar1*Bar2\r\n')

def _make_ini(options):
    from pkginfo.commandline import INI

    return INI(options)

def test_ini__call___duplicate():
    ini = _make_ini(_Options(fields=('foo',)))
    meta = _Meta(name='foo', version='0.1', foo='Foo')
    ini._parser.add_section('foo-0.1')
    with pytest.raises(ValueError):
        ini(meta)

def test_ini__call___wo_fields_wo_list():
    ini = _make_ini(_Options(fields=None))
    meta = _Meta(name='foo', version='0.1', foo='Foo')
    ini(meta)
    cp = ini._parser
    assert(cp.sections() == ['foo-0.1'])
    assert(sorted(cp.options('foo-0.1')) == ['foo', 'name', 'version'])
    assert(cp.get('foo-0.1', 'name') == 'foo')
    assert(cp.get('foo-0.1', 'version') == '0.1')
    assert(cp.get('foo-0.1', 'foo') == 'Foo')

def test_ini__call___w_fields_w_list():
    ini = _make_ini(_Options(fields=('foo', 'bar')))
    meta = _Meta(
        name='foo',
        version='0.1',
        foo='Foo',
        bar=['Bar1', 'Bar2'],
        baz='Baz',
    )
    ini(meta)
    cp = ini._parser
    assert(cp.sections() == ['foo-0.1'])
    assert(sorted(cp.options('foo-0.1')) == ['bar', 'foo'])
    assert(cp.get('foo-0.1', 'foo') == 'Foo')
    assert(cp.get('foo-0.1', 'bar') == 'Bar1\n\tBar2')

def _make_json(options):
    from pkginfo.commandline import JSON

    return JSON(options)

def test_json__call___duplicate_with_meta_and_fields():
    json = _make_json(_Options(fields=('name',)))
    meta = _Meta(name='foo', version='0.1', foo='Foo')
    json._mapping['name'] = 'foo'
    with pytest.raises(ValueError):
        json(meta)

def test_json__call___duplicate_with_meta_wo_fields():
    json = _make_json(_Options(fields=None))
    meta = _Meta(name='foo', version='0.1', foo='Foo')
    json._mapping['name'] = 'foo'
    with pytest.raises(ValueError):
        json(meta)

def test_json__call___wo_fields_wo_list():

    json = _make_json(_Options(fields=None))
    meta = _Meta(name='foo', version='0.1', foo='Foo')
    json(meta)
    expected = collections.OrderedDict([
        ('foo', 'Foo'), ('name', 'foo'), ('version', '0.1')])
    assert(expected == json._mapping)

def test_json__call___w_fields_w_list():
    json = _make_json(_Options(fields=('foo', 'bar')))
    meta = _Meta(name='foo', version='0.1',
                    foo='Foo', bar=['Bar1', 'Bar2'], baz='Baz')
    json(meta)
    expected = collections.OrderedDict([
        ('foo', 'Foo'), ('bar', ['Bar1', 'Bar2'])])
    assert(expected == json._mapping)

def test_json__call___output():
    json = _make_json(_Options(fields=None))
    meta = _Meta(name='foo', version='0.1', foo='Foo')
    json(meta)
    output = _capture_output(json.finish)
    output = json_module.loads(
        output, object_pairs_hook=collections.OrderedDict)
    expected = collections.OrderedDict([
        ('foo', 'Foo'), ('name', 'foo'), ('version', '0.1')])
    assert(expected == output)

def _call_main(args, monkey='simple'):
    from pkginfo.commandline import main

    formatter = mock.Mock(spec=["finish"])

    with mock.patch.dict(
        "pkginfo.commandline._FORMATTERS",
        simple=lambda *options: formatter,
    ):
        main(args)

    return formatter

def test_main_w_missing_dist():
    from pkginfo import commandline as MUT

    with mock.patch("pkginfo.commandline.get_metadata") as _get_metadata:
        _get_metadata.return_value = None
        formatter = _call_main(['foo'])

    formatter.assert_not_called()
    formatter.finish.assert_called_once_with()
    _get_metadata.assert_called_once_with("foo", None)

def test_main_w_dist_wo_download_url():
    from pkginfo import commandline as MUT

    meta = _Meta(download_url=None)

    with mock.patch("pkginfo.commandline.get_metadata") as _get_metadata:
        _get_metadata.return_value = meta
        formatter = _call_main(
            ['-d', 'http://example.com', '/path/to/foo'])

    formatter.assert_called_once_with(meta)
    formatter.finish.assert_called_once_with()
    _get_metadata.assert_called_once_with("/path/to/foo", None)

    assert(meta.download_url == 'http://example.com/foo')

def test_main_w_dist_w_download_url():
    from pkginfo import commandline as MUT

    meta = _Meta(download_url='http://example.com/dist/foo')

    with mock.patch("pkginfo.commandline.get_metadata") as _get_metadata:
        _get_metadata.return_value = meta
        formatter = _call_main(
            ['-d', 'http://example.com', '/path/to/foo'])

    formatter.assert_called_once_with(meta)
    formatter.finish.assert_called_once_with()
    _get_metadata.assert_called_once_with("/path/to/foo", None)

    assert(meta.download_url == 'http://example.com/dist/foo')

class _Options(object):

    def __init__(self, **kw):
        for k in kw:
            self.__dict__[k] = kw[k]

class _Meta(object):

    def __init__(self, **kw):
        for k in kw:
            self.__dict__[k] = kw[k]

    def __iter__(self):
        return iter(sorted(self.__dict__))
