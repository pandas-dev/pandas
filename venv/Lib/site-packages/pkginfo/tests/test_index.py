import pytest

def _make_index():
    from pkginfo.index import Index

    return Index()

def _make_dummy():
    from pkginfo.distribution import Distribution

    class DummyDistribution(Distribution):
        name = 'dummy'
        version = '1.0'

    return DummyDistribution()

class NotADistribution:
    name = 'not_a_distro'
    version = '1.0'

def test_index_empty():
    index = _make_index()
    assert(len(index) == 0)
    assert(len(index.keys()) == 0)
    assert(len(index.values()) == 0)
    assert(len(index.items()) == 0)

def test_index___getitem___miss():
    index = _make_index()

    with pytest.raises(KeyError):
        index['nonesuch']

def test_index___setitem___value_not_dist():

    not_a_dist = NotADistribution()
    index = _make_index()

    with pytest.raises(ValueError):
        index['dummy-1.0'] = not_a_dist

def test_index___setitem___bad_key():
    index = _make_index()
    dummy = _make_dummy()

    with pytest.raises(ValueError):
        index['nonesuch'] = dummy

def test_index___setitem___valid_key():
    index = _make_index()
    dummy = _make_dummy()

    index['dummy-1.0'] = dummy

    assert(index['dummy-1.0'] is dummy)
    assert(len(index) == 1)
    assert(len(index.keys()) == 1)
    assert(list(index.keys())[0] == 'dummy-1.0')
    assert(len(index.values()) == 1)
    assert(list(index.values())[0] == dummy)
    assert(len(index.items()) == 1)
    assert(list(index.items())[0] == ('dummy-1.0', dummy))

def test_index_add_not_dist():
    index = _make_index()
    dummy = NotADistribution()

    with pytest.raises(ValueError):
        index.add(dummy)

def test_index_add_valid_dist():
    index = _make_index()
    dummy = _make_dummy()

    index.add(dummy)

    assert(index['dummy-1.0'] is dummy)
    assert(len(index) == 1)
    assert(len(index.keys()) == 1)
    assert(list(index.keys())[0] == 'dummy-1.0')
    assert(len(index.values()) == 1)
    assert(list(index.values())[0] == dummy)
    assert(len(index.items()) == 1)
    assert(list(index.items())[0] == ('dummy-1.0', dummy))
