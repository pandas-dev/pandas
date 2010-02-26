import pandas.core.pytools as pytools

def test_rands():
    r = pytools.rands(10)
    assert(len(r) == 10)

def test_adjoin():
    data = [['a', 'b', 'c'],
            ['dd', 'ee', 'ff'],
            ['ggg', 'hhh', 'iii']]
    expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

    adjoined = pytools.adjoin(2, *data)

    assert(adjoined == expected)

def test_iterpairs():
    data = [1, 2, 3, 4]
    expected = [(1, 2),
                (2, 3),
                (3, 4)]

    result = list(pytools.iterpairs(data))

    assert(result == expected)

def test_indent():
    s = 'a b c\nd e f'
    result = pytools.indent(s, spaces=6)

    assert(result == '      a b c\n      d e f')

def test_banner():
    ban = pytools.banner('hi')
    assert(ban == ('%s\nhi\n%s' % ('=' * 80, '=' * 80)))

def test_map_indices_py():
    data = [4, 3, 2, 1]
    expected = {4 : 0, 3 : 1, 2 : 2, 1 : 3}

    result = pytools.map_indices_py(data)

    assert(result == expected)

def test_union():
    pass

def test_difference():
    pass

def test_intersection():
    pass
