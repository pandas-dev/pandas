import nose

from pandas import DataFrame
import numpy as np

import pandas.util.testing as tm

from pandas.io.json import json_normalize, nested_to_record

def _assert_equal_data(left, right):
    if not left.columns.equals(right.columns):
        left = left.reindex(columns=right.columns)

    tm.assert_frame_equal(left, right)


class TestJSONNormalize(tm.TestCase):

    def setUp(self):
        self.state_data = [
             {'counties': [{'name': 'Dade', 'population': 12345},
                           {'name': 'Broward', 'population': 40000},
                           {'name': 'Palm Beach', 'population': 60000}],
              'info': {'governor': 'Rick Scott'},
              'shortname': 'FL',
              'state': 'Florida'},
             {'counties': [{'name': 'Summit', 'population': 1234},
                           {'name': 'Cuyahoga', 'population': 1337}],
              'info': {'governor': 'John Kasich'},
              'shortname': 'OH',
              'state': 'Ohio'}]

    def test_simple_records(self):
        recs = [{'a': 1, 'b': 2, 'c': 3},
                {'a': 4, 'b': 5, 'c': 6},
                {'a': 7, 'b': 8, 'c': 9},
                {'a': 10, 'b': 11, 'c': 12}]

        result = json_normalize(recs)
        expected = DataFrame(recs)

        tm.assert_frame_equal(result, expected)

    def test_simple_normalize(self):
        result = json_normalize(self.state_data[0], 'counties')
        expected = DataFrame(self.state_data[0]['counties'])
        tm.assert_frame_equal(result, expected)

        result = json_normalize(self.state_data, 'counties')

        expected = []
        for rec in self.state_data:
            expected.extend(rec['counties'])
        expected = DataFrame(expected)

        tm.assert_frame_equal(result, expected)

        result = json_normalize(self.state_data, 'counties', meta='state')
        expected['state'] = np.array(['Florida', 'Ohio']).repeat([3, 2])

        tm.assert_frame_equal(result, expected)

    def test_more_deeply_nested(self):
        data = [{'country': 'USA',
                 'states': [{'name': 'California',
                             'cities': [{'name': 'San Francisco',
                                         'pop': 12345},
                                        {'name': 'Los Angeles',
                                         'pop': 12346}]
                            },
                            {'name': 'Ohio',
                             'cities': [{'name': 'Columbus',
                                         'pop': 1234},
                                        {'name': 'Cleveland',
                                         'pop': 1236}]}
                           ]
                 },
                {'country': 'Germany',
                 'states': [{'name': 'Bayern',
                             'cities': [{'name': 'Munich', 'pop': 12347}]
                            },
                            {'name': 'Nordrhein-Westfalen',
                             'cities': [{'name': 'Duesseldorf', 'pop': 1238},
                                        {'name': 'Koeln', 'pop': 1239}]}
                           ]
                 }
                ]

        result = json_normalize(data, ['states', 'cities'],
                                meta=['country', ['states', 'name']])
                                # meta_prefix={'states': 'state_'})

        ex_data = {'country': ['USA'] * 4 + ['Germany'] * 3,
                   'states.name': ['California', 'California', 'Ohio', 'Ohio',
                                   'Bayern', 'Nordrhein-Westfalen',
                                   'Nordrhein-Westfalen'],
                   'name': ['San Francisco', 'Los Angeles', 'Columbus',
                            'Cleveland', 'Munich', 'Duesseldorf', 'Koeln'],
                   'pop': [12345, 12346, 1234, 1236, 12347, 1238, 1239]}

        expected = DataFrame(ex_data, columns=result.columns)
        tm.assert_frame_equal(result, expected)

    def test_shallow_nested(self):
        data = [{'state': 'Florida',
                 'shortname': 'FL',
                 'info': {
                      'governor': 'Rick Scott'
                 },
                 'counties': [{'name': 'Dade', 'population': 12345},
                             {'name': 'Broward', 'population': 40000},
                             {'name': 'Palm Beach', 'population': 60000}]},
                {'state': 'Ohio',
                 'shortname': 'OH',
                 'info': {
                      'governor': 'John Kasich'
                 },
                 'counties': [{'name': 'Summit', 'population': 1234},
                              {'name': 'Cuyahoga', 'population': 1337}]}]

        result = json_normalize(data, 'counties',
                                ['state', 'shortname',
                                 ['info', 'governor']])
        ex_data = {'name': ['Dade', 'Broward', 'Palm Beach', 'Summit',
                            'Cuyahoga'],
                   'state': ['Florida'] * 3 + ['Ohio'] * 2,
                   'shortname': ['FL', 'FL', 'FL', 'OH', 'OH'],
                   'info.governor': ['Rick Scott'] * 3 + ['John Kasich'] * 2,
                   'population': [12345, 40000, 60000, 1234, 1337]}
        expected = DataFrame(ex_data, columns=result.columns)
        tm.assert_frame_equal(result, expected)

    def test_meta_name_conflict(self):
        data = [{'foo': 'hello',
                 'bar': 'there',
                 'data': [{'foo': 'something', 'bar': 'else'},
                          {'foo': 'something2', 'bar': 'else2'}]}]

        self.assertRaises(ValueError, json_normalize, data,
                          'data', meta=['foo', 'bar'])

        result = json_normalize(data, 'data', meta=['foo', 'bar'],
                                meta_prefix='meta')

        for val in ['metafoo', 'metabar', 'foo', 'bar']:
            self.assertTrue(val in result)

    def test_record_prefix(self):
        result = json_normalize(self.state_data[0], 'counties')
        expected = DataFrame(self.state_data[0]['counties'])
        tm.assert_frame_equal(result, expected)

        result = json_normalize(self.state_data, 'counties',
                                meta='state',
                                record_prefix='county_')

        expected = []
        for rec in self.state_data:
            expected.extend(rec['counties'])
        expected = DataFrame(expected)
        expected = expected.rename(columns=lambda x: 'county_' + x)
        expected['state'] = np.array(['Florida', 'Ohio']).repeat([3, 2])

        tm.assert_frame_equal(result, expected)


class TestNestedToRecord(tm.TestCase):

    def test_flat_stays_flat(self):
        recs = [dict(flat1=1,flat2=2),
                dict(flat1=3,flat2=4),
                ]

        result = nested_to_record(recs)
        expected = recs
        self.assertEqual(result, expected)

    def test_one_level_deep_flattens(self):
        data = dict(flat1=1,
                    dict1=dict(c=1,d=2))

        result = nested_to_record(data)
        expected =     {'dict1.c': 1,
             'dict1.d': 2,
             'flat1': 1}

        self.assertEqual(result,expected)

    def test_nested_flattens(self):
        data = dict(flat1=1,
                    dict1=dict(c=1,d=2),
                    nested=dict(e=dict(c=1,d=2),
                                d=2))

        result = nested_to_record(data)
        expected =     {'dict1.c': 1,
             'dict1.d': 2,
             'flat1': 1,
             'nested.d': 2,
             'nested.e.c': 1,
             'nested.e.d': 2}

        self.assertEqual(result,expected)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb',
                         '--pdb-failure', '-s'], exit=False)
