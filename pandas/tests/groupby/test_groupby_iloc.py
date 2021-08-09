""" Test positional grouped indexing with iloc GH#42864"""

import pandas as pd
import pandas._testing as tm
import random


def test_doc_examples():
    """ Test the examples in the documentation"""

    df = pd.DataFrame(
        [['a', 1], ['a', 2], ['a', 3], ['b', 4], ['b', 5]], columns=['A', 'B']
    )

    grouped = df.groupby('A')
    result = grouped.iloc[1:2, :]
    expected = pd.DataFrame([['a', 2], ['b', 5]], columns=['A', 'B'], index=[1, 4])

    tm.assert_frame_equal(result, expected)

    result = grouped.iloc[:-1, -1:]
    expected = pd.DataFrame([1, 2, 4], columns=['B'], index=[0, 1, 3])

    tm.assert_frame_equal(result, expected)


def test_multiindex():
    """ Test the multiindex mentioned as the use-case in the documentation """

    def make_df_from_data(data):
        rows = {}
        for date in dates:
            for level in data[date]:
                rows[(date, level[0])] = {'A': level[1], 'B': level[2]}

        df = pd.DataFrame.from_dict(rows, orient='index')
        df.index.names = ('Date', 'Item')
        return df

    ndates = 1000
    nitems = 40
    dates = pd.date_range("20130101", periods=ndates, freq='D')
    items = [f'item {i}' for i in range(nitems)]

    data = {}
    for date in dates:
        levels = [
            (item, random.randint(0, 10000) / 100, random.randint(0, 10000) / 100) for item in items
        ]
        levels.sort(key=lambda x: x[1])
        data[date] = levels

    df = make_df_from_data(data)
    result = df.groupby('Date').iloc[3:7]

    sliced = {date: data[date][3:7] for date in dates}
    expected = make_df_from_data(sliced)

    tm.assert_frame_equal(result, expected)


def test_against_head_and_tail():
    """ Test gives the same results as grouped head and tail"""

    n_groups = 100
    n_rows_per_group = 30

    data = {
        'group': [f'group {g}' for j in range(n_rows_per_group) for g in range(n_groups)],
        'value': [
            random.randint(0, 10000) / 100 for j in range(n_rows_per_group) for g in range(n_groups)
        ]
    }
    df = pd.DataFrame(data)
    grouped = df.groupby('group')

    for i in [1, 5, 29, 30, 31, 1000]:
        result = grouped.iloc[:i, :]
        expected = grouped.head(i)

        tm.assert_frame_equal(result, expected)

        result = grouped.iloc[-i:, :]
        expected = grouped.tail(i)

        tm.assert_frame_equal(result, expected)


def test_against_df_iloc():
    """ Test that a single group gives the same results as DataFame.iloc"""

    n_rows_per_group = 30

    data = {
        'group': [f'group 0' for j in range(n_rows_per_group)],
        'value': [random.randint(0, 10000) / 100 for j in range(n_rows_per_group)]
    }
    df = pd.DataFrame(data)
    grouped = df.groupby('group')

    for start in [None, 0, 1, 10, 29, 30, 1000, -1, -10, -29, -30, -1000]:
        for stop in [None, 0, 1, 10, 29, 30, 1000, -1, -10, -29, -30, -1000]:
            for step in [None, 1, 2, 3, 10, 29, 30, 100]:
                result = grouped.iloc[start:stop:step, :]
                expected = df.iloc[start:stop:step, :]

                tm.assert_frame_equal(result, expected)


def test_series():
    """ Test grouped Series"""

    ser = pd.Series([1, 2, 3, 4, 5], index=['a', 'a', 'a', 'b', 'b'])
    grouped = ser.groupby(level=0)
    result = grouped.iloc[1:2]
    expected = pd.Series([2, 5], index=['a', 'b'])

    tm.assert_series_equal(result, expected)


def test_step():
    """ Test grouped slice with step"""

    data = [['x', f'x{i}'] for i in range(5)]
    data += [['y', f'y{i}'] for i in range(4)]
    data += [['z', f'z{i}'] for i in range(3)]
    df = pd.DataFrame(data, columns=['A', 'B'])

    grouped = df.groupby('A')

    for step in [1, 2, 3, 4, 5]:
        result = grouped.iloc[::step, :]

        data = [['x', f'x{i}'] for i in range(0, 5, step)]
        data += [['y', f'y{i}'] for i in range(0, 4, step)]
        data += [['z', f'z{i}'] for i in range(0, 3, step)]

        index = [i for i in range(0, 5, step)]
        index += [5 + i for i in range(0, 4, step)]
        index += [9 + i for i in range(0, 3, step)]

        expected = pd.DataFrame(data, columns=['A', 'B'], index=index)

        tm.assert_frame_equal(result, expected)
    