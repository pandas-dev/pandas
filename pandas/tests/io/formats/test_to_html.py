# -*- coding: utf-8 -*-

import re
from datetime import datetime
from io import open

import pytest
import numpy as np
import pandas as pd
from pandas import compat, DataFrame, MultiIndex, option_context, Index
from pandas.compat import u, lrange, StringIO
from pandas.util import testing as tm
import pandas.io.formats.format as fmt


def expected_html(datapath, name):
    """
    Read HTML file from formats data directory.

    Parameters
    ----------
    datapath : pytest fixture
        The datapath fixture injected into a test by pytest.
    name : str
        The name of the HTML file without the suffix.

    Returns
    -------
    str : contents of HTML file.
    """
    filename = '.'.join([name, 'html'])
    filepath = datapath('io', 'formats', 'data', filename)
    with open(filepath, encoding='utf-8') as f:
        html = f.read()
    return html.rstrip()


class TestToHTML(object):

    def test_to_html_with_col_space(self):
        def check_with_width(df, col_space):
            # check that col_space affects HTML generation
            # and be very brittle about it.
            html = df.to_html(col_space=col_space)
            hdrs = [x for x in html.split(r"\n") if re.search(r"<th[>\s]", x)]
            assert len(hdrs) > 0
            for h in hdrs:
                assert "min-width" in h
                assert str(col_space) in h

        df = DataFrame(np.random.random(size=(1, 3)))

        check_with_width(df, 30)
        check_with_width(df, 50)

    def test_to_html_with_empty_string_label(self):
        # GH 3547, to_html regards empty string labels as repeated labels
        data = {'c1': ['a', 'b'], 'c2': ['a', ''], 'data': [1, 2]}
        df = DataFrame(data).set_index(['c1', 'c2'])
        result = df.to_html()
        assert "rowspan" not in result

    def test_to_html_unicode(self, datapath):
        df = DataFrame({u('\u03c3'): np.arange(10.)})
        expected = expected_html(datapath, 'unicode_1')
        assert df.to_html() == expected
        df = DataFrame({'A': [u('\u03c3')]})
        expected = expected_html(datapath, 'unicode_2')
        assert df.to_html() == expected

    def test_to_html_decimal(self, datapath):
        # GH 12031
        df = DataFrame({'A': [6.0, 3.1, 2.2]})
        result = df.to_html(decimal=',')
        expected = expected_html(datapath, 'gh12031_expected_output')
        assert result == expected

    def test_to_html_escaped(self, datapath):
        a = 'str<ing1 &amp;'
        b = 'stri>ng2 &amp;'

        test_dict = {'co<l1': {a: "<type 'str'>",
                               b: "<type 'str'>"},
                     'co>l2': {a: "<type 'str'>",
                               b: "<type 'str'>"}}
        result = DataFrame(test_dict).to_html()
        expected = expected_html(datapath, 'escaped')
        assert result == expected

    def test_to_html_escape_disabled(self, datapath):
        a = 'str<ing1 &amp;'
        b = 'stri>ng2 &amp;'

        test_dict = {'co<l1': {a: "<b>bold</b>",
                               b: "<b>bold</b>"},
                     'co>l2': {a: "<b>bold</b>",
                               b: "<b>bold</b>"}}
        result = DataFrame(test_dict).to_html(escape=False)
        expected = expected_html(datapath, 'escape_disabled')
        assert result == expected

    def test_to_html_multiindex_index_false(self, datapath):
        # GH 8452
        df = DataFrame({
            'a': range(2),
            'b': range(3, 5),
            'c': range(5, 7),
            'd': range(3, 5)
        })
        df.columns = MultiIndex.from_product([['a', 'b'], ['c', 'd']])
        result = df.to_html(index=False)
        expected = expected_html(datapath, 'gh8452_expected_output')
        assert result == expected

        df.index = Index(df.index.values, name='idx')
        result = df.to_html(index=False)
        assert result == expected

    def test_to_html_multiindex_sparsify_false_multi_sparse(self, datapath):
        with option_context('display.multi_sparse', False):
            index = MultiIndex.from_arrays([[0, 0, 1, 1], [0, 1, 0, 1]],
                                           names=['foo', None])

            df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], index=index)
            result = df.to_html()
            expected = expected_html(
                datapath, 'multiindex_sparsify_false_multi_sparse_1')
            assert result == expected

            df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]],
                           columns=index[::2], index=index)
            result = df.to_html()
            expected = expected_html(
                datapath, 'multiindex_sparsify_false_multi_sparse_2')
            assert result == expected

    def test_to_html_multiindex_sparsify(self, datapath):
        index = MultiIndex.from_arrays([[0, 0, 1, 1], [0, 1, 0, 1]],
                                       names=['foo', None])

        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], index=index)
        result = df.to_html()
        expected = expected_html(datapath, 'multiindex_sparsify_1')
        assert result == expected

        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], columns=index[::2],
                       index=index)
        result = df.to_html()
        expected = expected_html(datapath, 'multiindex_sparsify_2')
        assert result == expected

    def test_to_html_multiindex_odd_even_truncate(self, datapath):
        # GH 14882 - Issue on truncation with odd length DataFrame
        mi = MultiIndex.from_product([[100, 200, 300],
                                      [10, 20, 30],
                                      [1, 2, 3, 4, 5, 6, 7]],
                                     names=['a', 'b', 'c'])
        df = DataFrame({'n': range(len(mi))}, index=mi)
        result = df.to_html(max_rows=60)
        expected = expected_html(datapath, 'gh14882_expected_output_1')
        assert result == expected

        # Test that ... appears in a middle level
        result = df.to_html(max_rows=56)
        expected = expected_html(datapath, 'gh14882_expected_output_2')
        assert result == expected

    def test_to_html_index_formatter(self, datapath):
        df = DataFrame([[0, 1], [2, 3], [4, 5], [6, 7]], columns=['foo', None],
                       index=lrange(4))

        f = lambda x: 'abcd' [x]
        result = df.to_html(formatters={'__index__': f})
        expected = expected_html(datapath, 'index_formatter')
        assert result == expected

    def test_to_html_datetime64_monthformatter(self, datapath):
        months = [datetime(2016, 1, 1), datetime(2016, 2, 2)]
        x = DataFrame({'months': months})

        def format_func(x):
            return x.strftime('%Y-%m')
        result = x.to_html(formatters={'months': format_func})
        expected = expected_html(datapath, 'datetime64_monthformatter')
        assert result == expected

    def test_to_html_datetime64_hourformatter(self, datapath):

        x = DataFrame({'hod': pd.to_datetime(['10:10:10.100', '12:12:12.120'],
                                             format='%H:%M:%S.%f')})

        def format_func(x):
            return x.strftime('%H:%M')
        result = x.to_html(formatters={'hod': format_func})
        expected = expected_html(datapath, 'datetime64_hourformatter')
        assert result == expected

    def test_to_html_regression_GH6098(self):
        df = DataFrame({
            u('clé1'): [u('a'), u('a'), u('b'), u('b'), u('a')],
            u('clé2'): [u('1er'), u('2ème'), u('1er'), u('2ème'), u('1er')],
            'données1': np.random.randn(5),
            'données2': np.random.randn(5)})

        # it works
        df.pivot_table(index=[u('clé1')], columns=[u('clé2')])._repr_html_()

    def test_to_html_truncate(self, datapath):
        index = pd.date_range(start='20010101', freq='D', periods=20)
        df = DataFrame(index=index, columns=range(20))
        result = df.to_html(max_rows=8, max_cols=4)
        expected = expected_html(datapath, 'truncate')
        assert result == expected

    def test_to_html_truncate_multi_index(self, datapath):
        arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        df = DataFrame(index=arrays, columns=arrays)
        result = df.to_html(max_rows=7, max_cols=7)
        expected = expected_html(datapath, 'truncate_multi_index')
        assert result == expected

    @pytest.mark.xfail(reason='GH22887 TypeError')
    def test_to_html_truncate_multi_index_sparse_off(self, datapath):
        arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        df = DataFrame(index=arrays, columns=arrays)
        result = df.to_html(max_rows=7, max_cols=7, sparsify=False)
        expected = expected_html(datapath, 'truncate_multi_index_sparse_off')
        assert result == expected

    def test_to_html_border(self):
        df = DataFrame({'A': [1, 2]})
        result = df.to_html()
        assert 'border="1"' in result

    def test_to_html_border_option(self):
        df = DataFrame({'A': [1, 2]})
        with option_context('display.html.border', 0):
            result = df.to_html()
            assert 'border="0"' in result
            assert 'border="0"' in df._repr_html_()

    def test_to_html_border_zero(self):
        df = DataFrame({'A': [1, 2]})
        result = df.to_html(border=0)
        assert 'border="0"' in result

    def test_display_option_warning(self):
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            pd.options.html.border

    def test_to_html(self):
        # big mixed
        biggie = DataFrame({'A': np.random.randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))

        biggie.loc[:20, 'A'] = np.nan
        biggie.loc[:20, 'B'] = np.nan
        s = biggie.to_html()

        buf = StringIO()
        retval = biggie.to_html(buf=buf)
        assert retval is None
        assert buf.getvalue() == s

        assert isinstance(s, compat.string_types)

        biggie.to_html(columns=['B', 'A'], col_space=17)
        biggie.to_html(columns=['B', 'A'],
                       formatters={'A': lambda x: '{x:.1f}'.format(x=x)})

        biggie.to_html(columns=['B', 'A'], float_format=str)
        biggie.to_html(columns=['B', 'A'], col_space=12, float_format=str)

        frame = DataFrame(index=np.arange(200))
        frame.to_html()

    def test_to_html_filename(self):
        biggie = DataFrame({'A': np.random.randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))

        biggie.loc[:20, 'A'] = np.nan
        biggie.loc[:20, 'B'] = np.nan
        with tm.ensure_clean('test.html') as path:
            biggie.to_html(path)
            with open(path, 'r') as f:
                s = biggie.to_html()
                s2 = f.read()
                assert s == s2

        frame = DataFrame(index=np.arange(200))
        with tm.ensure_clean('test.html') as path:
            frame.to_html(path)
            with open(path, 'r') as f:
                assert frame.to_html() == f.read()

    def test_to_html_with_no_bold(self):
        x = DataFrame({'x': np.random.randn(5)})
        ashtml = x.to_html(bold_rows=False)
        assert '<strong' not in ashtml[ashtml.find("</thead>")]

    def test_to_html_columns_arg(self):
        frame = DataFrame(tm.getSeriesData())
        result = frame.to_html(columns=['A'])
        assert '<th>B</th>' not in result

    def test_to_html_multiindex(self, datapath):
        columns = MultiIndex.from_tuples(list(zip(np.arange(2).repeat(2),
                                                  np.mod(lrange(4), 2))),
                                         names=['CL0', 'CL1'])
        df = DataFrame([list('abcd'), list('efgh')], columns=columns)
        result = df.to_html(justify='left')
        expected = expected_html(datapath, 'multiindex_1')
        assert result == expected

        columns = MultiIndex.from_tuples(list(zip(
            range(4), np.mod(
                lrange(4), 2))))
        df = DataFrame([list('abcd'), list('efgh')], columns=columns)

        result = df.to_html(justify='right')
        expected = expected_html(datapath, 'multiindex_2')
        assert result == expected

    @pytest.mark.parametrize("justify", fmt._VALID_JUSTIFY_PARAMETERS)
    def test_to_html_justify(self, justify, datapath):
        df = DataFrame({'A': [6, 30000, 2],
                        'B': [1, 2, 70000],
                        'C': [223442, 0, 1]},
                       columns=['A', 'B', 'C'])
        result = df.to_html(justify=justify)
        expected = expected_html(datapath, 'justify').format(justify=justify)
        assert result == expected

    @pytest.mark.parametrize("justify", ["super-right", "small-left",
                                         "noinherit", "tiny", "pandas"])
    def test_to_html_invalid_justify(self, justify):
        # GH 17527
        df = DataFrame()
        msg = "Invalid value for justify parameter"

        with pytest.raises(ValueError, match=msg):
            df.to_html(justify=justify)

    def test_to_html_index(self, datapath):
        index = ['foo', 'bar', 'baz']
        df = DataFrame({'A': [1, 2, 3],
                        'B': [1.2, 3.4, 5.6],
                        'C': ['one', 'two', np.nan]},
                       columns=['A', 'B', 'C'],
                       index=index)
        expected_with_index = expected_html(datapath, 'index_1')
        assert df.to_html() == expected_with_index

        expected_without_index = expected_html(datapath, 'index_2')
        result = df.to_html(index=False)
        for i in index:
            assert i not in result
        assert result == expected_without_index
        df.index = Index(['foo', 'bar', 'baz'], name='idx')
        expected_with_index = expected_html(datapath, 'index_3')
        assert df.to_html() == expected_with_index
        assert df.to_html(index=False) == expected_without_index

        tuples = [('foo', 'car'), ('foo', 'bike'), ('bar', 'car')]
        df.index = MultiIndex.from_tuples(tuples)

        expected_with_index = expected_html(datapath, 'index_4')
        assert df.to_html() == expected_with_index

        result = df.to_html(index=False)
        for i in ['foo', 'bar', 'car', 'bike']:
            assert i not in result
        # must be the same result as normal index
        assert result == expected_without_index

        df.index = MultiIndex.from_tuples(tuples, names=['idx1', 'idx2'])
        expected_with_index = expected_html(datapath, 'index_5')
        assert df.to_html() == expected_with_index
        assert df.to_html(index=False) == expected_without_index

    def test_to_html_with_classes(self, datapath):
        df = DataFrame()
        result = df.to_html(classes="sortable draggable")
        expected = expected_html(datapath, 'with_classes')
        assert result == expected

        result = df.to_html(classes=["sortable", "draggable"])
        assert result == expected

    def test_to_html_no_index_max_rows(self, datapath):
        # GH 14998
        df = DataFrame({"A": [1, 2, 3, 4]})
        result = df.to_html(index=False, max_rows=1)
        expected = expected_html(datapath, 'gh14998_expected_output')
        assert result == expected

    def test_to_html_multiindex_max_cols(self, datapath):
        # GH 6131
        index = MultiIndex(levels=[['ba', 'bb', 'bc'], ['ca', 'cb', 'cc']],
                           codes=[[0, 1, 2], [0, 1, 2]],
                           names=['b', 'c'])
        columns = MultiIndex(levels=[['d'], ['aa', 'ab', 'ac']],
                             codes=[[0, 0, 0], [0, 1, 2]],
                             names=[None, 'a'])
        data = np.array(
            [[1., np.nan, np.nan], [np.nan, 2., np.nan], [np.nan, np.nan, 3.]])
        df = DataFrame(data, index, columns)
        result = df.to_html(max_cols=2)
        expected = expected_html(datapath, 'gh6131_expected_output')
        assert result == expected

    @pytest.mark.parametrize('index', [False, 0])
    def test_to_html_truncation_index_false_max_rows(self, datapath, index):
        # GH 15019
        data = [[1.764052, 0.400157],
                [0.978738, 2.240893],
                [1.867558, -0.977278],
                [0.950088, -0.151357],
                [-0.103219, 0.410599]]
        df = DataFrame(data)
        result = df.to_html(max_rows=4, index=index)
        expected = expected_html(datapath, 'gh15019_expected_output')
        assert result == expected

    @pytest.mark.parametrize('index', [False, 0])
    def test_to_html_truncation_index_false_max_cols(self, datapath, index):
        # GH 22783
        data = [[1.764052, 0.400157, 0.978738, 2.240893, 1.867558],
                [-0.977278, 0.950088, -0.151357, -0.103219, 0.410599]]
        df = DataFrame(data)
        result = df.to_html(max_cols=4, index=index)
        expected = expected_html(datapath, 'gh22783_expected_output')
        assert result == expected

    def test_to_html_notebook_has_style(self):
        df = DataFrame({"A": [1, 2, 3]})
        result = df.to_html(notebook=True)
        assert "tbody tr th:only-of-type" in result
        assert "vertical-align: middle;" in result
        assert "thead th" in result

    def test_to_html_notebook_has_no_style(self):
        df = DataFrame({"A": [1, 2, 3]})
        result = df.to_html()
        assert "tbody tr th:only-of-type" not in result
        assert "vertical-align: middle;" not in result
        assert "thead th" not in result

    def test_to_html_with_index_names_false(self):
        # GH 16493
        df = DataFrame({"A": [1, 2]}, index=Index(['a', 'b'],
                                                  name='myindexname'))
        result = df.to_html(index_names=False)
        assert 'myindexname' not in result

    def test_to_html_with_id(self):
        # GH 8496
        df = DataFrame({"A": [1, 2]}, index=Index(['a', 'b'],
                                                  name='myindexname'))
        result = df.to_html(index_names=False, table_id="TEST_ID")
        assert ' id="TEST_ID"' in result

    def test_to_html_float_format_no_fixed_width(self, datapath):

        # GH 21625
        df = DataFrame({'x': [0.19999]})
        expected = expected_html(datapath, 'gh21625_expected_output')
        assert df.to_html(float_format='%.3f') == expected

        # GH 22270
        df = DataFrame({'x': [100.0]})
        expected = expected_html(datapath, 'gh22270_expected_output')
        assert df.to_html(float_format='%.0f') == expected

    @pytest.mark.parametrize("render_links, file_name", [
        (True, 'render_links_true'),
        (False, 'render_links_false'),
    ])
    def test_to_html_render_links(self, render_links, file_name, datapath):
        # GH 2679
        data = [
            [0, 'http://pandas.pydata.org/?q1=a&q2=b', 'pydata.org'],
            [0, 'www.pydata.org', 'pydata.org']
        ]
        df = DataFrame(data, columns=['foo', 'bar', None])

        result = df.to_html(render_links=render_links)
        expected = expected_html(datapath, file_name)
        assert result == expected
