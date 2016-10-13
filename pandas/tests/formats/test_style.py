import os
from nose import SkipTest

import copy
import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.util.testing import TestCase
import pandas.util.testing as tm

# Getting failures on a python 2.7 build with
# whenever we try to import jinja, whether it's installed or not.
# so we're explicitly skipping that one *before* we try to import
# jinja. We still need to export the imports as globals,
# since importing Styler tries to import jinja2.
job_name = os.environ.get('JOB_NAME', None)
if job_name == '27_slow_nnet_LOCALE':
    raise SkipTest("No jinja")
try:
    # Do try except on just jinja, so the only reason
    # We skip is if jinja can't import, not something else
    import jinja2  # noqa
except ImportError:
    raise SkipTest("No Jinja2")
from pandas.formats.style import Styler, _get_level_lengths  # noqa


class TestStyler(TestCase):

    def setUp(self):
        np.random.seed(24)
        self.s = DataFrame({'A': np.random.permutation(range(6))})
        self.df = DataFrame({'A': [0, 1], 'B': np.random.randn(2)})
        self.f = lambda x: x
        self.g = lambda x: x

        def h(x, foo='bar'):
            return pd.Series(['color: %s' % foo], index=x.index, name=x.name)

        self.h = h
        self.styler = Styler(self.df)
        self.attrs = pd.DataFrame({'A': ['color: red', 'color: blue']})
        self.dataframes = [
            self.df,
            pd.DataFrame({'f': [1., 2.], 'o': ['a', 'b'],
                          'c': pd.Categorical(['a', 'b'])})
        ]

    def test_init_non_pandas(self):
        with tm.assertRaises(TypeError):
            Styler([1, 2, 3])

    def test_init_series(self):
        result = Styler(pd.Series([1, 2]))
        self.assertEqual(result.data.ndim, 2)

    def test_repr_html_ok(self):
        self.styler._repr_html_()

    def test_update_ctx(self):
        self.styler._update_ctx(self.attrs)
        expected = {(0, 0): ['color: red'],
                    (1, 0): ['color: blue']}
        self.assertEqual(self.styler.ctx, expected)

    def test_update_ctx_flatten_multi(self):
        attrs = DataFrame({"A": ['color: red; foo: bar',
                                 'color: blue; foo: baz']})
        self.styler._update_ctx(attrs)
        expected = {(0, 0): ['color: red', ' foo: bar'],
                    (1, 0): ['color: blue', ' foo: baz']}
        self.assertEqual(self.styler.ctx, expected)

    def test_update_ctx_flatten_multi_traliing_semi(self):
        attrs = DataFrame({"A": ['color: red; foo: bar;',
                                 'color: blue; foo: baz;']})
        self.styler._update_ctx(attrs)
        expected = {(0, 0): ['color: red', ' foo: bar'],
                    (1, 0): ['color: blue', ' foo: baz']}
        self.assertEqual(self.styler.ctx, expected)

    def test_copy(self):
        s2 = copy.copy(self.styler)
        self.assertTrue(self.styler is not s2)
        self.assertTrue(self.styler.ctx is s2.ctx)  # shallow
        self.assertTrue(self.styler._todo is s2._todo)

        self.styler._update_ctx(self.attrs)
        self.styler.highlight_max()
        self.assertEqual(self.styler.ctx, s2.ctx)
        self.assertEqual(self.styler._todo, s2._todo)

    def test_deepcopy(self):
        s2 = copy.deepcopy(self.styler)
        self.assertTrue(self.styler is not s2)
        self.assertTrue(self.styler.ctx is not s2.ctx)
        self.assertTrue(self.styler._todo is not s2._todo)

        self.styler._update_ctx(self.attrs)
        self.styler.highlight_max()
        self.assertNotEqual(self.styler.ctx, s2.ctx)
        self.assertEqual(s2._todo, [])
        self.assertNotEqual(self.styler._todo, s2._todo)

    def test_clear(self):
        s = self.df.style.highlight_max()._compute()
        self.assertTrue(len(s.ctx) > 0)
        self.assertTrue(len(s._todo) > 0)
        s.clear()
        self.assertTrue(len(s.ctx) == 0)
        self.assertTrue(len(s._todo) == 0)

    def test_render(self):
        df = pd.DataFrame({"A": [0, 1]})
        style = lambda x: pd.Series(["color: red", "color: blue"], name=x.name)
        s = Styler(df, uuid='AB').apply(style)
        s.render()
        # it worked?

    def test_render_double(self):
        df = pd.DataFrame({"A": [0, 1]})
        style = lambda x: pd.Series(["color: red; border: 1px",
                                     "color: blue; border: 2px"], name=x.name)
        s = Styler(df, uuid='AB').apply(style)
        s.render()
        # it worked?

    def test_set_properties(self):
        df = pd.DataFrame({"A": [0, 1]})
        result = df.style.set_properties(color='white',
                                         size='10px')._compute().ctx
        # order is deterministic
        v = ["color: white", "size: 10px"]
        expected = {(0, 0): v, (1, 0): v}
        self.assertEqual(result.keys(), expected.keys())
        for v1, v2 in zip(result.values(), expected.values()):
            self.assertEqual(sorted(v1), sorted(v2))

    def test_set_properties_subset(self):
        df = pd.DataFrame({'A': [0, 1]})
        result = df.style.set_properties(subset=pd.IndexSlice[0, 'A'],
                                         color='white')._compute().ctx
        expected = {(0, 0): ['color: white']}
        self.assertEqual(result, expected)

    def test_empty_index_name_doesnt_display(self):
        # https://github.com/pandas-dev/pandas/pull/12090#issuecomment-180695902
        df = pd.DataFrame({'A': [1, 2], 'B': [3, 4], 'C': [5, 6]})
        result = df.style._translate()

        expected = [[{'class': 'blank level0', 'type': 'th', 'value': '',
                      'is_visible': True, 'display_value': ''},
                     {'class': 'col_heading level0 col0',
                      'display_value': 'A',
                      'type': 'th',
                      'value': 'A',
                      'is_visible': True,
                      'attributes': ["colspan=1"],
                      },
                     {'class': 'col_heading level0 col1',
                      'display_value': 'B',
                      'type': 'th',
                      'value': 'B',
                      'is_visible': True,
                      'attributes': ["colspan=1"],
                      },
                     {'class': 'col_heading level0 col2',
                      'display_value': 'C',
                      'type': 'th',
                      'value': 'C',
                      'is_visible': True,
                      'attributes': ["colspan=1"],
                      }]]

        self.assertEqual(result['head'], expected)

    def test_index_name(self):
        # https://github.com/pandas-dev/pandas/issues/11655
        df = pd.DataFrame({'A': [1, 2], 'B': [3, 4], 'C': [5, 6]})
        result = df.set_index('A').style._translate()

        expected = [[{'class': 'blank level0', 'type': 'th', 'value': '',
                      'display_value': '', 'is_visible': True},
                     {'class': 'col_heading level0 col0', 'type': 'th',
                      'value': 'B', 'display_value': 'B',
                      'is_visible': True, 'attributes': ['colspan=1']},
                     {'class': 'col_heading level0 col1', 'type': 'th',
                      'value': 'C', 'display_value': 'C',
                      'is_visible': True, 'attributes': ['colspan=1']}],
                    [{'class': 'index_name level0', 'type': 'th',
                      'value': 'A'},
                     {'class': 'blank', 'type': 'th', 'value': ''},
                     {'class': 'blank', 'type': 'th', 'value': ''}]]

        self.assertEqual(result['head'], expected)

    def test_multiindex_name(self):
        # https://github.com/pandas-dev/pandas/issues/11655
        df = pd.DataFrame({'A': [1, 2], 'B': [3, 4], 'C': [5, 6]})
        result = df.set_index(['A', 'B']).style._translate()

        expected = [[
            {'class': 'blank', 'type': 'th', 'value': '',
             'display_value': '', 'is_visible': True},
            {'class': 'blank level0', 'type': 'th', 'value': '',
             'display_value': '', 'is_visible': True},
            {'class': 'col_heading level0 col0', 'type': 'th',
             'value': 'C', 'display_value': 'C',
             'is_visible': True, 'attributes': ['colspan=1'],
             }],
            [{'class': 'index_name level0', 'type': 'th',
              'value': 'A'},
             {'class': 'index_name level1', 'type': 'th',
              'value': 'B'},
             {'class': 'blank', 'type': 'th', 'value': ''}]]

        self.assertEqual(result['head'], expected)

    def test_numeric_columns(self):
        # https://github.com/pandas-dev/pandas/issues/12125
        # smoke test for _translate
        df = pd.DataFrame({0: [1, 2, 3]})
        df.style._translate()

    def test_apply_axis(self):
        df = pd.DataFrame({'A': [0, 0], 'B': [1, 1]})
        f = lambda x: ['val: %s' % x.max() for v in x]
        result = df.style.apply(f, axis=1)
        self.assertEqual(len(result._todo), 1)
        self.assertEqual(len(result.ctx), 0)
        result._compute()
        expected = {(0, 0): ['val: 1'], (0, 1): ['val: 1'],
                    (1, 0): ['val: 1'], (1, 1): ['val: 1']}
        self.assertEqual(result.ctx, expected)

        result = df.style.apply(f, axis=0)
        expected = {(0, 0): ['val: 0'], (0, 1): ['val: 1'],
                    (1, 0): ['val: 0'], (1, 1): ['val: 1']}
        result._compute()
        self.assertEqual(result.ctx, expected)
        result = df.style.apply(f)  # default
        result._compute()
        self.assertEqual(result.ctx, expected)

    def test_apply_subset(self):
        axes = [0, 1]
        slices = [pd.IndexSlice[:], pd.IndexSlice[:, ['A']],
                  pd.IndexSlice[[1], :], pd.IndexSlice[[1], ['A']],
                  pd.IndexSlice[:2, ['A', 'B']]]
        for ax in axes:
            for slice_ in slices:
                result = self.df.style.apply(self.h, axis=ax, subset=slice_,
                                             foo='baz')._compute().ctx
                expected = dict(((r, c), ['color: baz'])
                                for r, row in enumerate(self.df.index)
                                for c, col in enumerate(self.df.columns)
                                if row in self.df.loc[slice_].index and
                                col in self.df.loc[slice_].columns)
                self.assertEqual(result, expected)

    def test_applymap_subset(self):
        def f(x):
            return 'foo: bar'

        slices = [pd.IndexSlice[:], pd.IndexSlice[:, ['A']],
                  pd.IndexSlice[[1], :], pd.IndexSlice[[1], ['A']],
                  pd.IndexSlice[:2, ['A', 'B']]]

        for slice_ in slices:
            result = self.df.style.applymap(f, subset=slice_)._compute().ctx
            expected = dict(((r, c), ['foo: bar'])
                            for r, row in enumerate(self.df.index)
                            for c, col in enumerate(self.df.columns)
                            if row in self.df.loc[slice_].index and
                            col in self.df.loc[slice_].columns)
            self.assertEqual(result, expected)

    def test_empty(self):
        df = pd.DataFrame({'A': [1, 0]})
        s = df.style
        s.ctx = {(0, 0): ['color: red'],
                 (1, 0): ['']}

        result = s._translate()['cellstyle']
        expected = [{'props': [['color', ' red']], 'selector': 'row0_col0'},
                    {'props': [['', '']], 'selector': 'row1_col0'}]
        self.assertEqual(result, expected)

    def test_bar(self):
        df = pd.DataFrame({'A': [0, 1, 2]})
        result = df.style.bar()._compute().ctx
        expected = {
            (0, 0): ['width: 10em', ' height: 80%'],
            (1, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,#d65f5f 50.0%, transparent 0%)'],
            (2, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,#d65f5f 100.0%, transparent 0%)']
        }
        self.assertEqual(result, expected)

        result = df.style.bar(color='red', width=50)._compute().ctx
        expected = {
            (0, 0): ['width: 10em', ' height: 80%'],
            (1, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,red 25.0%, transparent 0%)'],
            (2, 0): ['width: 10em', ' height: 80%',
                     'background: linear-gradient('
                     '90deg,red 50.0%, transparent 0%)']
        }
        self.assertEqual(result, expected)

        df['C'] = ['a'] * len(df)
        result = df.style.bar(color='red', width=50)._compute().ctx
        self.assertEqual(result, expected)
        df['C'] = df['C'].astype('category')
        result = df.style.bar(color='red', width=50)._compute().ctx
        self.assertEqual(result, expected)

    def test_bar_0points(self):
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = df.style.bar()._compute().ctx
        expected = {(0, 0): ['width: 10em', ' height: 80%'],
                    (0, 1): ['width: 10em', ' height: 80%'],
                    (0, 2): ['width: 10em', ' height: 80%'],
                    (1, 0): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 50.0%,'
                             ' transparent 0%)'],
                    (1, 1): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 50.0%,'
                             ' transparent 0%)'],
                    (1, 2): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 50.0%,'
                             ' transparent 0%)'],
                    (2, 0): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 100.0%'
                             ', transparent 0%)'],
                    (2, 1): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 100.0%'
                             ', transparent 0%)'],
                    (2, 2): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 100.0%'
                             ', transparent 0%)']}
        self.assertEqual(result, expected)

        result = df.style.bar(axis=1)._compute().ctx
        expected = {(0, 0): ['width: 10em', ' height: 80%'],
                    (0, 1): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 50.0%,'
                             ' transparent 0%)'],
                    (0, 2): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 100.0%'
                             ', transparent 0%)'],
                    (1, 0): ['width: 10em', ' height: 80%'],
                    (1, 1): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 50.0%'
                             ', transparent 0%)'],
                    (1, 2): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 100.0%'
                             ', transparent 0%)'],
                    (2, 0): ['width: 10em', ' height: 80%'],
                    (2, 1): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 50.0%'
                             ', transparent 0%)'],
                    (2, 2): ['width: 10em', ' height: 80%',
                             'background: linear-gradient(90deg,#d65f5f 100.0%'
                             ', transparent 0%)']}
        self.assertEqual(result, expected)

    def test_highlight_null(self, null_color='red'):
        df = pd.DataFrame({'A': [0, np.nan]})
        result = df.style.highlight_null()._compute().ctx
        expected = {(0, 0): [''],
                    (1, 0): ['background-color: red']}
        self.assertEqual(result, expected)

    def test_nonunique_raises(self):
        df = pd.DataFrame([[1, 2]], columns=['A', 'A'])
        with tm.assertRaises(ValueError):
            df.style

        with tm.assertRaises(ValueError):
            Styler(df)

    def test_caption(self):
        styler = Styler(self.df, caption='foo')
        result = styler.render()
        self.assertTrue(all(['caption' in result, 'foo' in result]))

        styler = self.df.style
        result = styler.set_caption('baz')
        self.assertTrue(styler is result)
        self.assertEqual(styler.caption, 'baz')

    def test_uuid(self):
        styler = Styler(self.df, uuid='abc123')
        result = styler.render()
        self.assertTrue('abc123' in result)

        styler = self.df.style
        result = styler.set_uuid('aaa')
        self.assertTrue(result is styler)
        self.assertEqual(result.uuid, 'aaa')

    def test_table_styles(self):
        style = [{'selector': 'th', 'props': [('foo', 'bar')]}]
        styler = Styler(self.df, table_styles=style)
        result = ' '.join(styler.render().split())
        self.assertTrue('th { foo: bar; }' in result)

        styler = self.df.style
        result = styler.set_table_styles(style)
        self.assertTrue(styler is result)
        self.assertEqual(styler.table_styles, style)

    def test_table_attributes(self):
        attributes = 'class="foo" data-bar'
        styler = Styler(self.df, table_attributes=attributes)
        result = styler.render()
        self.assertTrue('class="foo" data-bar' in result)

        result = self.df.style.set_table_attributes(attributes).render()
        self.assertTrue('class="foo" data-bar' in result)

    def test_precision(self):
        with pd.option_context('display.precision', 10):
            s = Styler(self.df)
        self.assertEqual(s.precision, 10)
        s = Styler(self.df, precision=2)
        self.assertEqual(s.precision, 2)

        s2 = s.set_precision(4)
        self.assertTrue(s is s2)
        self.assertEqual(s.precision, 4)

    def test_apply_none(self):
        def f(x):
            return pd.DataFrame(np.where(x == x.max(), 'color: red', ''),
                                index=x.index, columns=x.columns)
        result = (pd.DataFrame([[1, 2], [3, 4]])
                  .style.apply(f, axis=None)._compute().ctx)
        self.assertEqual(result[(1, 1)], ['color: red'])

    def test_trim(self):
        result = self.df.style.render()  # trim=True
        self.assertEqual(result.count('#'), 0)

        result = self.df.style.highlight_max().render()
        self.assertEqual(result.count('#'), len(self.df.columns))

    def test_highlight_max(self):
        df = pd.DataFrame([[1, 2], [3, 4]], columns=['A', 'B'])
        # max(df) = min(-df)
        for max_ in [True, False]:
            if max_:
                attr = 'highlight_max'
            else:
                df = -df
                attr = 'highlight_min'
            result = getattr(df.style, attr)()._compute().ctx
            self.assertEqual(result[(1, 1)], ['background-color: yellow'])

            result = getattr(df.style, attr)(color='green')._compute().ctx
            self.assertEqual(result[(1, 1)], ['background-color: green'])

            result = getattr(df.style, attr)(subset='A')._compute().ctx
            self.assertEqual(result[(1, 0)], ['background-color: yellow'])

            result = getattr(df.style, attr)(axis=0)._compute().ctx
            expected = {(1, 0): ['background-color: yellow'],
                        (1, 1): ['background-color: yellow'],
                        (0, 1): [''], (0, 0): ['']}
            self.assertEqual(result, expected)

            result = getattr(df.style, attr)(axis=1)._compute().ctx
            expected = {(0, 1): ['background-color: yellow'],
                        (1, 1): ['background-color: yellow'],
                        (0, 0): [''], (1, 0): ['']}
            self.assertEqual(result, expected)

        # separate since we cant negate the strs
        df['C'] = ['a', 'b']
        result = df.style.highlight_max()._compute().ctx
        expected = {(1, 1): ['background-color: yellow']}

        result = df.style.highlight_min()._compute().ctx
        expected = {(0, 0): ['background-color: yellow']}

    def test_export(self):
        f = lambda x: 'color: red' if x > 0 else 'color: blue'
        g = lambda x, y, z: 'color: %s' if x > 0 else 'color: %s' % z
        style1 = self.styler
        style1.applymap(f)\
            .applymap(g, y='a', z='b')\
            .highlight_max()
        result = style1.export()
        style2 = self.df.style
        style2.use(result)
        self.assertEqual(style1._todo, style2._todo)
        style2.render()

    def test_display_format(self):
        df = pd.DataFrame(np.random.random(size=(2, 2)))
        ctx = df.style.format("{:0.1f}")._translate()

        self.assertTrue(all(['display_value' in c for c in row]
                            for row in ctx['body']))
        self.assertTrue(all([len(c['display_value']) <= 3 for c in row[1:]]
                            for row in ctx['body']))
        self.assertTrue(
            len(ctx['body'][0][1]['display_value'].lstrip('-')) <= 3)

    def test_display_format_raises(self):
        df = pd.DataFrame(np.random.randn(2, 2))
        with tm.assertRaises(TypeError):
            df.style.format(5)
        with tm.assertRaises(TypeError):
            df.style.format(True)

    def test_display_subset(self):
        df = pd.DataFrame([[.1234, .1234], [1.1234, 1.1234]],
                          columns=['a', 'b'])
        ctx = df.style.format({"a": "{:0.1f}", "b": "{0:.2%}"},
                              subset=pd.IndexSlice[0, :])._translate()
        expected = '0.1'
        self.assertEqual(ctx['body'][0][1]['display_value'], expected)
        self.assertEqual(ctx['body'][1][1]['display_value'], '1.1234')
        self.assertEqual(ctx['body'][0][2]['display_value'], '12.34%')

        raw_11 = '1.1234'
        ctx = df.style.format("{:0.1f}",
                              subset=pd.IndexSlice[0, :])._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], expected)
        self.assertEqual(ctx['body'][1][1]['display_value'], raw_11)

        ctx = df.style.format("{:0.1f}",
                              subset=pd.IndexSlice[0, :])._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], expected)
        self.assertEqual(ctx['body'][1][1]['display_value'], raw_11)

        ctx = df.style.format("{:0.1f}",
                              subset=pd.IndexSlice['a'])._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], expected)
        self.assertEqual(ctx['body'][0][2]['display_value'], '0.1234')

        ctx = df.style.format("{:0.1f}",
                              subset=pd.IndexSlice[0, 'a'])._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], expected)
        self.assertEqual(ctx['body'][1][1]['display_value'], raw_11)

        ctx = df.style.format("{:0.1f}",
                              subset=pd.IndexSlice[[0, 1], ['a']])._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], expected)
        self.assertEqual(ctx['body'][1][1]['display_value'], '1.1')
        self.assertEqual(ctx['body'][0][2]['display_value'], '0.1234')
        self.assertEqual(ctx['body'][1][2]['display_value'], '1.1234')

    def test_display_dict(self):
        df = pd.DataFrame([[.1234, .1234], [1.1234, 1.1234]],
                          columns=['a', 'b'])
        ctx = df.style.format({"a": "{:0.1f}", "b": "{0:.2%}"})._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], '0.1')
        self.assertEqual(ctx['body'][0][2]['display_value'], '12.34%')
        df['c'] = ['aaa', 'bbb']
        ctx = df.style.format({"a": "{:0.1f}", "c": str.upper})._translate()
        self.assertEqual(ctx['body'][0][1]['display_value'], '0.1')
        self.assertEqual(ctx['body'][0][3]['display_value'], 'AAA')

    def test_bad_apply_shape(self):
        df = pd.DataFrame([[1, 2], [3, 4]])
        with tm.assertRaises(ValueError):
            df.style._apply(lambda x: 'x', subset=pd.IndexSlice[[0, 1], :])

        with tm.assertRaises(ValueError):
            df.style._apply(lambda x: [''], subset=pd.IndexSlice[[0, 1], :])

        with tm.assertRaises(ValueError):
            df.style._apply(lambda x: ['', '', '', ''])

        with tm.assertRaises(ValueError):
            df.style._apply(lambda x: ['', '', ''], subset=1)

        with tm.assertRaises(ValueError):
            df.style._apply(lambda x: ['', '', ''], axis=1)

    def test_apply_bad_return(self):
        def f(x):
            return ''
        df = pd.DataFrame([[1, 2], [3, 4]])
        with tm.assertRaises(TypeError):
            df.style._apply(f, axis=None)

    def test_apply_bad_labels(self):
        def f(x):
            return pd.DataFrame(index=[1, 2], columns=['a', 'b'])
        df = pd.DataFrame([[1, 2], [3, 4]])
        with tm.assertRaises(ValueError):
            df.style._apply(f, axis=None)

    def test_get_level_lengths(self):
        index = pd.MultiIndex.from_product([['a', 'b'], [0, 1, 2]])
        expected = {(0, 0): 3, (0, 3): 3, (1, 0): 1, (1, 1): 1, (1, 2): 1,
                    (1, 3): 1, (1, 4): 1, (1, 5): 1}
        result = _get_level_lengths(index)
        tm.assert_dict_equal(result, expected)

    def test_get_level_lengths_un_sorted(self):
        index = pd.MultiIndex.from_arrays([
            [1, 1, 2, 1],
            ['a', 'b', 'b', 'd']
        ])
        expected = {(0, 0): 2, (0, 2): 1, (0, 3): 1,
                    (1, 0): 1, (1, 1): 1, (1, 2): 1, (1, 3): 1}
        result = _get_level_lengths(index)
        tm.assert_dict_equal(result, expected)

    def test_mi_sparse(self):
        df = pd.DataFrame({'A': [1, 2]},
                          index=pd.MultiIndex.from_arrays([['a', 'a'],
                                                           [0, 1]]))
        result = df.style._translate()
        body_0 = result['body'][0][0]
        expected_0 = {
            "value": "a", "display_value": "a", "is_visible": True,
            "type": "th", "attributes": ["rowspan=2"],
            "class": "row_heading level0 row0",
        }
        tm.assert_dict_equal(body_0, expected_0)

        body_1 = result['body'][0][1]
        expected_1 = {
            "value": 0, "display_value": 0, "is_visible": True,
            "type": "th", "attributes": ["rowspan=1"],
            "class": "row_heading level1 row0",
        }
        tm.assert_dict_equal(body_1, expected_1)

        body_10 = result['body'][1][0]
        expected_10 = {
            "value": 'a', "display_value": 'a', "is_visible": False,
            "type": "th", "attributes": ["rowspan=1"],
            "class": "row_heading level0 row1",
        }
        tm.assert_dict_equal(body_10, expected_10)

        head = result['head'][0]
        expected = [
            {'type': 'th', 'class': 'blank', 'value': '',
             'is_visible': True, "display_value": ''},
            {'type': 'th', 'class': 'blank level0', 'value': '',
             'is_visible': True, 'display_value': ''},
            {'attributes': ['colspan=1'], 'class': 'col_heading level0 col0',
             'is_visible': True, 'type': 'th', 'value': 'A',
             'display_value': 'A'}]
        self.assertEqual(head, expected)

    def test_mi_sparse_disabled(self):
        with pd.option_context('display.multi_sparse', False):
            df = pd.DataFrame({'A': [1, 2]},
                              index=pd.MultiIndex.from_arrays([['a', 'a'],
                                                              [0, 1]]))
            result = df.style._translate()
        body = result['body']
        for row in body:
            self.assertEqual(row[0]['attributes'], ['rowspan=1'])

    def test_mi_sparse_index_names(self):
        df = pd.DataFrame({'A': [1, 2]}, index=pd.MultiIndex.from_arrays(
            [['a', 'a'], [0, 1]],
            names=['idx_level_0', 'idx_level_1'])
        )
        result = df.style._translate()
        head = result['head'][1]
        expected = [{
            'class': 'index_name level0', 'value': 'idx_level_0',
            'type': 'th'},
            {'class': 'index_name level1', 'value': 'idx_level_1',
             'type': 'th'},
            {'class': 'blank', 'value': '', 'type': 'th'}]

        self.assertEqual(head, expected)

    def test_mi_sparse_column_names(self):
        df = pd.DataFrame(
            np.arange(16).reshape(4, 4),
            index=pd.MultiIndex.from_arrays(
                [['a', 'a', 'b', 'a'], [0, 1, 1, 2]],
                names=['idx_level_0', 'idx_level_1']),
            columns=pd.MultiIndex.from_arrays(
                [['C1', 'C1', 'C2', 'C2'], [1, 0, 1, 0]],
                names=['col_0', 'col_1']
            )
        )
        result = df.style._translate()
        head = result['head'][1]
        expected = [
            {'class': 'blank', 'value': '', 'display_value': '',
             'type': 'th', 'is_visible': True},
            {'class': 'index_name level1', 'value': 'col_1',
             'display_value': 'col_1', 'is_visible': True, 'type': 'th'},
            {'attributes': ['colspan=1'],
             'class': 'col_heading level1 col0',
             'display_value': 1,
             'is_visible': True,
             'type': 'th',
             'value': 1},
            {'attributes': ['colspan=1'],
             'class': 'col_heading level1 col1',
             'display_value': 0,
             'is_visible': True,
             'type': 'th',
             'value': 0},

            {'attributes': ['colspan=1'],
             'class': 'col_heading level1 col2',
             'display_value': 1,
             'is_visible': True,
             'type': 'th',
             'value': 1},

            {'attributes': ['colspan=1'],
             'class': 'col_heading level1 col3',
             'display_value': 0,
             'is_visible': True,
             'type': 'th',
             'value': 0},
        ]
        self.assertEqual(head, expected)


@tm.mplskip
class TestStylerMatplotlibDep(TestCase):

    def test_background_gradient(self):
        df = pd.DataFrame([[1, 2], [2, 4]], columns=['A', 'B'])
        for axis in [0, 1, 'index', 'columns']:
            for cmap in [None, 'YlOrRd']:
                result = df.style.background_gradient(cmap=cmap)._compute().ctx
                self.assertTrue(all("#" in x[0] for x in result.values()))
                self.assertEqual(result[(0, 0)], result[(0, 1)])
                self.assertEqual(result[(1, 0)], result[(1, 1)])

        result = (df.style.background_gradient(subset=pd.IndexSlice[1, 'A'])
                    ._compute().ctx)
        self.assertEqual(result[(1, 0)], ['background-color: #fff7fb'])
